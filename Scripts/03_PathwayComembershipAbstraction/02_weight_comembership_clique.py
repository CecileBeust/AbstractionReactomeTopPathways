from SPARQLWrapper import SPARQLWrapper, JSON
import subprocess
import time
import os
import requests
import rdflib
import pandas as pd
import numpy as np
import networkx as nx
import itertools
import glob
from utility import *

current_directory = os.getcwd()
BioPAX_Ontology_file_path = os.path.join(current_directory, 'Data/BioPAX/BioPAXOntology/biopax-level3.owl')
ReactomeBioPAX_file_path = os.path.join(current_directory, 'Data/BioPAX/ReactomeTopPathways')
results_dir = os.path.join(current_directory, 'Results/PathwayComembership')

os.makedirs(results_dir, exist_ok=True)
os.makedirs(os.path.join(results_dir, "ScoresPathways"), exist_ok=True)
os.makedirs(os.path.join(results_dir, "WeightedComembershipCliques"), exist_ok=True)

prefixes = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX bp3: <http://www.biopax.org/release/biopax-level3.owl#>
"""

query = """
SELECT (COUNT(DISTINCT ?id) AS ?nbID)
WHERE {
    VALUES ?db { "UniProt" "UniProt Isoform" }
    ?entityRef rdf:type/(rdfs:subClassOf*) bp3:ProteinReference .
    ?entityRef bp3:xref ?entityRefXref .
    ?entityRefXref bp3:db ?db .
    ?entityRefXref bp3:id ?id .
}
"""

endpoint = "http://localhost:3030/top_pathway"
counter = 1
#filelist = glob.glob(os.path.join(ReactomeBioPAX_file_path, '*.xml'))
filelist = glob.glob(os.path.join(ReactomeBioPAX_file_path, '23_ProgrammedCellDeath.xml'))

def wait_for_fuseki(base_endpoint, timeout=300, poll_interval=1):
    """Poll Fuseki until it's ready instead of a fixed sleep."""
    ping_url = base_endpoint.rsplit('/', 1)[0] + "/$/ping"
    start = time.time()
    while time.time() - start < timeout:
        try:
            r = requests.get(ping_url, timeout=2)
            if r.status_code == 200:
                return
        except requests.exceptions.RequestException:
            pass
        time.sleep(poll_interval)
    raise TimeoutError(f"Fuseki did not started in ({timeout}s)")


for owl_file in sorted(filelist):
    start_time = time.time()

    # ------------------------------------------------------------------
    # 0 - Load pathway abstraction
    # ------------------------------------------------------------------
    pathway_abstraction = pd.read_csv(
        f"Results/PathwayAbstraction/{counter:02d}_WeightedPathwayAbstraction.csv",
        sep=",", header=0
    )
    G = create_networkx_graph(pathway_abstraction)
    hierarchy_graph = extract_subgraph_by_interaction(G, "abstraction:IsAComponentOf")

    # ------------------------------------------------------------------
    # 1 - get number of Uniprot IDs per pathway
    # ------------------------------------------------------------------
    command = [
        '/home/cbeust/Softwares/JenaFuseki/apache-jena-fuseki-4.9.0/fuseki-server',
        '--file', owl_file,
        '--file', BioPAX_Ontology_file_path,
        '/top_pathway'
    ]
    print("Fuseki command:", command)
    process = subprocess.Popen(command)

    try:
        wait_for_fuseki(endpoint)
        sparql = SPARQLWrapper(endpoint)
        sparql.setQuery(prefixes + query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        nbEr_total = int(results["results"]["bindings"][0]["nbID"]["value"])

        # ------------------------------------------------------------------
        # 2 - get list of Uniprot IDs per pathway
        # ------------------------------------------------------------------
        er_per_pathway = pd.read_csv(
            f"Results/UtilityFiles/{counter:02d}_UpPerPathway.csv",
            sep=",", header=0
        )

        dico_er_per_pathway = {}
        for pathway, id_ in er_per_pathway.iloc[:, [0, 1]].values:
            dico_er_per_pathway.setdefault(pathway, []).append(id_)

        # ------------------------------------------------------------------
        # 3 - compute Score IC (with content of UniProt IDs) of pathways
        # ------------------------------------------------------------------
        dict_score_ic_pathways = {
            pathway: compute_ic_er_pathway(pathway, dico_er_per_pathway, nbEr_total)
            for pathway in dico_er_per_pathway
        }

        score_IC_pathways = pd.DataFrame({
            'pathway': list(dict_score_ic_pathways.keys()),
            'score': list(dict_score_ic_pathways.values())
        })
        score_IC_pathways.to_csv(
            f"Results/PathwayComembership/ScoresPathways/{counter:02d}_ScoreICPathways.csv",
            sep=",", header=True, index=False
        )

        # ------------------------------------------------------------------
        # 4 - get next step pathway scores (dict lookup instead of DataFrame .loc)
        # ------------------------------------------------------------------
        list_pathways = []
        seen_pathways = set()
        for pathway1, _, pathway2 in pathway_abstraction.iloc[:, [0, 1, 2]].values:
            if pathway1 not in seen_pathways:
                seen_pathways.add(pathway1)
                list_pathways.append(pathway1)
            if pathway2 not in seen_pathways:
                seen_pathways.add(pathway2)
                list_pathways.append(pathway2)

        mat_lookup = {}
        for pathway1, interaction, pathway2, weight_er in pathway_abstraction.iloc[:, [0, 1, 2, 3]].values:
            if interaction == "abstraction:NextStepPathway":
                mat_lookup[(pathway1, pathway2)] = weight_er
                mat_lookup[(pathway2, pathway1)] = weight_er

        # ------------------------------------------------------------------
        # 5 - get pathway ancestors
        # ------------------------------------------------------------------
        dico_ancestors = {}
        for pathway in list_pathways:
            anc, _ = find_ancestors_and_descendants(hierarchy_graph, pathway)
            dico_ancestors[pathway] = anc

        up_per_pathway = pd.read_csv(
            f"Results/UtilityFiles/{counter:02d}_UpPerPathway.csv",
            sep=",", header=0
        )

        dico_parents_up = {}
        for pathway_val, prot_val in up_per_pathway.iloc[:, [0, 1]].values:
            dico_parents_up.setdefault(prot_val, []).append(pathway_val)

        root = find_root(G, "abstraction:IsAComponentOf")

        # precompute shortest paths
        reversed_hierarchy = hierarchy_graph.reverse()
        shortest_path_len_from_root = {}
        for pathway in list_pathways:
            try:
                shortest_path_len_from_root[pathway] = len(
                    nx.shortest_path(reversed_hierarchy, root, pathway)
                )
            except nx.NetworkXNoPath:
                shortest_path_len_from_root[pathway] = 0

        dico_most_precise_parent = {}
        print(len(dico_parents_up))
        for prot, parents in dico_parents_up.items():
            if len(parents) == 1:
                dico_most_precise_parent[prot] = parents
            else:
                shortest_paths = {
                    parent: shortest_path_len_from_root.get(parent, 0)
                    for parent in parents
                }
                max_len = max(shortest_paths.values())
                dico_most_precise_parent[prot] = [
                    key for key, val in shortest_paths.items() if val == max_len
                ]

        # ------------------------------------------------------------------
        # 6 - weight comembership clique
        # ------------------------------------------------------------------
        comembership_clique = pd.read_csv(
            f"Results/PathwayComembership/ComembershipCliques/{counter:02d}_ComembershipClique.csv",
            sep=",", header=0
        )
        ancestors_cache = {}
        mica_cache = {}
        next_step_cache = {}

        def get_ancestors_for_entity(entity):
            cached = ancestors_cache.get(entity)
            if cached is None:
                parents = dico_most_precise_parent[entity]
                ancestors_pathways = get_pathways_and_ancestor(parents, dico_ancestors)
                cached = (parents, ancestors_pathways)
                ancestors_cache[entity] = cached
            return cached

        def get_mica(parents1, parents2):
            key = (tuple(sorted(parents1)), tuple(sorted(parents2)))
            cached = mica_cache.get(key)
            if cached is None:
                cached = find_mica_from_several_pathways(
                    hierarchy_graph, root, parents1, parents2,
                    dico_er_per_pathway, nbEr_total
                )
                mica_cache[key] = cached
            return cached

        def get_max_next_step_score(ancestors1, ancestors2):
            key = (tuple(sorted(ancestors1)), tuple(sorted(ancestors2)))
            cached = next_step_cache.get(key)
            if cached is None:
                best = 0
                for p1, p2 in itertools.product(ancestors1, ancestors2):
                    if p1 != p2:
                        val = mat_lookup.get((p1, p2))
                        if val and val > best:
                            best = val
                cached = best
                next_step_cache[key] = cached
            return cached

        rows = []
        for entity1, entity2 in comembership_clique.iloc[:, [0, 1]].values:
            if entity1 in dico_most_precise_parent and entity2 in dico_most_precise_parent:
                parents1, ancestors_pathways1 = get_ancestors_for_entity(entity1)
                parents2, ancestors_pathways2 = get_ancestors_for_entity(entity2)

                mica = get_mica(parents1, parents2)
                score_IC_UP_mica = dict_score_ic_pathways[mica]

                max_next_step_pathway_score = get_max_next_step_score(
                    ancestors_pathways1, ancestors_pathways2
                )

                score_comembership = max(score_IC_UP_mica, max_next_step_pathway_score)

                if score_comembership == score_IC_UP_mica:
                    rows.append([entity1, entity2, score_comembership, "scoreMICA"])
                else:
                    rows.append([entity1, entity2, score_comembership, "scoreNextStep"])

        weighted_comembership_clique = pd.DataFrame(
            rows, columns=["entity1", "entity2", "scoreComembership", "provenanceScore"]
        )
        weighted_comembership_clique = weighted_comembership_clique.sort_values(
            by="scoreComembership", ascending=False
        )
        weighted_comembership_clique.to_csv(
            f"Results/PathwayComembership/WeightedComembershipCliques/{counter:02d}_WeightedComembershipClique.csv",
            sep=",", header=True, index=False
        )

        elapsed_time = time.time() - start_time
        print(f"File {owl_file} processed in {elapsed_time:.2f}s")

    finally:
        process.kill()
        process.wait()

    counter += 1