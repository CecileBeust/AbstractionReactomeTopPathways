"""
weight_comembership_clique.py

For each Reactome top pathway BioPAX file:
    1. Launch a temporary Jena Fuseki SPARQL endpoint over the file
    2. Count total UniProt entity references = UniProt IDs and how many belong to each
        pathway
    3. Compute an Information Content (IC) score per pathway based on the number of UniProt IDs of the pathway
    4. Weight every pair of co-membered proteins (from a precomputed
        comembership clique) using the max of :
        - the IC score of their Most Informative Common Ancestor (MICA)
        - the highest NextStepPathway score connecting their ancestors
    5. Write the weighted comembership clique to CSV

    Expected inputs (per pathway file, indexed by a `counter`):
    - Results/PathwayAbstraction/01_TopPathways/{counter:02d}_WeightedPathwayAbstraction.csv
    - Results/UtilityFiles/{counter:02d}_UpPerPathway.csv
    - Results/PathwayComembership/01_ComembershipCliques/{counter:02d}_ComembershipClique.csv

    Outputs (per pathway file):
    - Results/PathwayComembership/04_ScoresPathways/{counter:02d}_ScoreICPathways.csv
    - Results/PathwayComembership/02_WeightedComembershipCliques/{counter:02d}_WeightedComembershipClique.csv
"""

# ----------------------------------------------------------------------------
# Imports
# ----------------------------------------------------------------------------
from __future__ import annotations
import glob
import csv
import itertools
import os
import subprocess
import time
import networkx as nx
import pandas as pd
import requests
from SPARQLWrapper import JSON, SPARQLWrapper
from utility import (
    compute_ic_er_pathway,
    create_networkx_graph,
    extract_subgraph_by_interaction,
    find_ancestors_and_descendants,
    find_mica_from_several_pathways,
    find_root,
    get_pathways_and_ancestor,
)


# ---------------------------------------------------------------------------
# File Paths and configurations
# ---------------------------------------------------------------------------
CURRENT_DIR = os.getcwd()
BIOPAX_ONTOLOGY_PATH = os.path.join(CURRENT_DIR, "Data/BioPAX/BioPAXOntology/biopax-level3.owl")
REACTOME_TOP_PATHWAYS_DIR = os.path.join(CURRENT_DIR, "Data/BioPAX/ReactomeTopPathways")
REACTOME_BIOPAX_FILE = os.path.join(CURRENT_DIR, "Data/BioPAX/ReactomeBioPAX/Homo_sapiens_v96.owl")
RESULTS_DIR = os.path.join(CURRENT_DIR, "Results/PathwayComembership")
SCORES_DIR = os.path.join(RESULTS_DIR, "04_ScoresPathways")
WEIGHTED_CLIQUES_DIR = os.path.join(RESULTS_DIR, "02_WeightedComembershipCliques")
MATRICES_TOP_PATHWAYS = os.path.join(RESULTS_DIR, "03_MatricesTopPathways")
MATRICES_REACTOME = os.path.join(RESULTS_DIR, "05_MatricesReactome")

FUSEKI_BINARY = "/home/cbeust/Softwares/JenaFuseki/apache-jena-fuseki-4.9.0/fuseki-server"
FUSEKI_ENDPOINT_TOP_PATHWAY = "http://localhost:3030/top_pathway"
FUSEKI_ENDPOINT_REACTOME = "http://localhost:3030/reactome"
FUSEKI_DATASET_TOP_PATHWAYS = "/top_pathway"
FUSEKI_DATASET_REACTOME = "/reactome"

SPARQL_PREFIXES = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX bp3: <http://www.biopax.org/release/biopax-level3.owl#>
"""

COUNT_UNIPROT_IDS_QUERY = """
SELECT (COUNT(DISTINCT ?id) AS ?nbID)
WHERE {
    VALUES ?db { "UniProt" "UniProt Isoform" }
    ?entityRef rdf:type/(rdfs:subClassOf*) bp3:ProteinReference .
    ?entityRef bp3:xref ?entityRefXref .
    ?entityRefXref bp3:db ?db .
    ?entityRefXref bp3:id ?id .
}
"""

HIERARCHY_INTERACTION = "abstraction:IsAComponentOf"
NEXT_STEP_INTERACTION = "abstraction:NextStepPathway"


# ---------------------------------------------------------------------------
# Fuseki server lifecycle
# ---------------------------------------------------------------------------

def wait_for_fuseki(base_endpoint: str, timeout: int = 300, poll_interval: float = 1) -> None:
    """Poll Fuseki's ping endpoint until it responds, instead of a fixed sleep."""
    ping_url = base_endpoint.rsplit("/", 1)[0] + "/$/ping"
    start = time.time()
    while time.time() - start < timeout:
        try:
            response = requests.get(ping_url, timeout=2)
            if response.status_code == 200:
                return
        except requests.exceptions.RequestException:
            pass
        time.sleep(poll_interval)
    raise TimeoutError(f"Fuseki did not start in ({timeout}s)")


def start_fuseki_server(owl_file: str, FUSEKI_DATASET) -> subprocess.Popen:
    """Launch a Fuseki server loaded with the pathway file and the BioPAX ontology."""
    command = [
        FUSEKI_BINARY,
        "--file", owl_file,
        "--file", BIOPAX_ONTOLOGY_PATH,
        FUSEKI_DATASET,
    ]
    print("Fuseki command:", command)
    return subprocess.Popen(command)


def count_total_uniprot_ids(endpoint: str) -> int:
    """Run the SPARQL query counting distinct UniProt IDs across the whole file."""
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(SPARQL_PREFIXES + COUNT_UNIPROT_IDS_QUERY)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    return int(results["results"]["bindings"][0]["nbID"]["value"])


# ---------------------------------------------------------------------------
# Per-file input loading
# ---------------------------------------------------------------------------

def load_pathway_abstraction(counter: int) -> pd.DataFrame:
    """Load the weighted pathway-abstraction edge table for this pathway file."""
    return pd.read_csv(
        f"Results/PathwayAbstraction/01_TopPathways/{counter:02d}_WeightedPathwayAbstraction.csv",
        sep=",", header=0,
    )


def load_entity_refs_per_pathway(counter: int) -> pd.DataFrame:
    """Load the (pathway, UniProt ID) table for this pathway file."""
    return pd.read_csv(
        f"Results/UtilityFiles/{counter:02d}_UpPerPathway.csv",
        sep=",", header=0,
    )


def load_comembership_clique(counter: int) -> pd.DataFrame:
    """Load the precomputed entity co-membership clique table for this pathway file."""
    return pd.read_csv(
        f"Results/PathwayComembership/01_ComembershipCliques/{counter:02d}_ComembershipClique.csv",
        sep=",", header=0,
    )


def build_entity_refs_per_pathway_dict(er_per_pathway: pd.DataFrame) -> dict[str, list[str]]:
    """Build {pathway: [uniprot_id, ...]} from the raw (pathway, id) table."""
    dico: dict[str, list[str]] = {}
    for pathway, entity_id in er_per_pathway.iloc[:, [0, 1]].values:
        dico.setdefault(pathway, []).append(entity_id)
    return dico


# ---------------------------------------------------------------------------
# Functions to compute IC scores per pathway based on the number of UniProt IDs
# ---------------------------------------------------------------------------

def compute_ic_scores(dico_er_per_pathway: dict, nb_er_total: int) -> dict[str, float]:
    """Compute the Information Content score of every pathway."""
    return {
        pathway: compute_ic_er_pathway(pathway, dico_er_per_pathway, nb_er_total)
        for pathway in dico_er_per_pathway
    }


def save_ic_scores(dict_score_ic_pathways: dict[str, float], counter: int) -> None:
    """Write the per-pathway IC scores to CSV."""
    score_ic_pathways = pd.DataFrame({
        "pathway": list(dict_score_ic_pathways.keys()),
        "score": list(dict_score_ic_pathways.values()),
    })
    score_ic_pathways.to_csv(
        os.path.join(SCORES_DIR, f"{counter:02d}_ScoreICPathways.csv"),
        sep=",", header=True, index=False,
    )


# ---------------------------------------------------------------------------
# Pathway hierarchy helpers
# ---------------------------------------------------------------------------

def get_ordered_pathway_list(pathway_abstraction: pd.DataFrame) -> list[str]:
    """Return every pathway appearing as source or target, in first-seen order."""
    ordered_pathways: list[str] = []
    seen: set[str] = set()
    for pathway1, _, pathway2 in pathway_abstraction.iloc[:, [0, 1, 2]].values:
        for pathway in (pathway1, pathway2):
            if pathway not in seen:
                seen.add(pathway)
                ordered_pathways.append(pathway)
    return ordered_pathways


def build_next_step_lookup(pathway_abstraction: pd.DataFrame) -> dict[tuple[str, str], float]:
    """Build a symmetric {(pathway1, pathway2): weight} lookup for NextStepPathway edges."""
    lookup: dict[tuple[str, str], float] = {}
    for pathway1, interaction, pathway2, weight_er in pathway_abstraction.iloc[:, [0, 1, 2, 3]].values:
        if interaction == NEXT_STEP_INTERACTION:
            lookup[(pathway1, pathway2)] = weight_er
            lookup[(pathway2, pathway1)] = weight_er
    return lookup


def build_ancestors_dict(list_pathways: list[str], hierarchy_graph: nx.DiGraph) -> dict[str, set]:
    """Precompute {pathway: ancestors} for every pathway in the hierarchy."""
    dico_ancestors: dict[str, set] = {}
    for pathway in list_pathways:
        ancestors, _ = find_ancestors_and_descendants(hierarchy_graph, pathway)
        dico_ancestors[pathway] = ancestors
    return dico_ancestors


def build_protein_to_pathways_dict(up_per_pathway: pd.DataFrame) -> dict[str, list[str]]:
    """Build {protein: [pathway, ...]} from the (pathway, protein) table."""
    dico_parents_up: dict[str, list[str]] = {}
    for pathway_val, prot_val in up_per_pathway.iloc[:, [0, 1]].values:
        dico_parents_up.setdefault(prot_val, []).append(pathway_val)
    return dico_parents_up


def compute_shortest_path_lengths_from_root(
    hierarchy_graph: nx.DiGraph, root: str, list_pathways: list[str]
) -> dict[str, int]:
    """
    Compute, for every pathway, its shortest-path distance from the
    hierarchy root (used to determine which parent pathway is "most
    precise" i.e. deepest, when a protein belongs to several pathways).
    """
    reversed_hierarchy = hierarchy_graph.reverse()
    lengths: dict[str, int] = {}
    for pathway in list_pathways:
        try:
            lengths[pathway] = len(nx.shortest_path(reversed_hierarchy, root, pathway))
        except nx.NetworkXNoPath:
            lengths[pathway] = 0
    return lengths


def build_most_precise_parent_map(
    dico_parents_up: dict[str, list[str]],
    shortest_path_len_from_root: dict[str, int],
) -> dict[str, list[str]]:
    """
    For each protein, keep only its "most precise" parent pathway(s):
    the pathway(s) furthest from the root. Ties are kept as a list.
    """
    dico_most_precise_parent: dict[str, list[str]] = {}
    for protein, parents in dico_parents_up.items():
        if len(parents) == 1:
            dico_most_precise_parent[protein] = parents
            continue

        depths = {parent: shortest_path_len_from_root.get(parent, 0) for parent in parents}
        max_depth = max(depths.values())
        dico_most_precise_parent[protein] = [
            parent for parent, depth in depths.items() if depth == max_depth
        ]
    return dico_most_precise_parent

# ---------------------------------------------------------------------------
# Get number of total Uniprot IDs in Reactome
# ---------------------------------------------------------------------------
fuseki_process_reactome = start_fuseki_server(REACTOME_BIOPAX_FILE, FUSEKI_DATASET_REACTOME)
try:
    wait_for_fuseki(FUSEKI_ENDPOINT_REACTOME)
    nb_er_reactome = count_total_uniprot_ids(FUSEKI_ENDPOINT_REACTOME)
    print(nb_er_reactome)
finally:
    fuseki_process_reactome.kill()
    fuseki_process_reactome.wait()


# ---------------------------------------------------------------------------
# Co-membership clique scoring
# ---------------------------------------------------------------------------

class ComembershipScorer:
    """
    Scores pairs of co-membered entities, caching intermediate results
    (per-entity ancestor lookups, per-parent-set MICA, and per-ancestor-set
    best next-step score) since the same entities/parent sets recur across
    many pairs in the clique table.
    """

    def __init__(
        self,
        hierarchy_graph: nx.DiGraph,
        root: str,
        dico_most_precise_parent: dict[str, list[str]],
        dico_ancestors: dict[str, set],
        dico_er_per_pathway: dict[str, list[str]],
        nb_er_total: int,
        next_step_lookup: dict[tuple[str, str], float],
    ) -> None:
        self.hierarchy_graph = hierarchy_graph
        self.root = root
        self.dico_most_precise_parent = dico_most_precise_parent
        self.dico_ancestors = dico_ancestors
        self.dico_er_per_pathway = dico_er_per_pathway
        self.nb_er_total = nb_er_total
        self.next_step_lookup = next_step_lookup

        self._ancestors_cache: dict[str, tuple[list[str], set]] = {}
        self._mica_cache: dict[tuple, str] = {}
        self._next_step_cache: dict[tuple, float] = {}

    def _get_parents_and_ancestor_pathways(self, entity: str) -> tuple[list[str], set]:
        if entity not in self._ancestors_cache:
            parents = self.dico_most_precise_parent[entity]
            ancestor_pathways = get_pathways_and_ancestor(parents, self.dico_ancestors)
            self._ancestors_cache[entity] = (parents, ancestor_pathways)
        return self._ancestors_cache[entity]

    def _get_mica_score(self, parents1: list[str], parents2: list[str], nb_er_reactome) -> float:
        key = (tuple(sorted(parents1)), tuple(sorted(parents2)))
        if key not in self._mica_cache:
            mica = find_mica_from_several_pathways(
                self.hierarchy_graph, self.root, parents1, parents2,
                self.dico_er_per_pathway, self.nb_er_total,
            )
            self._mica_cache[key] = mica
        mica = self._mica_cache[key]
        return mica, compute_ic_er_pathway(mica, self.dico_er_per_pathway, self.nb_er_total), compute_ic_er_pathway(mica, self.dico_er_per_pathway, nb_er_reactome)

    def _get_max_next_step_score(self, ancestors1: set, ancestors2: set) -> float:
        key = (tuple(sorted(ancestors1)), tuple(sorted(ancestors2)))
        if key not in self._next_step_cache:
            best = 0.0
            for p1, p2 in itertools.product(ancestors1, ancestors2):
                if p1 != p2:
                    weight = self.next_step_lookup.get((p1, p2))
                    if weight and weight > best:
                        best = weight
            self._next_step_cache[key] = best
        if self._next_step_cache[key] != 0:
            return self._next_step_cache[key]
        else:
            return 0

    def score_pair(self, entity1: str, entity2: str, matrix_mica: dict, matrix_nsp: dict, matrix_mica_reactome:dict) -> tuple[float, str] | None:
        """
        Return (score, provenance) for a pair of co-membered entities, or
        None if either entity has no known pathway parent.
        """
        if entity1 not in self.dico_most_precise_parent or entity2 not in self.dico_most_precise_parent:
            return None

        parents1, ancestor_pathways1 = self._get_parents_and_ancestor_pathways(entity1)
        parents2, ancestor_pathways2 = self._get_parents_and_ancestor_pathways(entity2)

        mica, mica_score, mica_score_reactome = self._get_mica_score(parents1, parents2, nb_er_reactome)
        next_step_score = self._get_max_next_step_score(ancestor_pathways1, ancestor_pathways2)
        matrix_mica[(entity1, entity2)] = mica_score
        matrix_mica[(entity2, entity1)] = mica_score
        matrix_mica_reactome[(entity1, entity2)] = mica_score_reactome
        matrix_mica_reactome[(entity2, entity1)] = mica_score_reactome
        matrix_nsp[(entity1, entity2)] = next_step_score
        matrix_nsp[(entity2, entity1)] = next_step_score

        if mica_score >= next_step_score:
            return mica_score, "scoreMICA"
        return next_step_score, "scoreNextStep"


def weight_comembership_clique(
    comembership_clique: pd.DataFrame, scorer: ComembershipScorer, matrix_mica:dict, matrix_nsp:dict, matrix_mica_reactome:dict
) -> pd.DataFrame:
    """Score every entity pair in the co-membership clique table."""
    rows = []
    for entity1, entity2 in comembership_clique.iloc[:, [0, 1]].values:
        result = scorer.score_pair(entity1, entity2, matrix_mica, matrix_nsp, matrix_mica_reactome)
        if result is not None:
            score, provenance = result
            rows.append([entity1, entity2, score, provenance])

    weighted = pd.DataFrame(rows, columns=["entity1", "entity2", "scoreComembership", "provenanceScore"])
    return weighted.sort_values(by="scoreComembership", ascending=False)


def save_weighted_clique(weighted_comembership_clique: pd.DataFrame, counter: int, matrix_mica:dict, matrix_nsp:dict, matrix_mica_reactome:dict) -> None:
    weighted_comembership_clique.to_csv(
        os.path.join(WEIGHTED_CLIQUES_DIR, f"{counter:02d}_WeightedComembershipClique.csv"),
        sep=",", header=True, index=False,
    )
    with open(os.path.join(MATRICES_TOP_PATHWAYS, f"{counter:02d}_MICA_matrix.csv"), 'w', newline='') as fichier_csv:
        writer = csv.writer(fichier_csv)
        writer.writerow(['Source', 'Destination', 'Valeur'])
        for (source, destination), valeur in matrix_mica.items():
            writer.writerow([source, destination, valeur])
    with open(os.path.join(MATRICES_TOP_PATHWAYS, f"{counter:02d}_NSP_matrix.csv"), 'w', newline='') as fichier_csv:
        writer = csv.writer(fichier_csv)
        writer.writerow(['Source', 'Destination', 'Valeur'])
        for (source, destination), valeur in matrix_nsp.items():
            writer.writerow([source, destination, valeur])
    with open(os.path.join(MATRICES_REACTOME, f"{counter:02d}_MICA_matrix_REACTOME.csv"), 'w', newline='') as fichier_csv:
        writer = csv.writer(fichier_csv)
        writer.writerow(['Source', 'Destination', 'Valeur'])
        for (source, destination), valeur in matrix_mica_reactome.items():
            writer.writerow([source, destination, valeur])


# ---------------------------------------------------------------------------
# Per-file pipeline
# ---------------------------------------------------------------------------

def process_pathway_file(owl_file: str, counter: int) -> None:
    """Run the full co-membership weighting pipeline for a single pathway file."""
    start_time = time.time()

    matrix_mica = dict()
    matrix_nsp = dict()
    matrix_mica_reactome = dict()
    pathway_abstraction = load_pathway_abstraction(counter)
    graph = create_networkx_graph(pathway_abstraction)
    hierarchy_graph = extract_subgraph_by_interaction(graph, HIERARCHY_INTERACTION)

    fuseki_process = start_fuseki_server(owl_file, FUSEKI_DATASET_TOP_PATHWAYS)
    try:
        wait_for_fuseki(FUSEKI_ENDPOINT_TOP_PATHWAY)

        # 1-2. Total UniProt count, and per-pathway entity reference lists.
        nb_er_total = count_total_uniprot_ids(FUSEKI_ENDPOINT_TOP_PATHWAY)
        er_per_pathway = load_entity_refs_per_pathway(counter)
        dico_er_per_pathway = build_entity_refs_per_pathway_dict(er_per_pathway)
        
        # 3. IC score per pathway.
        dict_score_ic_pathways = compute_ic_scores(dico_er_per_pathway, nb_er_total)
        save_ic_scores(dict_score_ic_pathways, counter)

        # 4. Next-step pathway score lookup.
        list_pathways = get_ordered_pathway_list(pathway_abstraction)
        next_step_lookup = build_next_step_lookup(pathway_abstraction)

        # 5. Pathway ancestors and "most precise parent" per protein.
        dico_ancestors = build_ancestors_dict(list_pathways, hierarchy_graph)
        up_per_pathway = load_entity_refs_per_pathway(counter)
        dico_parents_up = build_protein_to_pathways_dict(up_per_pathway)

        root = find_root(graph, HIERARCHY_INTERACTION)
        shortest_path_len_from_root = compute_shortest_path_lengths_from_root(
            hierarchy_graph, root, list_pathways
        )
        dico_most_precise_parent = build_most_precise_parent_map(
            dico_parents_up, shortest_path_len_from_root
        )
        print(dico_most_precise_parent)

        # 6. Weight the co-membership clique.
        comembership_clique = load_comembership_clique(counter)
        scorer = ComembershipScorer(
            hierarchy_graph=hierarchy_graph,
            root=root,
            dico_most_precise_parent=dico_most_precise_parent,
            dico_ancestors=dico_ancestors,
            dico_er_per_pathway=dico_er_per_pathway,
            nb_er_total=nb_er_total,
            next_step_lookup=next_step_lookup,
        )
        weighted_comembership_clique = weight_comembership_clique(comembership_clique, scorer, matrix_mica, matrix_nsp, matrix_mica_reactome)
        save_weighted_clique(weighted_comembership_clique, counter, matrix_mica, matrix_nsp, matrix_mica_reactome)

        elapsed_time = time.time() - start_time
        print(f"File {owl_file} processed in {elapsed_time:.2f}s")

    finally:
        fuseki_process.kill()
        fuseki_process.wait()


def main() -> None:
    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(SCORES_DIR, exist_ok=True)
    os.makedirs(WEIGHTED_CLIQUES_DIR, exist_ok=True)

    filelist = sorted(glob.glob(os.path.join(REACTOME_TOP_PATHWAYS_DIR, "*.xml")))
    #filelist = sorted(glob.glob(os.path.join(REACTOME_TOP_PATHWAYS_DIR, "01_Autophagy.xml")))
    for counter, owl_file in enumerate(filelist, start=1):
        process_pathway_file(owl_file, counter)
        counter += 1


if __name__ == "__main__":
    main()