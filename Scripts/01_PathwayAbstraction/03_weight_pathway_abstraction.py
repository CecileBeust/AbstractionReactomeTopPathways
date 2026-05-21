from SPARQLWrapper import SPARQLWrapper, CSV, JSON
import subprocess
import time
import os
import pandas as pd
import glob
from utility import *
from xml.etree import ElementTree as ET  # Pour parser le XML

current_directory = os.getcwd()
results_dir = os.path.join(current_directory, 'Results/PathwayAbstraction')
BioPAX_Ontology_file_path = os.path.join(current_directory, 'Data/BioPAX/BioPAXOntology/biopax-level3.owl')
ReactomeBioPAX_file_path = os.path.join(current_directory, 'Data/BioPAX/ReactomeTopPathways')

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

biopax_filelist = glob.glob(os.path.join(ReactomeBioPAX_file_path, '*.xml'))
abstraction_filelist = glob.glob(os.path.join(results_dir, '*IsAComponentOf.csv'))

dico_top_pathway_ids = {
    "1": "Autophagy",
    "2": "CellCycle",
    "3": "CellCellCommunication",
    "4": "CellularResponseToStimuli",
    "5": "ChromatinOrganization",
    "6": "CircadianClock",
    "7": "DevelopmentalBiology",
    "8": "DigestionAndAbsorption",
    "9": "Disease",
    "10": "DNARepair",
    "11": "DNAReplication",
    "12": "DrugADME",
    "13": "ExtracellularMatrixOrganization",
    "14": "GeneExpression(Transcription)",
    "15": "Hemostasis",
    "16": "ImmuneSystem",
    "17": "Metabolism",
    "18": "MetabolismOfProteins",
    "19": "MetabolismOfRNA",
    "20": "MuscleContraction",
    "21": "NeuronalSystem",
    "22": "OrganelleBiogenesisAndMaintenance",
    "23": "ProgrammedCellDeath",
    "24": "ProteinLocalization",
    "25": "Reproduction",
    "26": "SensoryPerception",
    "27": "SignalTransduction",
    "28": "TransportOfSmallMolecules",
    "29": "VesicleMediatedTransport"
}

for counter in range(29,30):
    biopax = f"Data/BioPAX/ReactomeTopPathways/{counter:02d}_{dico_top_pathway_ids[str(counter)]}.xml"
    isacomponentof = pd.read_csv(f"Results/PathwayAbstraction/{counter:02d}_IsAComponentOf.csv", sep=',', header=0)
    nextsteppathway = pd.read_csv(f"Results/PathwayAbstraction/{counter:02d}_NextStepPathway.csv", sep=',', header=0)
    G = create_networkx_graph(isacomponentof)
    root = [node for node in G.nodes() if G.out_degree(node) == 0]

    up_per_pathway = pd.read_csv(f"Results/UtilityFiles/{counter:02d}_UpPerPathway.csv", sep=",", header=0)
    print(up_per_pathway.head())

    dico_er_per_pathway = dict()
    for index, row in up_per_pathway.iterrows():
        pathway = row[0]
        id = row[1]
        if not pathway in dico_er_per_pathway.keys():
            dico_er_per_pathway[pathway] = [id]
        else:
            dico_er_per_pathway[pathway] += [id]

    command = [
        '/home/cbeust/Softwares/JenaFuseki/apache-jena-fuseki-4.9.0/fuseki-server',
        '--file', biopax,
        '--file', BioPAX_Ontology_file_path,
        '/top_pathway'
    ]
    print("Fuseki command:", command)

    process = subprocess.Popen(command)
    time.sleep(30)
    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(prefixes + query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    nbEr_total = int(results["results"]["bindings"][0]["nbID"]["value"])
    print(nbEr_total)

    dico_isacomponentof_resnik_er = dict()
    for index, row in isacomponentof.iterrows():
        pathway1 = row[0]
        pathway2 = row[2]
        interaction = row[1]
        if interaction == "abstraction:IsAComponentOf":
            score = resnik_similarity_ER(G, root[0], pathway1, pathway2, dico_er_per_pathway, nbEr_total)
            key = f"{str(pathway1)}-<>{str(pathway2)}"
            dico_isacomponentof_resnik_er[key] = score
    print(dico_isacomponentof_resnik_er)

    weighted_graph = pd.DataFrame(columns=['source', 'interaction', 'target', 'weight'])
    counter_rows = 0
    for index, row in isacomponentof.iterrows():
        weighted_graph.at[counter_rows, 'source'] = row[0]
        weighted_graph.at[counter_rows, 'interaction'] = row[1]
        weighted_graph.at[counter_rows, 'target'] = row[2]
        key = f'{str(row[0])}-<>{str(row[2])}'
        weight = dico_isacomponentof_resnik_er[key]
        weighted_graph.at[counter_rows, 'weight'] = weight
        counter_rows += 1

    weighted_graph.to_csv(f"Results/PathwayAbstraction/{counter:02d}_IsAComponentOf_ResnikER.csv", sep=",", header=True, index=False)

    process.kill()
    time.sleep(30)

    # WEIGHT NEXTSTEPPATHWAY EDGES
    dico_next_step_er = dict()
    distribution_weight_er = list()

    for item, row in nextsteppathway.iterrows():
        if row[1] == "abstraction:NextStepPathway":
            Pathway1 = row[0]
            Pathway2 = row[2]
            key = f"{str(Pathway1)}->{str(Pathway2)}"
            weight_er = compute_similarity(Pathway1, Pathway2, dico_er_per_pathway, nbEr_total)
            distribution_weight_er.append(weight_er)
            dico_next_step_er[key] = weight_er

    weighted_graph = pd.DataFrame(columns=['source', 'interaction', 'target', 'weight'])
    counter_rows = 0
    for item, row in nextsteppathway.iterrows():
        if row[1] == "abstraction:NextStepPathway":
            weighted_graph.at[counter_rows, 'source'] = row[0]
            weighted_graph.at[counter_rows, 'interaction'] = row[1]
            weighted_graph.at[counter_rows, 'target'] = row[2]
            key = f"{str(row[0])}->{str(row[2])}"
            weight = dico_next_step_er[key]
            weighted_graph.at[counter_rows, 'weight'] = weight
        counter_rows += 1
    weighted_graph.to_csv(f"Results/PathwayAbstraction/{counter:02d}_NextStepPathway_ERcontent.csv", sep=",", header=True, index=False)

    # CONCATENATE GRAPHS
    weighted_is_a_component_of = pd.read_csv(f"Results/PathwayAbstraction/{counter:02d}_IsAComponentOf_ResnikER.csv", sep=",", header=0)
    print(weighted_is_a_component_of.head())

    weighted_next_step_pathway = pd.read_csv(f"Results/PathwayAbstraction/{counter:02d}_NextStepPathway_ERcontent.csv", sep=",", header=0)
    print(weighted_next_step_pathway.head())

    global_abstraction = pd.concat([weighted_is_a_component_of, weighted_next_step_pathway], axis=0)
    global_abstraction.to_csv(f"Results/PathwayAbstraction/{counter:02d}_WeightedPathwayAbstraction.csv", sep=",", index=False)
