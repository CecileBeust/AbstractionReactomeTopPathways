from SPARQLWrapper import SPARQLWrapper, JSON
import subprocess
import time
import os
import pandas as pd
import glob

current_directory = os.getcwd()
BioPAX_Ontology_file_path = os.path.join(current_directory, '../Data/BioPAX/BioPAXOntology/biopax-level3.owl')
ReactomeBioPAX_file_path = os.path.join(current_directory, '../Data/BioPAX/ReactomeTopPathways')
results_dir = os.path.join(current_directory, '../Results/')

prefixes = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX bp3: <http://www.biopax.org/release/biopax-level3.owl#>
"""

query_nb_triples = """
SELECT (COUNT(DISTINCT *) AS ?nbTriples)
WHERE {
    ?s ?p ?o.
}
"""

query_nb_pathways = """
SELECT (COUNT(DISTINCT ?pathway) AS ?nbPathway)
WHERE {
    ?pathway rdf:type bp3:Pathway .
}
"""

dico_top_pathway_ids = {
    "01": "Autophagy",
    "02": "CellCycle",
    "03": "CellCellCommunication",
    "04": "CellularResponseToStimuli",
    "05": "ChromatinOrganization",
    "06": "CircadianClock",
    "07": "DevelopmentalBiology",
    "08": "DigestionAndAbsorption",
    "09": "Disease",
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

endpoint = "http://localhost:3030/top_pathway"
counter = 1
df_metrics_pathways = pd.DataFrame(columns=["top_pathway", "number_of_pathways", "number_of_triples"])
filelist = glob.glob(os.path.join(ReactomeBioPAX_file_path, '*.xml'))

counter = 1
for owl_file in sorted(filelist):
    print(owl_file, counter)
    print(f"Processing file: {owl_file}")
    command = [
        '/home/cbeust/Softwares/JenaFuseki/apache-jena-fuseki-4.9.0/fuseki-server',
        '--file', owl_file,
        '--file', BioPAX_Ontology_file_path,
        '/top_pathway'
    ]
    print("Fuseki command:", command)

    process = subprocess.Popen(command)
    time.sleep(30)

    sparql = SPARQLWrapper(endpoint)
    sparql.setQuery(prefixes + query_nb_triples)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    nb_triples = int(results["results"]["bindings"][0]["nbTriples"]["value"])
    print(nb_triples)

    sparql.setQuery(prefixes + query_nb_pathways)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    nb_pathways = int(results["results"]["bindings"][0]["nbPathway"]["value"])
    print(nb_pathways)

    df_metrics_pathways.at[counter, "top_pathway"] = dico_top_pathway_ids[f"{counter:02d}"]
    df_metrics_pathways.at[counter, "number_of_pathways"] = nb_pathways
    df_metrics_pathways.at[counter, "number_of_triples"] = nb_triples

    process.kill()
    time.sleep(30)
    counter += 1

df_metrics_pathways.to_csv("../Results/PathwayAbstraction/MetricsTopPathways.csv", sep=",", header=0)
