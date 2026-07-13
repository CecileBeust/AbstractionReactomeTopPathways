from SPARQLWrapper import SPARQLWrapper, TURTLE, JSON, CSV
import subprocess
import time
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import glob

current_directory = os.getcwd()
BioPAX_Ontology_file_path = os.path.join(current_directory, 'Data/BioPAX/BioPAXOntology/biopax-level3.owl')
ReactomeBioPAX_file_path = os.path.join(current_directory, 'Data/BioPAX/ReactomeTopPathways')
results_dir = os.path.join(current_directory, 'Results/PathwayComembership/01_ComembershipCliques')

os.makedirs(results_dir, exist_ok=True)

prefixes = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX bp3: <http://www.biopax.org/release/biopax-level3.owl#>
"""

query="""
SELECT DISTINCT ?id1 ?id2
WHERE {
    VALUES ?db1 { "UniProt" "UniProt Isoform" }
    ?entityRef1 rdf:type/(rdfs:subClassOf*) bp3:ProteinReference .
    ?entityRef1 bp3:xref ?entityRef1Xref .
    ?entityRefXref1 rdf:type bp3:UnificationXref .
    ?entityRefXref1 bp3:db ?db1 .
    ?entityRefXref1 bp3:id ?id1 .
    
    VALUES ?db2 { "UniProt" "UniProt Isoform" }
    ?entityRef2 rdf:type/(rdfs:subClassOf*) bp3:ProteinReference .
    ?entityRef2 bp3:xref ?entityRef2Xref .
    ?entityRef2Xref rdf:type bp3:UnificationXref .
    ?entityRef2Xref bp3:db ?db2 .
    ?entityRef2Xref bp3:id ?id2 .
    FILTER (?id1 < ?id2)
}
"""

endpoint = "http://localhost:3030/top_pathway"
counter = 19
filelist = glob.glob(os.path.join(ReactomeBioPAX_file_path, '*.xml'))
#filelist = glob.glob(os.path.join(ReactomeBioPAX_file_path, '01_Autophagy.xml'))

for owl_file in sorted(filelist):
    start_time = time.time()
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

    sparql.setQuery(prefixes + query)
    sparql.setReturnFormat(CSV)
    try:
        results = sparql.query().convert()
        output_filename = os.path.join(results_dir, f"{counter:02d}_ComembershipClique.csv")
        with open(output_filename, 'wb') as f:
            f.write(results)
        print(f"Results saved in {output_filename}")
    except Exception as e:
        print(f"Error in SPARQL query: {e}")

    elapsed_time = time.time() - start_time
    print(f"Time taken for {owl_file}: {elapsed_time:.2f} seconds")
    process.kill()
    time.sleep(30)
    counter += 1