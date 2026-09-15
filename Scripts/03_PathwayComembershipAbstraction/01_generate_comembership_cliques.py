from SPARQLWrapper import SPARQLWrapper, TURTLE, JSON, CSV
import subprocess
import time
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import glob
from itertools import combinations

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

# query="""
# SELECT DISTINCT ?id1 ?id2
# WHERE {
#     VALUES ?db1 { "UniProt" "UniProt Isoform" }
#     ?entityRef1 rdf:type/(rdfs:subClassOf*) bp3:ProteinReference .
#     ?entityRef1 bp3:xref ?entityRef1Xref .
#     ?entityRefXref1 rdf:type bp3:UnificationXref .
#     ?entityRefXref1 bp3:db ?db1 .
#     ?entityRefXref1 bp3:id ?id1 .
    
#     VALUES ?db2 { "UniProt" "UniProt Isoform" }
#     ?entityRef2 rdf:type/(rdfs:subClassOf*) bp3:ProteinReference .
#     ?entityRef2 bp3:xref ?entityRef2Xref .
#     ?entityRef2Xref rdf:type bp3:UnificationXref .
#     ?entityRef2Xref bp3:db ?db2 .
#     ?entityRef2Xref bp3:id ?id2 .
#     FILTER (?id1 < ?id2)
# }
# """

query = """
SELECT DISTINCT ?entityID
WHERE {
    VALUES ?db { "UniProt" "UniProt Isoform" }

    ?pathway rdf:type bp3:Pathway .
    ?pathway (bp3:pathwayComponent | bp3:pathwayOrder/bp3:stepProcess)* ?interaction .
    ?interaction rdf:type/(rdfs:subClassOf*) bp3:Interaction .
    ?pathway bp3:xref ?pathwayXref .
    ?pathwayXref rdf:type bp3:UnificationXref ;
                    bp3:db "Reactome" ;
                    bp3:id ?pathwayID .

    VALUES ?relation { bp3:left bp3:right bp3:participant bp3:controller }
    ?interaction ?relation ?entity .
    
    ?entity (bp3:component | bp3:memberPhysicalEntity)* ?entityCompo .

    ?entityCompo rdf:type/(rdfs:subClassOf*) bp3:PhysicalEntity .
    ?entityCompo bp3:entityReference ?entityRef .
    ?entityRef bp3:xref ?entityRefXref .
    ?entityRefXref rdf:type bp3:UnificationXref ;
                    bp3:db ?db ;
                    bp3:id ?entityID .
}
"""

endpoint = "http://localhost:3030/top_pathway"
counter = 1
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
        #output_filename = os.path.join(results_dir, f"{counter:02d}_ComembershipClique.csv")
        output_filename = os.path.join(results_dir, f"{counter:02d}_ListProteins.csv")
        with open(output_filename, 'wb') as f:
            f.write(results)
        print(f"Results saved in {output_filename}")
    except Exception as e:
        print(f"Error in SPARQL query: {e}")

    # added
    list_proteins = os.path.join(current_directory, f"Results/PathwayComembership/01_ComembershipCliques/{counter:02d}_ListProteins.csv")
    df = pd.read_csv(list_proteins, sep=",", header=0)
    proteins = df['entityID'].tolist()
    pairs = list(combinations(proteins, 2))
    df_pairs = pd.DataFrame(pairs, columns=['id1', 'id2'])
    output_file = os.path.join(results_dir,  f"{counter:02d}_ComembershipCliqueCleaned.csv")
    df_pairs.to_csv(output_file, index=False)
    print(len(df_pairs))

    elapsed_time = time.time() - start_time
    print(f"Time taken for {owl_file}: {elapsed_time:.2f} seconds")
    process.kill()
    time.sleep(30)
    counter += 1