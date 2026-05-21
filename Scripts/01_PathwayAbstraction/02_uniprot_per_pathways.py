from SPARQLWrapper import SPARQLWrapper, CSV
import subprocess
import time
import os
import pandas as pd
import glob
from utility import *

current_directory = os.getcwd()
BioPAX_Ontology_file_path = os.path.join(current_directory, 'Data/BioPAX/BioPAXOntology/biopax-level3.owl')
ReactomeBioPAX_file_path = os.path.join(current_directory, 'Data/BioPAX/ReactomeTopPathways')
output_dir = os.path.join(current_directory, 'Results/UtilityFiles')

prefixes = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX bp3: <http://www.biopax.org/release/biopax-level3.owl#>
"""

query_up_per_pathway = """ 
SELECT DISTINCT ?pathwayID ?entityID
WHERE {
    VALUES ?db { "UniProt" "UniProt Isoform" }
    ?pathway rdf:type bp3:Pathway .
    ?pathway (bp3:pathwayComponent)*|(bp3:pathwayOrder/bp3:stepProcess)* ?interaction .
    ?interaction rdf:type/(rdfs:subClassOf*) bp3:Interaction .
    ?pathway bp3:xref ?pathwayXref .
    ?pathwayXref rdf:type bp3:UnificationXref .
    ?pathwayXref bp3:db "Reactome" .
    ?pathwayXref bp3:id ?pathwayID .

    VALUES ?relation { bp3:left bp3:right bp3:participant bp3:controller }
    ?interaction ?relation ?entity .

    ?entity (bp3:component*)|(bp3:memberPhysicalEntity*)|(bp3:component*/bp3:memberPhysicalEntity*)|(bp3:memberPhysicalEntity*/bp3:component*)|(bp3:component*/bp3:memberPhysicalEntity*/bp3:component*) ?entityCompo .
    ?entityCompo rdf:type/(rdfs:subClassOf*) bp3:PhysicalEntity .
    ?entityCompo bp3:entityReference ?entityRef .
    ?entityRef bp3:xref ?entityRefXref .
    ?entityRefXref rdf:type bp3:UnificationXref .
    ?entityRefXref bp3:db ?db .
    ?entityRefXref bp3:id ?entityID .
} 
"""

query_nb_up_per_pathway = """ 
SELECT ?pathwayID (COUNT(DISTINCT ?entityID) AS ?nbID)
WHERE {
    VALUES ?db { "UniProt" "UniProt Isoform" }
    ?pathway rdf:type bp3:Pathway .
    ?pathway bp3:displayName ?pathwayName .
    ?pathway (bp3:pathwayComponent)*|(bp3:pathwayOrder/bp3:stepProcess)* ?interaction .
    ?interaction rdf:type/(rdfs:subClassOf*) bp3:Interaction .
    ?pathway bp3:xref ?pathwayXref .
    ?pathwayXref rdf:type bp3:UnificationXref .
    ?pathwayXref bp3:db "Reactome" .
    ?pathwayXref bp3:id ?pathwayID .

    VALUES ?relation { bp3:left bp3:right bp3:participant bp3:controller }
    ?interaction ?relation ?entity .

    ?entity (bp3:component*)|(bp3:memberPhysicalEntity*)|(bp3:component*/bp3:memberPhysicalEntity*)|(bp3:memberPhysicalEntity*/bp3:component*)|(bp3:component*/bp3:memberPhysicalEntity*/bp3:component*) ?entityCompo .
    ?entityCompo rdf:type/(rdfs:subClassOf*) bp3:PhysicalEntity .
    ?entityCompo bp3:entityReference ?entityRef .
    ?entityRef bp3:xref ?entityRefXref .
    ?entityRefXref rdf:type bp3:UnificationXref .
    ?entityRefXref bp3:db ?db .
    ?entityRefXref bp3:id ?entityID .
}
GROUP BY ?pathwayID
"""

endpoint = "http://localhost:3030/top_pathway"
counter = 1
filelist = glob.glob(os.path.join(ReactomeBioPAX_file_path, '*.xml'))

for owl_file in sorted(filelist):
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

    sparql.setQuery(prefixes + query_up_per_pathway)
    sparql.setReturnFormat(CSV)
    try:
        results = sparql.query().convert()
        output_filename = os.path.join(output_dir, f"{counter:02d}_UpPerPathway.csv")
        with open(output_filename, 'wb') as f:
            f.write(results)
        print(f"Results saved in {output_filename}")
    except Exception as e:
        print(f"Error: {e}")

    sparql.setQuery(prefixes + query_nb_up_per_pathway)
    sparql.setReturnFormat(CSV)
    try:
        results = sparql.query().convert()
        output_filename = os.path.join(output_dir, f"{counter:02d}_NbUpPerPathway.csv")
        with open(output_filename, 'wb') as f:
            f.write(results)
        print(f"Results saved in {output_filename}")
    except Exception as e:
        print(f"Error: {e}")
    
    process.kill()
    time.sleep(30)
    counter += 1