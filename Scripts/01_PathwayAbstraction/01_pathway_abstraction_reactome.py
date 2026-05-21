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
results_dir = os.path.join(current_directory, 'Results/PathwayAbstraction')

os.makedirs(results_dir, exist_ok=True)

prefixes = """
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX bp3: <http://www.biopax.org/release/biopax-level3.owl#>
"""

query_is_a_component_of = """
SELECT DISTINCT ?subPathwayID ?pathwayID
WHERE {
    ?pathway rdf:type bp3:Pathway .
    ?pathway bp3:pathwayComponent ?subPathway .
    ?subPathway rdf:type bp3:Pathway .

    ?pathway bp3:xref ?pathwayXref .
    ?pathwayXref rdf:type bp3:UnificationXref .
    ?pathwayXref bp3:db "Reactome" .
    ?pathwayXref bp3:id ?pathwayID .
    
    ?subPathway bp3:xref ?subPathwayXref .
    ?subPathwayXref rdf:type bp3:UnificationXref .
    ?subPathwayXref bp3:db "Reactome" .
    ?subPathwayXref bp3:id ?subPathwayID .
}
"""

query_next_step_pathway = """
SELECT DISTINCT ?pathwayID ?nextPathwayID
WHERE {
    ?pathway rdf:type bp3:Pathway .
    ?nextPathway rdf:type bp3:Pathway .
    
    ?pathway bp3:pathwayOrder ?pathwayStep .
    ?nextPathway bp3:pathwayOrder ?nextStep .
    
    ?pathwayStep bp3:nextStep ?nextStep .
    
    FILTER (?pathway != ?nextPathway)

    ?pathway bp3:xref ?pathwayXref .
    ?pathwayXref rdf:type bp3:UnificationXref .
    ?pathwayXref bp3:db "Reactome" .
    ?pathwayXref bp3:id ?pathwayID .
    
    ?nextPathway bp3:xref ?nextPathwayXref .
    ?nextPathwayXref rdf:type bp3:UnificationXref .
    ?nextPathwayXref bp3:db "Reactome" .
    ?nextPathwayXref bp3:id ?nextPathwayID .
}
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

    sparql.setQuery(prefixes + query_is_a_component_of)
    sparql.setReturnFormat(CSV)
    try:
        results = sparql.query().convert()
        output_filename = os.path.join(results_dir, f"{counter:02d}_IsAComponentOf.csv")
        with open(output_filename, 'wb') as f:
            f.write(results)
        print(f"Results saved in {output_filename}")
    except Exception as e:
        print(f"Error in IsAComponentOf SPARQL query: {e}")
    iscomponentof = pd.read_csv(f"Results/PathwayAbstraction/{counter:02d}_IsAComponentOf.csv", sep=",", header=0)
    new_col = ["abstraction:IsAComponentOf"]*len(iscomponentof)
    iscomponentof.insert(loc=1, column="interaction", value=new_col)
    iscomponentof.to_csv(f"Results/PathwayAbstraction/{counter:02d}_IsAComponentOf.csv", sep=",", header=True, index=False)

    sparql.setQuery(prefixes + query_next_step_pathway)
    sparql.setReturnFormat(CSV)
    try:
        results = sparql.query().convert()
        output_filename = os.path.join(results_dir, f"{counter:02d}_NextStepPathway.csv")
        with open(output_filename, 'wb') as f:
            f.write(results)
        print(f"Results saved in {output_filename}")
    except Exception as e:
        print(f"Error in NextStepPathway query: {e}")

    nextsteppathway = pd.read_csv(f"Results/PathwayAbstraction/{counter:02d}_NextStepPathway.csv", sep=",", header=0)
    new_col = ["abstraction:NextStepPathway"]*len(nextsteppathway)
    nextsteppathway.insert(loc=1, column="interaction", value=new_col)
    nextsteppathway.to_csv(f"Results/PathwayAbstraction/{counter:02d}_NextStepPathway.csv", sep=",", header=True, index=False)

    try:
        q1 = pd.read_csv(os.path.join(results_dir, f"{counter:02d}_IsAComponentOf.csv"), header=None)
        q2 = pd.read_csv(os.path.join(results_dir, f"{counter:02d}_NextStepPathway.csv"), header=None)
        q1 = q1.drop(0).reset_index(drop=True)
        q2 = q2.drop(0).reset_index(drop=True)
        concat_df = pd.concat([q1, q2], ignore_index=True)
        output_csv = os.path.join(results_dir, f"{counter:02d}_PathwayAbstraction.tsv")
        concat_df.to_csv(output_csv, sep=',', header=False, index=False)
        print(f"CSV output saved in {output_csv}")
    except Exception as e:
        print(f"Error: {e}")

    process.kill()
    time.sleep(30)
    counter += 1


