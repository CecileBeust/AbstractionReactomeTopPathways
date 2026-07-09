from SPARQLWrapper import SPARQLWrapper, CSV, JSON
import subprocess
import time
import os
import pandas as pd
import glob
from utility import *

current_directory = os.getcwd()
BioPAX_Ontology_file_path = os.path.join(current_directory, 'Data/BioPAX/BioPAXOntology/biopax-level3.owl')
results_dir = os.path.join(current_directory, 'Results/PathwayAbstraction/02_Reactome96Global')

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

# query_up_per_pathway = """ 
# SELECT DISTINCT ?pathwayID ?entityID
# WHERE {
#     VALUES ?db { "UniProt" "UniProt Isoform" }
#     ?pathway rdf:type bp3:Pathway .
#     ?pathway (bp3:pathwayComponent)*|(bp3:pathwayOrder/bp3:stepProcess)* ?interaction .
#     ?interaction rdf:type/(rdfs:subClassOf*) bp3:Interaction .
#     ?pathway bp3:xref ?pathwayXref .
#     ?pathwayXref rdf:type bp3:UnificationXref .
#     ?pathwayXref bp3:db "Reactome" .
#     ?pathwayXref bp3:id ?pathwayID .

#     VALUES ?relation { bp3:left bp3:right bp3:participant bp3:controller }
#     ?interaction ?relation ?entity .

#     ?entity (bp3:component*)|(bp3:memberPhysicalEntity*)|(bp3:component*/bp3:memberPhysicalEntity*)|(bp3:memberPhysicalEntity*/bp3:component*)|(bp3:component*/bp3:memberPhysicalEntity*/bp3:component*) ?entityCompo .
#     ?entityCompo rdf:type/(rdfs:subClassOf*) bp3:PhysicalEntity .
#     ?entityCompo bp3:entityReference ?entityRef .
#     ?entityRef bp3:xref ?entityRefXref .
#     ?entityRefXref rdf:type bp3:UnificationXref .
#     ?entityRefXref bp3:db ?db .
#     ?entityRefXref bp3:id ?entityID .
# } 
# """

# query_nb_up_per_pathway = """ 
# SELECT ?pathwayID (COUNT(DISTINCT ?entityID) AS ?nbID)
# WHERE {
#     VALUES ?db { "UniProt" "UniProt Isoform" }
#     ?pathway rdf:type bp3:Pathway .
#     ?pathway bp3:displayName ?pathwayName .
#     ?pathway (bp3:pathwayComponent)*|(bp3:pathwayOrder/bp3:stepProcess)* ?interaction .
#     ?interaction rdf:type/(rdfs:subClassOf*) bp3:Interaction .
#     ?pathway bp3:xref ?pathwayXref .
#     ?pathwayXref rdf:type bp3:UnificationXref .
#     ?pathwayXref bp3:db "Reactome" .
#     ?pathwayXref bp3:id ?pathwayID .

#     VALUES ?relation { bp3:left bp3:right bp3:participant bp3:controller }
#     ?interaction ?relation ?entity .

#     ?entity (bp3:component*)|(bp3:memberPhysicalEntity*)|(bp3:component*/bp3:memberPhysicalEntity*)|(bp3:memberPhysicalEntity*/bp3:component*)|(bp3:component*/bp3:memberPhysicalEntity*/bp3:component*) ?entityCompo .
#     ?entityCompo rdf:type/(rdfs:subClassOf*) bp3:PhysicalEntity .
#     ?entityCompo bp3:entityReference ?entityRef .
#     ?entityRef bp3:xref ?entityRefXref .
#     ?entityRefXref rdf:type bp3:UnificationXref .
#     ?entityRefXref bp3:db ?db .
#     ?entityRefXref bp3:id ?entityID .
# }
# GROUP BY ?pathwayID
# """

query_up_per_pathway = """ 
SELECT DISTINCT ?pathwayID ?entityID
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

query_nb_up_per_pathway = """ 
SELECT ?pathwayID (COUNT(DISTINCT ?entityID) AS ?nbID)
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
GROUP BY ?pathwayID
"""

query_ids = """
SELECT (COUNT(DISTINCT ?id) AS ?nbID)
WHERE {
VALUES ?db { "UniProt" "UniProt Isoform" }
?entityRef rdf:type/(rdfs:subClassOf*) bp3:ProteinReference .
?entityRef bp3:xref ?entityRefXref .
?entityRefXref bp3:db ?db .
?entityRefXref bp3:id ?id .
}
"""

endpoint = "http://localhost:3030/reactome"
file = "Data/BioPAX/ReactomeBioPAX/Homo_sapiens_v96.owl"
name = "Reactome96"

command = [
    '/home/cbeust/Softwares/JenaFuseki/apache-jena-fuseki-4.9.0/fuseki-server',
    '--file', file,
    '--file', BioPAX_Ontology_file_path,
    '/reactome'
]
print("Fuseki command:", command)

process = subprocess.Popen(command)
time.sleep(30)

sparql = SPARQLWrapper(endpoint)

sparql.setQuery(prefixes + query_is_a_component_of)
sparql.setReturnFormat(CSV)
try:
    results = sparql.query().convert()
    output_filename = os.path.join(results_dir, f"{name}_IsAComponentOf.csv")
    with open(output_filename, 'wb') as f:
        f.write(results)
    print(f"Results saved in {output_filename}")
except Exception as e:
    print(f"Error in IsAComponentOf SPARQL query: {e}")
iscomponentof = pd.read_csv(f"Results/PathwayAbstraction/02_Reactome96Global/{name}_IsAComponentOf.csv", sep=",", header=0)
new_col = ["abstraction:IsAComponentOf"]*len(iscomponentof)
iscomponentof.insert(loc=1, column="interaction", value=new_col)
iscomponentof.to_csv(f"Results/PathwayAbstraction/02_Reactome96Global/{name}_IsAComponentOf.csv", sep=",", header=True, index=False)

sparql.setQuery(prefixes + query_next_step_pathway)
sparql.setReturnFormat(CSV)
try:
    results = sparql.query().convert()
    output_filename = os.path.join(results_dir, f"{name}_NextStepPathway.csv")
    with open(output_filename, 'wb') as f:
        f.write(results)
    print(f"Results saved in {output_filename}")
except Exception as e:
    print(f"Error in NextStepPathway query: {e}")

nextsteppathway = pd.read_csv(f"Results/PathwayAbstraction/02_Reactome96Global/{name}_NextStepPathway.csv", sep=",", header=0)
new_col = ["abstraction:NextStepPathway"]*len(nextsteppathway)
nextsteppathway.insert(loc=1, column="interaction", value=new_col)
nextsteppathway.to_csv(f"Results/PathwayAbstraction/02_Reactome96Global/{name}_NextStepPathway.csv", sep=",", header=True, index=False)

try:
    q1 = pd.read_csv(os.path.join(results_dir, f"{name}_IsAComponentOf.csv"), header=None)
    q2 = pd.read_csv(os.path.join(results_dir, f"{name}_NextStepPathway.csv"), header=None)
    q1 = q1.drop(0).reset_index(drop=True)
    q2 = q2.drop(0).reset_index(drop=True)
    concat_df = pd.concat([q1, q2], ignore_index=True)
    output_csv = os.path.join(results_dir, f"{name}_PathwayAbstraction.tsv")
    concat_df.to_csv(output_csv, sep=',', header=False, index=False)
    print(f"CSV output saved in {output_csv}")
except Exception as e:
    print(f"Error: {e}")


sparql.setQuery(prefixes + query_up_per_pathway)
sparql.setReturnFormat(CSV)
try:
    results = sparql.query().convert()
    output_filename = os.path.join(os.path.join(current_directory, 'Results/UtilityFiles'), f"{name}_UpPerPathway.csv")
    with open(output_filename, 'wb') as f:
        f.write(results)
    print(f"Results saved in {output_filename}")
except Exception as e:
    print(f"Error: {e}")

sparql.setQuery(prefixes + query_nb_up_per_pathway)
sparql.setReturnFormat(CSV)
try:
    results = sparql.query().convert()
    output_filename = os.path.join(os.path.join(current_directory, 'Results/UtilityFiles'), f"{name}_NbUpPerPathway.csv")
    with open(output_filename, 'wb') as f:
        f.write(results)
    print(f"Results saved in {output_filename}")
except Exception as e:
    print(f"Error: {e}")

isacomponentof = pd.read_csv(f"Results/PathwayAbstraction/02_Reactome96Global/{name}_IsAComponentOf.csv", sep=',', header=0)
nextsteppathway = pd.read_csv(f"Results/PathwayAbstraction/02_Reactome96Global/{name}_NextStepPathway.csv", sep=',', header=0)
G = create_networkx_graph(isacomponentof)
root = [node for node in G.nodes() if G.out_degree(node) == 0]

up_per_pathway = pd.read_csv(f"Results/UtilityFiles/{name}_UpPerPathway.csv", sep=",", header=0)
print(up_per_pathway.head())

dico_er_per_pathway = dict()
for index, row in up_per_pathway.iterrows():
    pathway = row[0]
    id = row[1]
    if not pathway in dico_er_per_pathway.keys():
        dico_er_per_pathway[pathway] = [id]
    else:
        dico_er_per_pathway[pathway] += [id]

sparql.setQuery(prefixes + query_ids)
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

weighted_graph.to_csv(f"Results/PathwayAbstraction/02_Reactome96Global/{name}_IsAComponentOf_ResnikER.csv", sep=",", header=True, index=False)

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
weighted_graph.to_csv(f"Results/PathwayAbstraction/02_Reactome96Global/{name}_NextStepPathway_ERcontent.csv", sep=",", header=True, index=False)

# CONCATENATE GRAPHS
weighted_is_a_component_of = pd.read_csv(f"Results/PathwayAbstraction/02_Reactome96Global/{name}_IsAComponentOf_ResnikER.csv", sep=",", header=0)
print(weighted_is_a_component_of.head())

weighted_next_step_pathway = pd.read_csv(f"Results/PathwayAbstraction/02_Reactome96Global/{name}_NextStepPathway_ERcontent.csv", sep=",", header=0)
print(weighted_next_step_pathway.head())

global_abstraction = pd.concat([weighted_is_a_component_of, weighted_next_step_pathway], axis=0)
global_abstraction.to_csv(f"Results/PathwayAbstraction/02_Reactome96Global/{name}_WeightedPathwayAbstraction.csv", sep=",", index=False)


process.kill()
time.sleep(30)

