import pandas as pd
import math
import networkx as nx

def compute_similarity(Pathway1, Pathway2, dicoErperPathway, nbERTotal):
    if Pathway1 in dicoErperPathway.keys() and Pathway2 in dicoErperPathway.keys():
        ER_pathway1 = dicoErperPathway[Pathway1]
        ER_pathway2 = dicoErperPathway[Pathway2]
        nb_ER_MICA = (len(set(ER_pathway1)) + len(set(ER_pathway2))) - len(list(set(ER_pathway1) & set(ER_pathway2)))
    else:
        nb_ER_MICA = 0
    if nb_ER_MICA != 0 and nb_ER_MICA != nbERTotal:
        ic = -math.log10(nb_ER_MICA/nbERTotal)
    else:
        ic = 0
    return ic

def create_networkx_graph(abstraction):
    # Créer le graphe en gardant le type d'interaction comme attribut
    PathwayAbstractionGraph = nx.from_pandas_edgelist(
        abstraction,
        source="subPathwayID",
        target="pathwayID",
        edge_attr="interaction",
        create_using=nx.MultiDiGraph()
    )
    print(f"Pathway abstraction loaded in networkx: {PathwayAbstractionGraph}")
    return PathwayAbstractionGraph

def  find_ancestors_and_descendants(G, node, interaction):
    #Since edges are child->parent, we need to reverse the graph to find descendants 
    reversed_graph = G.reverse()
    ancestors = nx.ancestors(reversed_graph, node)
    descendants = nx.descendants(reversed_graph, node)
    return ancestors, descendants

def compute_ic_er_pathway(pathway:str, dico_er_per_pathway:dict, nb_er_total:int):
    nb_er_pathway = len(dico_er_per_pathway[pathway])
    if len(dico_er_per_pathway[pathway]) != nb_er_total:
        return -math.log10(nb_er_pathway/nb_er_total)
    else:
        return 0

def  find_mica_from_several_pathways(G, root, list_pathways1, list_pathways2, dico_er_per_pathway, nbEr_total):
    for p1 in list_pathways1:
        ancestors1, _ = find_ancestors_and_descendants(G, p1, "abstraction:IsAComponentOf")
    for p2 in list_pathways2:
        ancestors2, _ = find_ancestors_and_descendants(G, p2, "abstraction:IsAComponentOf")
    common_ancestors = list(set(ancestors1) & set(ancestors2))
    # direct parent is the MICA
    if list_pathways1 == list_pathways2 and len(list_pathways1) == 1:
        mica = list_pathways1[0]
    elif common_ancestors:
        if p1 in ancestors2 or p2 in ancestors1:
            if p1 in ancestors2:
                mica = p1
            elif p2 in ancestors1:
                mica = p2
        else:
            maxIC = -1
            for ancestor in common_ancestors:
                #_, desc = find_ancestors_and_descendants(G , ancestor, "abstraction:IsAComponentOf")
                #IC = information_content(len(desc), len(G .nodes()))
                IC = compute_ic_er_pathway(ancestor, dico_er_per_pathway, nbEr_total)
                if IC > maxIC:
                    maxIC = IC
                    mica = ancestor
        return mica
    else:
        mica = root
    return mica 

def  resnik_similarity_ER(G, root, node1, node2, dico_er_per_pathway, nbEr_total):
    mica = find_mica_from_several_pathways(G , root, [node1] , [node2], dico_er_per_pathway, nbEr_total)
    if mica:
        return compute_ic_er_pathway(mica, dico_er_per_pathway, nbEr_total)
    return 0 
