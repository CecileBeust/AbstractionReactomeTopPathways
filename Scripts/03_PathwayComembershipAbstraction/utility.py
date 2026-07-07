import networkx as nx
import math

def  find_ancestors_and_descendants(G, node):
    #Since edges are child->parent, we need to reverse the graph to find descendants 
    reversed_graph = G.reverse()
    ancestors = nx.ancestors(reversed_graph, node)
    descendants = nx.descendants(reversed_graph, node)
    return ancestors, descendants

def create_networkx_graph(abstraction):
    # Créer le graphe en gardant le type d'interaction comme attribut
    PathwayAbstractionGraph = nx.from_pandas_edgelist(
        abstraction,
        source="source",
        target="target",
        edge_attr="interaction",
        create_using=nx.MultiDiGraph()
    )
    print(f"Pathway abstraction loaded in networkx: {PathwayAbstractionGraph}")
    return PathwayAbstractionGraph

def get_pathways_and_ancestor(pathways:list, dico_ancestors:dict):
    new_set = set()
    all_ancestors = list()
    for pathway in pathways:
        ancestors = dico_ancestors[pathway]
        all_ancestors += ancestors
    all_ancestors = set(all_ancestors)
    new_set = all_ancestors.union(pathways)
    return new_set

def get_edge_interaction(G, u, v) -> str:
    data = G.get_edge_data(u, v, default={})
    if 0 in data:
        return data[0].get("interaction", "")
    if "edge_type" in data:
        return data["edge_type"]
    return data.get("interaction", "")

def extract_subgraph_by_interaction(G: nx.Graph, interaction: str) -> nx.DiGraph:
    subgraph = nx.DiGraph()
    for u, v, data in G.edges(data=True):
        if get_edge_interaction(G, u, v) == interaction:
            subgraph.add_edge(u, v, edge_type=interaction)
    return subgraph

def find_root(G: nx.Graph, interaction: str):
    root = None
    for node in G.nodes():
        out_degree = sum(1 for u, v, edge_data in G.out_edges(node, data=True) if edge_data.get('interaction') == interaction)
        if out_degree == 0:
            root = node
            break
    return root

def compute_ic_er_pathway(pathway:str, dico_er_per_pathway:dict, nb_er_total:int):
    nb_er_pathway = len(dico_er_per_pathway[pathway])
    if len(dico_er_per_pathway[pathway]) != nb_er_total:
        return -math.log(nb_er_pathway/nb_er_total)
    else:
        return 0

def  find_mica_from_several_pathways(G, root, list_pathways1, list_pathways2, dico_er_per_pathway, nbEr_total):
    for p1 in list_pathways1:
        ancestors1, _ = find_ancestors_and_descendants(G, p1)
    for p2 in list_pathways2:
        ancestors2, _ = find_ancestors_and_descendants(G, p2)
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