"""
utility.py

Helper functions for pathway-hierarchy graph analysis and
Information-Content (IC) based semantic similarity scoring, used to
weight pathway comembership cliques built from Reactome BioPAX data.

The pathway hierarchy is represented as a NetworkX graph where an edge
u -> v with interaction "abstraction:IsAComponentOf" means
"u is a component of v" (i.e. child -> parent).
"""

from __future__ import annotations
import math
from typing import Any
import networkx as nx
import pandas as pd


# ---------------------------------------------------------------------------
# Graph construction
# ---------------------------------------------------------------------------

def create_networkx_graph(abstraction: pd.DataFrame) -> nx.MultiDiGraph:
    """
    Build a MultiDiGraph from a pathway-abstraction edge table.

    Parameters
    ----------
    abstraction : pandas.DataFrame
        Must contain "source", "target" and "interaction" columns.

    Returns
    -------
    networkx.MultiDiGraph
    """
    graph = nx.from_pandas_edgelist(
        abstraction,
        source="source",
        target="target",
        edge_attr="interaction",
        create_using=nx.MultiDiGraph(),
    )
    print(f"Pathway abstraction loaded in networkx: {graph}")
    return graph


def get_edge_interaction(graph: nx.Graph, u: Any, v: Any) -> str:
    """
    Return the interaction type stored on edge (u, v).

    Handles both a MultiDiGraph (interaction stored under key 0) and a
    plain DiGraph (interaction stored as "edge_type" or "interaction").
    """
    data = graph.get_edge_data(u, v, default={})
    if 0 in data:
        return data[0].get("interaction", "")
    if "edge_type" in data:
        return data["edge_type"]
    return data.get("interaction", "")


def extract_subgraph_by_interaction(graph: nx.Graph, interaction: str) -> nx.DiGraph:
    """
    Extract a simple DiGraph containing only the edges whose interaction
    type matches `interaction` (e.g. "abstraction:IsAComponentOf").
    """
    subgraph = nx.DiGraph()
    for u, v in graph.edges():
        if get_edge_interaction(graph, u, v) == interaction:
            subgraph.add_edge(u, v, edge_type=interaction)
    return subgraph


def find_root(graph: nx.Graph, interaction: str) -> Any | None:
    """
    Find the root pathway of the hierarchy: the node that has no
    outgoing edge of the given interaction type (nothing it is a
    component of).
    """
    for node in graph.nodes():
        out_degree = sum(
            1
            for _, _, edge_data in graph.out_edges(node, data=True)
            if edge_data.get("interaction") == interaction
        )
        if out_degree == 0:
            return node
    return None


# ---------------------------------------------------------------------------
# Ancestors / descendants fetcher functions
# ---------------------------------------------------------------------------

def find_ancestors_and_descendants(graph: nx.Graph, node: Any) -> tuple[set, set]:
    """
    Return (ancestors, descendants) of `node` in the hierarchy.

    Edges are child -> parent, so the graph must be reversed to walk
    "downwards" (towards descendants).
    """
    reversed_graph = graph.reverse()
    ancestors = nx.ancestors(reversed_graph, node)
    descendants = nx.descendants(reversed_graph, node)
    return ancestors, descendants


def get_pathways_and_ancestor(pathways: list, dico_ancestors: dict) -> set:
    """
    Return the union of `pathways` with all of their known ancestors,
    using a precomputed {pathway: ancestors} lookup table.
    """
    all_ancestors: list = []
    for pathway in pathways:
        all_ancestors.extend(dico_ancestors[pathway])
    return set(all_ancestors).union(pathways)


# ---------------------------------------------------------------------------
# Information Content (IC) scoring
# ---------------------------------------------------------------------------

def compute_ic_er_pathway(pathway: str, dico_er_per_pathway: dict, nb_er_total: int) -> float:
    """
    Compute the Information Content (IC) of a pathway, based on how many
    entity references (UniProt IDs) it contains relative to the total
    number of entity references across the whole pathway set.

    IC = -log(n_pathway / n_total)
    Returns 0 if the pathway contains all entity references (i.e. it is
    as uninformative as the whole universe / root).
    """
    nb_er_pathway = len(dico_er_per_pathway[pathway])
    if nb_er_pathway != nb_er_total:
        return -math.log10(nb_er_pathway / nb_er_total)
    return 0


# def find_mica_from_several_pathways(
#     graph: nx.Graph,
#     root: Any,
#     list_pathways1: list,
#     list_pathways2: list,
#     dico_er_per_pathway: dict,
#     nb_er_total: int,
# ) -> Any:
#     """
#         Find the Most Informative Common Ancestor (MICA) between two sets of
#         candidate pathways (each entity can have several "most precise
#         parents" when there is a tie in hierarchy depth, hence the lists).

#         For every pair (p1, p2) taken from (list_pathways1, list_pathways2):
#         1. If both lists reduce to the very same single pathway, that
#             pathway is trivially its own MICA.
#         2. Else, if p1 and p2 share at least one common ancestor:
#             - if one is a direct ancestor of the other, that one is the MICA;
#             - otherwise, the common ancestor with the highest IC is the MICA.
#         3. Else (no common ancestor for this pair), fall back to the
#             hierarchy root.

#         The function evaluates every (p1, p2) pair and returns the MICA
#         computed for the last pair examined.
#     """
#     mica = root
#     for p1 in list_pathways1:
#         ancestors1, _ = find_ancestors_and_descendants(graph, p1)
#         for p2 in list_pathways2:
#             ancestors2, _ = find_ancestors_and_descendants(graph, p2)
#             common_ancestors = list(set(ancestors1) & set(ancestors2))

#             if list_pathways1 == list_pathways2 and len(list_pathways1) == 1:
#                 mica = list_pathways1[0]
#             elif common_ancestors:
#                 if p1 in ancestors2:
#                     mica = p1
#                 elif p2 in ancestors1:
#                     mica = p2
#                 else:
#                     max_ic = -1.0
#                     for ancestor in common_ancestors:
#                         ic = compute_ic_er_pathway(ancestor, dico_er_per_pathway, nb_er_total)
#                         if ic > max_ic:
#                             max_ic = ic
#                             mica = ancestor
#             else:
#                 mica = root

#     return mica

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