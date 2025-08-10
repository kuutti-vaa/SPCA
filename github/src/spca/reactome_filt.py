import pandas as pd
import networkx as nx
from io import StringIO
import requests
import matplotlib.pyplot as plt


# functions for Reactome DAG and filtering
def drop_prevalent(reactome_df, num  = 5):
    """ 
    returns top 5 (or num) most populated pathways (least None values in metabolites)

    Args:
        reactome_df (pd.DataFrame): dataframe containing pathway ids (index), pathway names and ChEBI ids for moleculas   
        """
    reactome_df["NaN_Count"] = reactome_df.isna().sum(axis=1)
    reactome_df = reactome_df.sort_values(by="NaN_Count").drop(columns=["NaN_Count"])
    top_pathways_to_drop = reactome_df["Pathway_name"].iloc[:num]
    return top_pathways_to_drop

def build_reactome_dag(org_id = False):
    """ 
    Builds reactome digraph  with parent_id and child_id from reactome API
    
    Args:
        org_ig: Organism id (e.g. for human "HSA"). if left blank all organism will be returned
    """

    url = "https://reactome.org/download/current/ReactomePathwaysRelation.txt"
    response = requests.get(url)
    df = pd.read_csv(StringIO(response.text), sep="\t", header=None)
    df.columns = [ "parent_id","child_id"]
    if org_id:
        df = df[df["parent_id"].str.contains(org_id)]
    G = nx.DiGraph()
    for _, row in df.iterrows():
        G.add_edge(row['parent_id'], row['child_id'])

    return G

def get_pathways_at_level(G, n, from_roots=True):
    """
    Returns all pathway IDs that are exactly n steps away 
    from the root(s) or leaf/leaves of the Reactome DAG.
    
    Args:
        G (networkx.DiGraph): The Reactome DAG.
        n (int): Distance level from root or leaf.
        from_roots (bool): Whether to start from root (True) or leaf (False).
        
    Returns:
        set: Set of pathway IDs at level n
    """
    if from_roots:
        start_nodes = [node for node in G.nodes if G.in_degree(node) == 0]
        direction = G
    else:
        start_nodes = [node for node in G.nodes if G.out_degree(node) == 0]
        direction = G.reverse()

    result = set()
    for start in start_nodes:
        lengths = nx.single_source_shortest_path_length(direction, start, cutoff=n)
        result.update([node for node, dist in lengths.items() if dist == n])
    
    return result

def get_descendant_subgraph(G, root_node):
    """
    given a root node (pathway id) returns the descendant subgraph for plotting

    Args:
        G (networkx.DiGraph): The Reactome DAG.
        root_node: 
    """
    # Get all descendants from the root
    descendants = nx.descendants(G, root_node)
    # Include the root itself
    descendants.add(root_node)
    # Create subgraph
    return G.subgraph(descendants).copy()

def draw_subgraph(G_sub):
    """
    visualizes subgraph
    """
    pos = nx.spring_layout(G_sub)  # or use graphviz_layout for better hierarchy
    plt.figure(figsize=(10, 8))
    nx.draw(
        G_sub,
        pos,
        with_labels=True,
        node_color='lightblue',
        node_size=1000,
        arrows=True,
        font_size=10
    )
    plt.title("Subgraph from Root Node")
    plt.show()