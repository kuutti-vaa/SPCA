#chemical class 
from bioservices import ChEBI

def chebi_parent_mol(chebi_id):
    chebi = ChEBI()
    try:
        relation = chebi.getOntologyParents(chebi_id)
        if "ListElement" in relation:
            parents = relation["ListElement"]
            if parents[0][2] == "is a":
                return parents[0][0]
    except Exception as e:
        return None
    
def chebi_parent_id(chebi_id):
    chebi = ChEBI()

    try:
        relation = chebi.getOntologyParents(chebi_id)
        if "ListElement" in relation:
            parents = relation["ListElement"]
            if parents[0][2] == "is a":
                return parents[0][1]
    except Exception as e:
        return None
    
def chebi_parent(chebi_id):
    """
    queries the bioservices ChEBI API to get parent molecule and ChEBI id based off of ChEBI ontology
    Args:

        chebi_id (string): ChEBI id precedeed by "CHEBI:

    Returns:
        parent_info (touple): touple with parent ChEBI id[0] and parent molecule [1]
    """
    chebi = ChEBI()

    try:
        relation = chebi.getOntologyParents(chebi_id)
        if "ListElement" in relation:
            parents = relation["ListElement"]
            if parents[0][2] == "is a":
                return parents[0][1], parents[0][0]
    except Exception as e:
        return None