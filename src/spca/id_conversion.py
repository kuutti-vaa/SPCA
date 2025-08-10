import pandas as pd
import xml.etree.ElementTree as ET
import numpy as np

def parse_hmdb_xml(xml_file, dict_pathway):
    """
    Parse HMDB XML file and store HMDB → ChEBI ID mappings.
    input:
    xml_file: pathway to HMDB xml file which is large but can be downloaded from (https://hmdb.ca/downloads)
    dict_pathway: pathway where to store the hmdb to chebi dict
    
    """
    """
    Efficiently parse HMDB XML using iterparse().
    Extracts {HMDB_ID: ChEBI_ID} mappings while freeing memory.
    """
    hmdb_to_chebi = {}
    context = ET.iterparse(xml_file, events=("start", "end"))
    
    # Get namespace (first element in XML)
    _, root = next(context)

    for event, elem in context:
        if event == "end" and elem.tag.endswith("metabolite"):
            # Extract HMDB ID
            hmdb_id = elem.find("{http://www.hmdb.ca}accession")
            chebi_id = elem.find("{http://www.hmdb.ca}chebi_id")

            if hmdb_id is not None and chebi_id is not None and chebi_id.text:
                hmdb_to_chebi[hmdb_id.text] = chebi_id.text
            
            # Free memory by clearing processed elements
            root.clear()

    # Convert dictionary to DataFrame
    df = pd.DataFrame(hmdb_to_chebi.items(), columns=["HMDB_ID", "ChEBI_ID"])

    # Save as CSV
    df.to_csv(dict_pathway, index=False)

# functions for "clean" reactome. i.e. removes any pathway with a child
def HMDB_to_chebi(results, dict_pathway):
    """ 
    Converts dataset.results object with HMDB annotations to ChEBI annotations. 
    Also removes all non-convertable features (no HMDB ids have ChEBI matches) and annotations which do not match.
    Input:
    results:        A dataset.results dataframe accessed through the Metaspace API
    dict_pathway:   pathway to a .csv file with HMDB to ChEBI id pairs
    returns:
    processed results as ds.results df

    """
    # load in csv and convert to dict
    conv_df = pd.read_csv(dict_pathway, sep = ",")
    hmdb_to_chebi = dict(zip(conv_df["HMDB_ID"], conv_df["ChEBI_ID"] ))

    # as we are loading in HMDB annotations, we have to create a ChEBI id column 
    results["ChEBI_id"] = results["moleculeIds"].apply(lambda x: [hmdb_to_chebi.get(k, np.nan) for k in x]) #convert to string for dict

    results["moleculeIds"] = results["ChEBI_id"]

    # remove all ids without a chebi counterpart (NaN)
    results["moleculeIds"] = results["moleculeIds"].apply(lambda x: [i for i in x if pd.notna(i)])
    results = results[results["moleculeIds"].apply(lambda x: len(x) > 0)]
    #convert to string
    results["moleculeIds"] = results["moleculeIds"].apply(lambda x: list(map(str, x)))
    results.drop(columns=["ChEBI_id"], inplace=True)

    return results