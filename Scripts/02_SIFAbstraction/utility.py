from SPARQLWrapper import SPARQLWrapper, TURTLE, JSON, CSV
import subprocess
import time
import os
from requests.utils import requote_uri
from urllib.parse import quote
import re
import rdflib
import pandas as pd

def extract_prefix_mappings(prefixes_string):
    """
    Extract prefix mappings from the SPARQL prefixes 
    string.
    
    Parameters:
    prefixes_string (str): String containing PREFIX declarations
    
    Returns:
    dict: Mapping of full URIs to their prefixes
    """
    # Extract prefix declarations using regex
    prefix_pattern = re.compile(r'PREFIX\s+(\w+):\s*<([^>]+)>', re.IGNORECASE)
    return {uri: prefix for prefix, uri in prefix_pattern.findall(prefixes_string)}

def convert_to_prefixed_uri(uri_string, prefix_mappings):
    """
    Convert a full URI to prefixed format.
    
    Parameters:
    uri_string (str): Full URI string
    prefix_mappings (dict): Mapping of URIs to prefixes
    
    Returns:
    str: URI in prefixed format (e.g., 'reactome:Protein')
    """
    for uri_base, prefix in prefix_mappings.items():
        if uri_string.startswith(uri_base):
            local_part = uri_string[len(uri_base):]
            return f"{prefix}:{local_part}"
    return uri_string  # Return original if no prefix matches

def save_for_cytoscape(sparql, prefixes_string, output_file, format='csv', separator=','):
    """
    Save SPARQL CONSTRUCT results in a format compatible with Cytoscape,
    using prefix notation for URIs.

    Parameters:
    sparql (SPARQLWrapper): Configured SPARQLWrapper instance with query
    prefixes_string (str): String containing PREFIX declarations
    output_file (str): Path to save the output file
    format (str): Output format ('csv' or 'tsv')
    separator (str): Column separator (',' for CSV, '\t' for TSV)
    """
    # Extract prefix mappings
    prefix_mappings = extract_prefix_mappings(prefixes_string)

    # Get the results as an RDF graph
    sparql.setReturnFormat(TURTLE)
    results = sparql.queryAndConvert()

    # Create an RDFlib graph
    g = rdflib.Graph()
    if isinstance(results, bytes):
        g.parse(data=results.decode('utf-8'), format='turtle')
    else:
        g.parse(data=results, format='turtle')

    # Convert triples to a list of dictionaries with prefixed URIs
    triples_data = []
    for s, p, o in g:
        # Convert each URI to prefixed format
        subject = convert_to_prefixed_uri(str(s), prefix_mappings)
        predicate = convert_to_prefixed_uri(str(p), prefix_mappings)
        object_ = convert_to_prefixed_uri(str(o), prefix_mappings)

        triples_data.append({
            'Source': subject,
            'Interaction': predicate,
            'Target': object_
        })

    # Convert to DataFrame for easy CSV/TSV export
    df = pd.DataFrame(triples_data)

    # Only save if there are results
    if not df.empty:
        # Save to file
        if format == 'csv':
            df.to_csv(output_file, index=False, sep=',')
        else:  # tsv
            df.to_csv(output_file, index=False, sep='\t')

        print(f"Saved {len(triples_data)} interactions to {output_file}")
    else:
        print("No results found - file not saved")

    return df


def preview_network_data(df, n=5):
    """
    Preview the network data before importing into Cytoscape.
    
    Parameters:
    df (pandas.DataFrame): DataFrame containing the network data
    n (int): Number of rows to preview
    """
    print(f"\nPreview of network data ({len(df)} total interactions):")
    print(f"\nFirst {n} interactions:")
    print(df.head(n))
    
    # Print some basic network statistics
    unique_nodes = set(df['Source'].unique()) | set(df['Target'].unique())
    print(f"\nNetwork statistics:")
    print(f"Number of unique nodes: {len(unique_nodes)}")
    print(f"Number of interactions: {len(df)}")
    print(f"Unique interaction types:")
    for interaction in sorted(df['Interaction'].unique()):
        print(f"  - {interaction}")