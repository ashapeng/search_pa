# functions related to data retrieval from opensource databases
import requests
import pandas as pd
from Bio import Entrez, SeqIO, AlignIO
from Bio.Blast import NCBIWWW, NCBIXML

# Retrieve protein data from UniProt
def fetch_protein_data(PROTEIN_DB_API:str,protein_id:str):
    """
    Fetch protein data from UniProt using its API.
    """
    url = f"{PROTEIN_DB_API}{protein_id}.fasta"
    response = requests.get(url)
    # Check if the request was successful
    if response.status_code == 200: # 200 means the request was successful
        return response.text# Return the protein data in FASTA format
    else:
        print(f"Error fetching protein data for {protein_id}")
        return None

# Retrieve transcription data from ENA
def fetch_transcription_data(TRANSCRIPTION_DB_API:str, gene_id:str):
    """
    Fetch transcription data from ENA using its API.
    """
    url = f"{TRANSCRIPTION_DB_API}{gene_id}&display=fasta"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        print(f"Error fetching transcription data for {gene_id}")
        return None

# Find orthologs using OrthoDB
def fetch_orthologs(ORTHODB_API_URL:str, gene_id:str, level="2759"):
    """
    Fetch orthologs for a given gene ID from OrthoDB.
    
    Parameters:
        gene_id (str): The gene ID (e.g., ENSG00000139618).
        level (int): Taxonomic level (default: 2 for vertebrates).
    
    Returns:
        list: A list of orthologs with their details.
    """
    # Define the query parameters
    params = {
        "query": gene_id,
        "level": level
    }
    
    # Make the API request
    response = requests.get(ORTHODB_API_URL, params=params)
    
    # Check if the request was successful
    if response.status_code == 200:
        # Parse the JSON response
        data = response.json()
        
        # Extract orthologs from the response
        if "data" in data:
            return data["data"]
        else:
            print(f"No orthologs found for gene ID: {gene_id}")
            return []
    else:
        print(f"Error fetching orthologs for gene ID: {gene_id}")
        print(response.text)
        return None

# Perform sequence alignment (e.g., BLAST)
def perform_blast(sequence, database="nr"):
    """
    Perform a BLAST search for a given sequence.
    """
    result_handle = NCBIWWW.qblast("blastp", database, sequence)
    blast_record = NCBIXML.read(result_handle)
    return blast_record