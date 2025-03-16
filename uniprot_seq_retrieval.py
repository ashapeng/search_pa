# functions related to data retrieval from opensource databases
import requests
import argparse
import pandas as pd
import os
from Bio import Entrez, SeqIO, AlignIO
from Bio.Blast import NCBIWWW, NCBIXML

# Retrieve protein data from UniProt
def fetch_protein_data(uniprot_id:str):
    """
    Fetch protein data from UniProt using its API.
    """
    # Define the API URL
    PROTEIN_DB_API = "https://www.uniprot.org/uniprot/"
    
    # Make the API request
    response = requests.get(f"{PROTEIN_DB_API}{uniprot_id}.fasta")
    
    # Check if the request was successful
    if response.status_code == 200: # 200 means the request was successful
        fasta_text = response.text# Return the protein data in FASTA format
        lines = fasta_text.strip().split("\n") # Split the text into lines
        sequence = ''.join(lines[1:])
        
        return fasta_text, sequence

    else:
        print(f"Error fetching protein data for {uniprot_id}")
        return None

# Save protein data to a file as FASTA format to a specified directory
def save_protein_data_to_file(fasta_text:str, filename:str):
    """
    Save protein data to a file in FASTA format.
    """
    # check if the filename has the .fasta extension
    if not filename.endswith(".fasta"):
        filename += ".fasta"
        print(f"Filename changed to {filename}")

    with open(filename, "w") as file:
        file.write(fasta_text)
 
def main():
    parser = argparse.ArgumentParser(description="Extract protein sequence from UniProt")
    parser.add_argument("uniprot_id", help="UniProt accession number (e.g., P12345)")
    parser.add_argument("-o", "--output", help="Output directory and file name (default: the current py directory uniprot_id.fasta)")
    
    args = parser.parse_args()
    
    # Set default output filename if not specified
    if not args.output:
        args.output = f"{args.uniprot_id}.fasta"
    
    print(f"Retrieving protein sequence for {args.uniprot_id}...")
    fasta_text, sequence = fetch_protein_data(args.uniprot_id)
    
    if sequence:
        print(f"Retrieved sequence of length {len(sequence)} amino acids")
        save_protein_data_to_file(fasta_text, args.output)
        
        if not (args.output).endswith(".fasta"):
            args.output += ".fasta" 
        print(f"Sequence saved to {args.output}")
        
        # Print first 50 amino acids as preview
        preview_length = min(50, len(sequence))
        print(f"\nSequence preview (first {preview_length} amino acids):")
        print(sequence[:preview_length] + "..." if len(sequence) > preview_length else sequence)
    else:
        print("no sequence extracted")  # Print error message

if __name__ == "__main__":
    main()

# # Retrieve transcription data from ENA
# def fetch_transcription_data(TRANSCRIPTION_DB_API:str, gene_id:str):
#     """
#     Fetch transcription data from ENA using its API.
#     """
#     url = f"{TRANSCRIPTION_DB_API}{gene_id}&display=fasta"
#     response = requests.get(url)
#     if response.status_code == 200:
#         return response.text
#     else:
#         print(f"Error fetching transcription data for {gene_id}")
#         return None

# # Find orthologs using OrthoDB
# def fetch_orthologs(ORTHODB_API_URL:str, gene_id:str, level="2759"):
#     """
#     Fetch orthologs for a given gene ID from OrthoDB.
    
#     Parameters:
#         gene_id (str): The gene ID (e.g., ENSG00000139618).
#         level (int): Taxonomic level (default: 2 for vertebrates).
    
#     Returns:
#         list: A list of orthologs with their details.
#     """
#     # Define the query parameters
#     params = {
#         "query": gene_id,
#         "level": level
#     }
    
#     # Make the API request
#     response = requests.get(ORTHODB_API_URL, params=params)
    
#     # Check if the request was successful
#     if response.status_code == 200:
#         # Parse the JSON response
#         data = response.json()
        
#         # Extract orthologs from the response
#         if "data" in data:
#             return data["data"]
#         else:
#             print(f"No orthologs found for gene ID: {gene_id}")
#             return []
#     else:
#         print(f"Error fetching orthologs for gene ID: {gene_id}")
#         print(response.text)
#         return None

# # Perform sequence alignment (e.g., BLAST)
# def perform_blast(sequence, database="nr"):
#     """
#     Perform a BLAST search for a given sequence.
#     """
#     result_handle = NCBIWWW.qblast("blastp", database, sequence)
#     blast_record = NCBIXML.read(result_handle)
#     return blast_record