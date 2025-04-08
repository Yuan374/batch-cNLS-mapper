import requests
from bs4 import BeautifulSoup
import csv

# Define the URL of the cNLS Mapper
url = "https://nls-mapper.iab.keio.ac.jp/cgi-bin/NLS_Mapper_y.cgi"

# Function to parse a FASTA file and extract sequence IDs and sequences
def parse_fasta(fasta_file):
    sequences = []
    sequence_id = ""
    sequence = ""
    with open(fasta_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                if sequence:
                    sequences.append((sequence_id, sequence))
                    sequence = ""
                sequence_id = line.split("|")[1]  # Extract the second field as sequence ID
            else:
                sequence += line.strip()
        if sequence:
            sequences.append((sequence_id, sequence))
    return sequences

# Function to submit a single sequence and retrieve the result
def submit_sequence(sequence):
    # Define the form data
    form_data = {
        "typedseq": sequence,  # Input sequence
        "cut_off": "0",  # Cutoff score
        "linker": "Entire region",  # Protein region search
    }
    
    # Send the POST request to the server
    response = requests.post(url, data=form_data)
    
    # Check if the request was successful
    if response.status_code == 200:
        # Parse the response with BeautifulSoup
        soup = BeautifulSoup(response.text, "html.parser")
        
        # Extract mono-partite NLS results
        monopartite_results = []
        table_2 = soup.find_all("table", {"border": "3"})[0]  # Locate the second table
        if table_2:
            rows = table_2.find_all("tr")[2:]  # Skip header rows
            for row in rows:
                cols = row.find_all("td")
                for i in range(len(cols[0].find_all("code"))):
                    position = cols[0].find_all("code")[i].text.strip()
                    sequence = cols[1].find_all("code")[i].text.strip()
                    score = cols[2].find_all("code")[i].text.strip()
                    if position or sequence or score:  # Ensure at least one field is non-empty
                        monopartite_results.append(["Mono-partite", position, sequence, score])
        
        # Extract bi-partite NLS results
        bipartite_results = []
        table_3 = soup.find_all("table", {"border": "3"})[1]  # Locate the third table
        if table_3:
            rows = table_3.find_all("tr")[2:]  # Skip header rows
            for row in rows:
                cols = row.find_all("td")
                for i in range(len(cols[0].find_all("code"))):
                    position = cols[0].find_all("code")[i].text.strip()
                    sequence = cols[1].find_all("code")[i].text.strip()
                    score = cols[2].find_all("code")[i].text.strip()
                    if position or sequence or score:  # Ensure at least one field is non-empty
                        bipartite_results.append(["Bi-partite", position, sequence, score])
        
        return monopartite_results, bipartite_results
    else:
        raise Exception(f"Failed to submit sequence. HTTP status code: {response.status_code}")

# Write results to a CSV file
def write_results_to_csv(filename, results):
    with open(filename, mode="w", newline="") as file:
        writer = csv.writer(file)

        # Write combined results
        writer.writerow(["Sequence ID", "Mono/Bi-partite", "Position", "Sequence", "Score"])
        writer.writerows(results)

# Process multiple sequences
def process_multiple_sequences(fasta_file, output_csv):
    sequences = parse_fasta(fasta_file)
    all_results = []

    for seq_id, (sequence_id, sequence) in enumerate(sequences, start=1):
        print(f"Processing Sequence {seq_id}/{len(sequences)}: {sequence_id}...")
        monopartite_results, bipartite_results = submit_sequence(sequence)

        # Add sequence ID to each result
        for result in monopartite_results:
            all_results.append([sequence_id] + result)
        for result in bipartite_results:
            all_results.append([sequence_id] + result)

    # Write results to CSV
    write_results_to_csv(output_csv, all_results)

# Input FASTA file and output CSV file
fasta_file = "/cwork/yy374/nuclear_localization/idmapping_2024_11_27.fasta"
output_csv = "cNLS_mapper_results.csv"

# Process all sequences and save to CSV
process_multiple_sequences(fasta_file, output_csv)

print(f"Results saved to {output_csv}")

