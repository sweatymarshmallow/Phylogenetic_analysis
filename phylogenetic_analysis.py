import os
import time
import subprocess
from Bio import Entrez, SeqIO, Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt

def get_user_input():
    # Get user inputs
    email = input("Enter your email (required by NCBI): ").strip()
    Entrez.email = email
    
    species_list_input = input("Enter species names separated by commas (or leave empty to use default): ").strip()
    species_list = [species.strip() for species in species_list_input.split(",")] if species_list_input else [
        "Streptomyces albus", "Streptomyces coelicolor", 
        "Streptomyces griseus", "Streptomyces lividans", 
        "Streptomyces scabiei", "Streptomyces avermitilis", 
        "Streptomyces rimosus", "Streptomyces venezuelae", 
        "Streptomyces hygroscopicus", "Streptomyces clavuligerus"
    ]
    
    max_sequences = int(input("Enter the maximum number of sequences to download: ").strip())
    
    return email, species_list, max_sequences

def download_sequences(species_list, max_sequences, output_dir):
    total_downloaded = 0
    for species in species_list:
        if total_downloaded >= max_sequences:
            break
        search_term = f"{species}[Organism]"
        handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=max_sequences-total_downloaded)
        record = Entrez.read(handle)
        handle.close()

        id_list = record["IdList"]
        print(f"Found {len(id_list)} sequences for {species}")

        for seq_id in id_list:
            if total_downloaded >= max_sequences:
                break
            try:
                seq_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
                record = SeqIO.read(seq_handle, "fasta")
                seq_handle.close()
                
                output_file = os.path.join(output_dir, f"{species.replace(' ', '_')}_{seq_id}.fasta")
                with open(output_file, "w") as f:
                    SeqIO.write(record, f, "fasta")
                
                print(f"Downloaded {seq_id} for {species}")
                total_downloaded += 1
            except Exception as e:
                print(f"Error downloading sequence with ID {seq_id}: {e}")
            time.sleep(0.5)  # Be nice to NCBI servers

    print(f"Total sequences downloaded: {total_downloaded}")

def combine_and_filter_fasta_files(input_folder, output_file, min_length=50, max_length=5000):
    with open(output_file, 'w') as outfile:
        for filename in os.listdir(input_folder):
            if filename.endswith('.fasta'):
                with open(os.path.join(input_folder, filename)) as infile:
                    for record in SeqIO.parse(infile, 'fasta'):
                        if min_length <= len(record.seq) <= max_length:
                            SeqIO.write(record, outfile, 'fasta')

def verify_combined_fasta(output_file):
    sequences = list(SeqIO.parse(output_file, "fasta"))
    print(f"Number of sequences in combined file: {len(sequences)}")
    for seq_record in sequences:
        print(f"ID: {seq_record.id}, Length: {len(seq_record.seq)}")
    return len(sequences)

def perform_msa_and_tree_construction(fasta_file, dirs):
    # Perform MSA using Clustal Omega
    clustalomega_command = f"clustalo -i {fasta_file} -o {os.path.join(dirs['alignment'], 'aligned_sequences.fasta')} --force"
    subprocess.run(clustalomega_command, shell=True, check=True)

    # Load the alignment file
    alignment = AlignIO.read(os.path.join(dirs["alignment"], "aligned_sequences.fasta"), "fasta")

    # Calculate distances and create tree
    calculator = DistanceCalculator('identity')
    distances = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor(calculator)
    tree = constructor.build_tree(alignment)

    # Save the tree in Newick format
    tree_file = "phylogenetic_tree.streptomyces.newick"
    Phylo.write(tree, os.path.join(dirs["trees"], tree_file), "newick")

    # Draw and save the tree
    fig = plt.figure(figsize=(100, 100))  # Adjust size as needed
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=ax, show_confidence=False)
    fig.savefig(os.path.join(dirs["trees"], "phylogenetic_tree.png"), dpi=300)
    plt.close(fig)

def main():
    # Get user input
    email, species_list, max_sequences = get_user_input()
    
    # Output directories
    base_dir = "phylogenetic_analysis"
    dirs = {
        "sequences": os.path.join(base_dir, "sequences"),
        "alignment": os.path.join(base_dir, "alignment"),
        "trees": os.path.join(base_dir, "trees")
    }

    # Create directories if they do not exist
    for dir_path in dirs.values():
        os.makedirs(dir_path, exist_ok=True)

    # Run the functions
    download_sequences(species_list, max_sequences, dirs["sequences"])

    combined_fasta_file = os.path.join(base_dir, 'combined_sequences_filtered.fasta')
    combine_and_filter_fasta_files(dirs["sequences"], combined_fasta_file)

    # Verify the combined FASTA file
    num_sequences = verify_combined_fasta(combined_fasta_file)

    # Perform MSA and construct phylogenetic tree
    perform_msa_and_tree_construction(combined_fasta_file, dirs)

if __name__ == "__main__":
    main()
