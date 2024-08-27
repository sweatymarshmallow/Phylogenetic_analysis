# Phylogenetic Analysis Script

This script downloads sequences from NCBI, combines and filters them, performs Multiple Sequence Alignment (MSA), and constructs a phylogenetic tree.

## Requirements

- Python 3.x
- Clustal Omega (for MSA, must be installed separately)

## Installation

1. Clone the repository:

    ```bash
    git clone "paste link to clone the repository here"
    cd yourrepository
    ```

2. Install Python dependencies:

    ```bash
    pip install biopython matplotlib
    ```

3. Install Clustal Omega:

    - **On Ubuntu/Debian:**
        ```bash
        sudo apt-get install clustal-omega
        ```

    - **On macOS:**
        ```bash
        brew install clustal-omega
        ```

## Usage

1. Run the script:

    ```bash
    python phylogenetic_analysis.py
    ```

2. Follow the prompts to enter your email, species names, and maximum number of sequences.

The script will create necessary directories and perform the analysis. Results will be saved in the `phylogenetic_analysis` directory.

## Notes

- Ensure your email is correctly entered as it is required by NCBI Entrez.
- Adjust paths and parameters if needed.
