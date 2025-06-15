# SARS-CoV-2 Spike Protein Analyzer

A Python tool for fetching, analyzing, and comparing Spike protein sequences from major SARS-CoV-2 variants.

## Features

- Fetches sequence data directly from the NCBI database.
- Calculates and compares key physicochemical properties (e.g., molecular weight, pI, instability index).
- Exports analysis results to a clean `.csv` file.
- Performs multiple sequence alignment using ClustalW to identify variations.

## Tech Stack

- **Language:** Python 3
- **Core Libraries:** Biopython, Pandas, NumPy
- **Alignment Tool:** ClustalW (external dependency)

## Quickstart

1.  **Clone this repository:**
    ```bash
    git clone [https://github.com/your-username/Spike-Protein-Analyzer.git](https://github.com/your-username/Spike-Protein-Analyzer.git)
    cd Spike-Protein-Analyzer
    ```

2.  **Installing Python libraries:**
    ```bash
    pip install -r requirements.txt 
    # (Assuming you create a requirements.txt with biopython, pandas, numpy)
    # Or manually:
    pip install biopython pandas numpy
    ```

3.  **Running the script:**
    ```bash
    python analyzer.py
    ```
    *(Note: The alignment feature requires ClustalW to be installed and accessible in your system's PATH.)*

## Results Preview

The script generates `spike_protein_analysis.csv` with detailed property comparisons:

| Variant | Accession      | Sequence Length (AA) | Molecular Weight (Da) | Isoelectric Point (pI) |
|---------|----------------|----------------------|-----------------------|------------------------|
| Wuhan   | YP_009724390.1 | 1273                 | 141178.40             | 6.79                   |
| Delta   | QWK65231.1     | 1273                 | 141154.38             | 7.03                   |
| Omicron | UFO69279.1     | 1270                 | 140816.03             | 7.56                   |

---
