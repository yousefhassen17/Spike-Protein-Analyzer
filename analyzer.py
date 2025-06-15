# analyzer.py

import os
import sys
from Bio import Entrez, SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd

# --- Configuration ---
Entrez.email = "yhassen134@gmail.com" 
ACCESSION_NUMBERS = {
    'Wuhan': 'YP_009724390.1',
    'Delta': 'QWK65231.1',
    'Omicron': 'UFO69279.1'
}

# --- Functions ---

def fetch_sequences(accession_dict):
    print("Fetching sequences from NCBI...")
    sequences = []
    for variant, acc_num in accession_dict.items():
        try:
            with Entrez.efetch(db="protein", id=acc_num, rettype="fasta", retmode="text") as handle:
                record = SeqIO.read(handle, "fasta")
                record.id = f"{variant}|{acc_num}"
                record.description = ""
                sequences.append(record)
                print(f"  > Successfully fetched {variant} variant.")
        except Exception as e:
            print(f"  > ERROR: Could not fetch sequence for {variant} ({acc_num}). Details: {e}")
    return sequences

def analyze_and_export_properties(sequences):
    print("\nAnalyzing protein properties...")
    analysis_results = []
    for record in sequences:
        variant_name = record.id.split('|')[0]
        protein_analysis = ProteinAnalysis(str(record.seq))
        
        properties = {
            'Variant': variant_name,
            'Accession': record.id.split('|')[1],
            'Sequence Length (AA)': len(record.seq),
            'Molecular Weight (Da)': f"{protein_analysis.molecular_weight():.2f}",
            'Isoelectric Point (pI)': f"{protein_analysis.isoelectric_point():.2f}",
            'Aromaticity': f"{protein_analysis.aromaticity():.4f}",
            'Instability Index': f"{protein_analysis.instability_index():.2f}"
        }
        analysis_results.append(properties)
    
    df = pd.DataFrame(analysis_results)
    df.to_csv('spike_protein_analysis.csv', index=False)
    print("\nAnalysis complete. Results saved to 'spike_protein_analysis.csv'")
    return df

def align_sequences(sequences):
    print("\nAttempting multiple sequence alignment with ClustalW...")
    input_fasta = "temp_sequences.fasta"
    SeqIO.write(sequences, input_fasta, "fasta")
    
    try:
        clustalw_executable = "clustalw2"
        if sys.platform == "win32":
            clustalw_executable = "clustalw2.exe"
            
        clustalw_cline = ClustalwCommandline(clustalw_executable, infile=input_fasta)
        stdout, stderr = clustalw_cline()
        
        alignment = AlignIO.read(input_fasta.replace(".fasta", ".aln"), "clustal")
        with open("sequence_alignment.txt", "w") as f:
            f.write(str(alignment))
        print("  > Alignment successful. Results saved to 'sequence_alignment.txt'")

    except FileNotFoundError:
        print("\n  > WARNING: ClustalW not found. Alignment step skipped.")
        print("  > To enable alignment, install ClustalW and add it to your system's PATH.")
    finally:
        for ext in [".fasta", ".aln", ".dnd"]:
            if os.path.exists(input_fasta.replace(".fasta", ext)):
                os.remove(input_fasta.replace(".fasta", ext))

# --- Main Execution Block ---
def main():
    print("--- SARS-CoV-2 Spike Protein Analyzer ---")
    sequences = fetch_sequences(ACCESSION_NUMBERS)
    
    if not sequences:
        print("No sequences were fetched. Exiting program.")
        return

    analysis_df = analyze_and_export_properties(sequences)
    print("\n--- Analysis DataFrame ---")
    print(analysis_df)
    
    align_sequences(sequences)
    
    print("\n--- Project Execution Finished ---")

if __name__ == "__main__":
    main()