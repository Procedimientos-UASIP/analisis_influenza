import re
import sys
import argparse
import pandas as pd
import warnings
from Bio import SeqIO
from Bio import BiopythonWarning

# Suppress Biopython warnings
warnings.simplefilter('ignore', BiopythonWarning)

def search_motif(seq, sign):
    # The logic of the analysis implies that results must be placed in a list
    results = []
    for i in range(3):
        prot_seq = str(seq[i:].translate())
        
        # Establish pattern and search for matches
        pattern = r"P[^P]+?G[LI]F"
        matches = re.findall(pattern, prot_seq)
        
        # Extract shortest motif or keep result empty
        shortest = min(matches, key=len) if matches else ""
        
        # Predict pathogenicity by count of basic residues (Arg or Lys) in cleavage site
        pato = ""
        if shortest:
            basic_count = sum(shortest.count(aa) for aa in "RK")
            pato = "HPAI" if basic_count >= 4 else "LPAI"

        # Add frame and save data as tupple
        label = f"Frame {sign}{i+1}"
        results.append((label, shortest, pato))
    return results

def main():
    # Config parser
    parser = argparse.ArgumentParser(
        description="Cleavage site analysis for avian influenza pathogenicity."
    )
    
    parser.add_argument("--faa", required=True, help="Path to fasta with HA protein sequence (.faa)")
    parser.add_argument("--fna", required=True, help="Path to fasta with HA nucleotide sequence (.fna)")
    parser.add_argument("--db", required=True, help="Path to csv containing the database of HA cleavage sites.")

    args = parser.parse_args()

    # --- Search in OFFLU cleavage site lists ---
    try:
        seq_reg = SeqIO.read(args.faa, "fasta")
    except Exception as e:
        sys.exit(f"Error reading protein fasta file (--faa): {e}")
        
    seq_protein = str(seq_reg.seq)

    # Load OFFLU data
    df = pd.read_csv(args.db)
    patterns = df["cleavage_site"].unique().tolist()

    # Search for pattern
    pattern_in_list = False
    for pattern in patterns:
        if pattern in seq_protein:
            pattern_in_list = True
            break

    # --- If cleavage site is found in OFFLU list, print data ---
    if pattern_in_list:
        # Filter df according to pattern found
        df_filtrado = df[df["cleavage_site"] == pattern]

        # Header of results
        print(f"\nCleavage site found in OFFLU list: {pattern}")
        print("\n" + "="*40)
        print("DATA REPORT")
        print("="*40)

        # Print each column of the df in a readable format
        for col in df_filtrado.columns:
            data_offlu = []
            for val in df_filtrado[col].tolist():
                # Convert float to int
                if isinstance(val, float) and val.is_integer():
                    data_offlu.append(int(val))
                else:
                    data_offlu.append(val)
            print(f"{col.capitalize()}: {data_offlu}")
        print("="*40)

    # --- If cleavage site is not found in OFFLU list, perform a regex search ---
    if not pattern_in_list:
        print("\nCleavage site not found in OFFLU list.")
        print("\nMaking hard search using regular expressions.")
        try:
            seq_reg = SeqIO.read(args.fna, "fasta")
        except Exception as e:
            sys.exit(f"Error reading nucleotide fasta file (--fna): {e}")
        
        # Extract sequence and reverse complement
        dna_forward = seq_reg.seq
        dna_reverse = seq_reg.seq.reverse_complement()

        # Perform search in each ORF
        table_results = search_motif(dna_forward, "+")
        table_results += search_motif(dna_reverse, "-")

        # Tabulate results
        print("\nPathogenicity analysis based on the pattern: P*G[LI]F")
        print("=" * 50)
        print(f"{'ORF':<8} | {'Pathogenicity':<15} | {'Motif found'}")
        print("=" * 50)
        for frame, motif, pato in table_results:
            print(f"{frame:<8} | {pato:<15} | {motif}")
        print("=" * 50 + "\n")

if __name__ == "__main__":
    main()
