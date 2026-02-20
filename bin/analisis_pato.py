import sys
import os
from Bio import SeqIO
import re
import warnings
from Bio import BiopythonWarning

# Suprimir advertencias de Biopython
warnings.simplefilter('ignore', BiopythonWarning)

def search_motif(seq, sign):
    # Los resultados se depositaran como lista de tuplas
    results = []
    for i in range(3):
        prot_seq = str(seq[i:].translate())
        
        # Establecer patron y buscar coincidencias
        pattern = r"P[^P]+?G[LI]F"
        matches = re.findall(pattern, prot_seq)
        
        # Extraer el match más corto, si no hay se queda como cadena vacía
        shortest = min(matches, key=len) if matches else ""
        
        # Clasificar como HPAI o LPAI basado en el número de aminoácidos básicos
        pato = ""
        if shortest:
            # Contar aminoácidos básicos: Arginina (R), Lisina (K), Histidina (H)
            basic_count = sum(shortest.count(aa) for aa in "RKH")
            pato = "HPAI" if basic_count >= 4 else "LPAI"

        # Etiquetar marco y guardar como tupla
        label = f"Frame {sign}{i+1}"
        results.append((label, shortest, pato))
    return results

def analysis_pato(fasta_file):
    # --- Leer y verificar fasta ---
    if not os.path.exists(fasta_file):
        print(f"\nError: El archivo '{fasta_file}' no existe.\n")
        return

    try:
        # Intentar leer el archivo como FASTA
        record = SeqIO.read(fasta_file, "fasta")
    except ValueError as e:
        print(f"\nError de validación: El archivo no es un FASTA válido o tiene más de una secuencia.\n")
        return
    # ---------------------------

    # Extraer secuencia y generar complementaria
    dna_forward = record.seq
    dna_reverse = record.seq.reverse_complement()

    # Hacer análisis en cada motivo. Usar funcion definida anteriormente
    table_results = search_motif(dna_forward, "+")
    table_results += search_motif(dna_reverse, "-")

    # Tabular resultados
    print("\nAnálisis de patogenicidad basado en el motivo P*G[LI]F")
    print("-" * 50)
    print(f"{'ORF':<8} | {'Patogenicidad':<15} | {'Motivo encontrado'}")
    print("-" * 50)
    for frame, motif, pato in table_results:
        print(f"{frame:<8} | {pato:<15} | {motif}")
    print("-" * 50 + "\n")

if __name__ == "__main__":
    # Verificar que se haya pasado un argumento
    if len(sys.argv) > 1:
        analysis_pato(sys.argv[1])
    else:
        print("Error: Por favor proporciona la ruta del archivo FASTA.")
        print("Uso: python script.py ruta/al/archivo.fa")

