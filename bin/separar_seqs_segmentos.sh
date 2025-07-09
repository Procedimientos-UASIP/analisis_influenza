#!/bin/bash
set -euo pipefail

# Valor por defecto para la identidad
IDENTITY=0.7

# -----------------------------------------
# 1. Parseo de argumentos
# -----------------------------------------
# Usamos getopts para reconocer -i o --identity
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--identity) IDENTITY=$2 shift 2 ;;
        -h|--help)
            echo "Uso: $0 [-i|--identity valor]"
            echo "   -i, --identity    Valor de identidad mínima (0 < x <= 1). Default: 0.7"
            exit 0
            ;;
        *) echo "Argumento desconocido: $1"; exit 1 ;;
    esac
done

# -----------------------------------------
# 2. Validación de rango
# -----------------------------------------
# Validación combinada: formato y rango
if ! [[ "$IDENTITY" =~ ^[0-9]*\.?[0-9]+$ ]]; then
  echo "Error: El valor de identidad debe ser numérico."
  exit 1
fi

if (( $(echo "$IDENTITY < 0" | bc -l) )) || (( $(echo "$IDENTITY > 1" | bc -l) )); then
  echo "Error: El valor de identidad debe estar entre 0 y 1."
  exit 1
fi

# -----------------------------------------
# 3. Ejecutar el filtrado
# -----------------------------------------
for i in {1..8}
do
    echo "Procesando segmento $i"

    samtools view ./ALINEAMIENTO/BWA/ALL_SEGMENTS_MAPPING/all_sorted.bam | 
    awk -v segmento="segmento_$i" -v identity="$IDENTITY" '
        $0 !~ /SA:/ && tolower($0) ~ segmento {
            read_length = length($10);
            cigar = $6;
            sum = 0;
            while (match(cigar, /[0-9]+M/)) {
                sum += substr(cigar, RSTART, RLENGTH - 1);
                cigar = substr(cigar, RSTART + RLENGTH);
            }
            if (sum >= read_length * identity) {
                print $1;
            }
        }' |
    sort |
    uniq | sed 's/ //' >"./ALINEAMIENTO/BWA/S${i}/s${i}_lecturas.txt"

    echo "Extrayendo lecturas R1 del segmento $i"

    # Extraer las secuencias R1
    seqkit grep -f "./ALINEAMIENTO/BWA/S${i}/s${i}_lecturas.txt" reads_r1_tr | gzip >"./ALINEAMIENTO/BWA/S${i}/s${i}_reads_r1.fq.gz"

    echo "Extrayendo lecturas R2 del segmento $i"

    # Extraer las secuencias R2
    seqkit grep -f "./ALINEAMIENTO/BWA/S${i}/s${i}_lecturas.txt" reads_r2_tr | gzip >"./ALINEAMIENTO/BWA/S${i}/s${i}_reads_r2.fq.gz"

    echo "Extrayendo lecturas no pareadas de R1 del segmento $i"

    # Extraer secuencias no pareadas de secuencias R1
    seqkit grep -f "./ALINEAMIENTO/BWA/S${i}/s${i}_lecturas.txt" reads_u1_tr | gzip >"./ALINEAMIENTO/BWA/S${i}/s${i}_reads_u1.fq.gz"

    echo "Extrayendo lecturas no pareadas de R2 del segmento $i"

    # Extraer secuencias no pareadas de secuencias R2
    seqkit grep -f "./ALINEAMIENTO/BWA/S${i}/s${i}_lecturas.txt" reads_u2_tr | gzip >"./ALINEAMIENTO/BWA/S${i}/s${i}_reads_u2.fq.gz"

done
