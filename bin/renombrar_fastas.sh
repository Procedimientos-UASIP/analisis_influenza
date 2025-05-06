#!/bin/bash
set -euo pipefail

# Obtener nombre de la muestra desde dos niveles arriba del directorio actual
SAMPLE_NAME=$(basename "$(dirname "$(dirname "$PWD")")")

# Archivo de salida
OUTPUT="${SAMPLE_NAME}_segmentos.fna"
[ -f "$OUTPUT" ] && rm -f "$OUTPUT"
touch "$OUTPUT"

# Función para mapear número de segmento a proteína
get_protein_name() {
    case "$1" in
        1) echo "PB2" ;;
        2) echo "PB1" ;;
        3) echo "PA" ;;
        4) echo "HA" ;;
        5) echo "NP" ;;
        6) echo "NA" ;;
        7) echo "M1_M2" ;;
        8) echo "NS1_NS2" ;;
        *) echo "Desconocido" ;;
    esac
}

# Iterar sobre directorios S1 a S8
for dir in S[1-8]; do
    if [ -d "$dir" ]; then
        numero_segmento=${dir#S}
        protein_name=$(get_protein_name "$numero_segmento")
        fasta_file="${dir}/best_result_${dir}_sequence_corrected.fna"

        if [ -f "$fasta_file" ]; then
            # Extraer header original
            original_header=$(grep '^>' "$fasta_file")

            # Extraer longitud y cobertura con regex
            longitud=$(echo "$original_header" | grep -oP 'length_\K[0-9]+')
            cobertura=$(echo "$original_header" | grep -oP 'cov_\K[0-9.]+')

            # Crear nuevo header
            nuevo_header=">${SAMPLE_NAME}|Segment_${numero_segmento}_${protein_name}|A/HOST/ORIGIN/YEAR(HaNb)|LEN=${longitud}|COV=${cobertura}"

            # Escribir nuevo header y secuencia en una sola línea
            awk -v header="$nuevo_header" '
                BEGIN { ORS="" }
                /^>/ { next }  # Saltar el header original
                {
                    seq = seq $0
                }
                END {
                    print header "\n" seq "\n"
                }
            ' "$fasta_file" >> "$OUTPUT"

        else
            echo "⚠️  Archivo no encontrado: $fasta_file" >&2
        fi
    fi
done

echo "✅ Archivo '$OUTPUT' generado correctamente con secuencias en una sola línea."

