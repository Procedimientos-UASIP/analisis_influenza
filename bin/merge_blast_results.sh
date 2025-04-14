#!/bin/bash

# Recorrer las carpetas S1 a S8
for sample_dir in S*/; do
    sample_name="${sample_dir%/}"
    echo "🔍 Procesando $sample_name..."

    # Archivo de salida para la muestra actual
    output_file="${sample_name}_blastn-careful_merged.txt"
    touch "$output_file"  # Vaciar o crear el archivo

    # Buscar archivos blastn-careful.txt dentro de OUT_*/BLAST_RESULTS/
    for blast_file in "$sample_dir"/OUT_*/BLAST_RESULTS/blastn-careful.txt; do
        if [ -f "$blast_file" ]; then
            echo "📄 Añadiendo: $blast_file"
            cat "$blast_file" >> "$output_file"
            #echo -e "\n" >> "$output_file"  # Separador opcional
        fi
    done

    echo "✅ Archivo combinado creado: $output_file"
done