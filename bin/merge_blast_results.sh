#!/bin/bash

# Obtener el nombre de la carpeta actual
sample_name=$(basename "$PWD")
echo "🔍 Procesando $sample_name en $(pwd)"

# Definir el archivo de salida dentro de la carpeta actual
output_file="${sample_name}_blastn-careful_merged.txt"

# Verificar si el archivo ya existe y eliminarlo si es necesario
if [ -f "$output_file" ]; then
    echo -e "\t⚠️  Archivo existente encontrado: $output_file. Será sobreescrito."
    rm "$output_file"
fi

# Crear un nuevo archivo vacío
touch "$output_file"

# Obtener lista ordenada por el primer número de cada carpeta OUT_*
sorted_files=$(ls OUT_*/BLAST_RESULTS/blastn-careful.txt 2>/dev/null | sed -E 's|OUT_([0-9]+)_.*|\1 &|' | sort -n | cut -d' ' -f2-)


# Recorrer y concatenar archivos en el orden deseado
for blast_file in $sorted_files; do
    if [ -f "$blast_file" ]; then
        echo -e "\t📄 Añadiendo: $blast_file"
        cat "$blast_file" >> "$output_file"
    fi
done

echo -e "\t✅ Archivo combinado creado: $output_file"