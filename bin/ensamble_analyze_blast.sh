#!/bin/bash
set -euo pipefail

# Verificar que se proporcione un argumento
if [ "$#" -ne 1 ]; then
    echo "❌ Uso: $0 <archivo_blastn-careful.txt>" >&2
    exit 1
fi

# Guardar argumento en variable y verificar que exista el archivo
INPUT_FILE="$1"
if [ ! -f "$INPUT_FILE" ]; then
    echo "❌ Error: El archivo '$INPUT_FILE' no existe." >&2
    exit 1
fi

# Extraer número de segmento del nombre del directorio actual
BASENAME=$(basename "$INPUT_FILE")
RESULT_FILE="best_result_${BASENAME}"
SEGMENT_NUMBER=$(basename "$PWD" | sed 's/^S//')

# Inicializar variables del mejor hit
best_kmer="" best_node="" best_blast_reference=""
best_length=0 best_slen=0 best_percentage=0 max_cov=0

# Redirigir el archivo al principio del ciclo while
exec 3< "$INPUT_FILE"
    
# Leer el archivo línea por línea
while IFS= read -r line <&3; do

if [[ $line == K-mer* ]]; then
    CURRENT_KMER=$(awk '{print $2}' <<< "$line")

    elif [[ $line == *"segmento_${SEGMENT_NUMBER}"* ]]; then
        # Extraer valores clave
        cov=$(awk -F"_cov_" '{print $2}' <<< "$line" | awk '{print $1}')
        length=$(awk '{print $4}' <<< "$line")
        slen=$(awk '{print $5}' <<< "$line")
        node=$(awk '{print $1}' <<< "$line")
        blast_reference=$(awk -F'|' '{print $2}' <<< "$line")
        similarity_percentage=$(awk "BEGIN { printf \"%.4f\", ($length / $slen) * 100 }")

        # Seleccionar mejor hit basado en porcentaje y cobertura
        if awk "BEGIN {exit !($similarity_percentage > $best_percentage || \
            ($similarity_percentage == $best_percentage && $cov > $max_cov))}"; then

            best_kmer="$CURRENT_KMER"
            best_node="$node"
            best_blast_reference="$blast_reference"
            best_length="$length"
            best_slen="$slen"
            best_percentage="$similarity_percentage"
            max_cov="$cov"
        fi
    fi
done 
    
# Cerrar el descriptor de archivo cuando ya no se necesite
exec 3<&-

# Mostrar resultados
echo "NODE del hit: $best_node"
echo "Referencia del BLAST: $best_blast_reference"
echo "Cobertura (cov): $max_cov"
echo "K-mer: $best_kmer"
echo "Segmento: $SEGMENT_NUMBER"
echo "Length = $best_length, Slen = $best_slen"
echo "Porcentaje de similitud: $best_percentage%"
    
#Guardar resultados en archivo con encabezado y tabulaciones
{
echo -e "KMER\tSEGMENT\tLENGTH\tSLEN\tPERCENTAGE\tCOV\tNODE\tBLAST_REFERENCE"
echo -e "${best_kmer}\t${SEGMENT_NUMBER}\t${best_length}\t${best_slen}\t${best_percentage}\t${max_cov}\t${best_node}\t${best_blast_reference}"
} > "$RESULT_FILE"
echo "✅ Resultados guardados en: $RESULT_FILE"
