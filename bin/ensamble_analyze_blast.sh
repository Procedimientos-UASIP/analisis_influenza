#!/bin/bash
set -euo pipefail

analyze_blast() {
    # Extraer número de segmento del nombre del directorio actual
    SEGMENT_NUMBER=$(basename "$PWD" | sed 's/^S//')

    # Inicializar variables del mejor hit
    best_kmer="" best_node="" best_blast_reference=""
    best_length=0 best_slen=0 best_percentage=0 max_cov=0

    # Redirigir el archivo al principio del ciclo while
    exec 3< "$BLAST_OUTPUT"
    
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
    } > result.txt
}



# Declarar array para checkpoints
read -ra KMER_ARRAY <<< "$(seq "$KINI" 2 "$KFIN")"
# Si es un k-mer de checkpoint y no es ejecución completa
    if printf '%s\n' "${KMER_ARRAY[@]}" | grep -q -x "$KMER" && [[ "$COMPLETE" == "FALSE" ]]; then
        echo "🧪 [$(date '+%Y-%m-%d %H:%M:%S')] Analizando resultados del BLAST para k-mer = $KMER..."
        
        analyze_blast

        # Preguntar al usuario si desea continuar
        while true; do
            read -rp "Estás en k=$KMER. ¿Quieres continuar? (si/no): " RESPUESTA
            case "$RESPUESTA" in
                si) break ;;
                no) echo "🛑 Deteniendo el script en k=$KMER."; exit 0 ;;
                *) echo "❗ Respuesta no válida. Por favor responde 'si' o 'no'." ;;
            esac
        done
    fi