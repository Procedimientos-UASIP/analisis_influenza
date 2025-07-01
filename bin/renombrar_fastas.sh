#!/bin/bash
set -euo pipefail

# Inicializar variables con valor por defecto
HOST="undefined"
ORIGIN="undefined"
YEAR="undefined"

# Procesar opciones de línea de comandos
while [[ $# -gt 0 ]]; do
    case "$1" in
        --host)
            if [[ $# -lt 2 || "$2" == --* ]]; then
                echo "❌ Error: La opción --host requiere un argumento." >&2
                exit 1
            fi
            HOST="${2// /_}"  # Reemplazar espacios por guion bajo
            shift 2
            ;;
        --origin)
            if [[ $# -lt 2 || "$2" == --* ]]; then
                echo "❌ Error: La opción --origin requiere un argumento." >&2
                exit 1
            fi
            ORIGIN="${2// /_}"
            shift 2
            ;;
        --year)
            if [[ $# -lt 2 || "$2" == --* ]]; then
                echo "❌ Error: La opción --year requiere un argumento." >&2
                exit 1
            fi
            YEAR="$2"
            if [[ ! "$YEAR" =~ ^[0-9]{4}$ ]]; then
                echo "❌ Error: '--year' debe ser un número de 4 dígitos." >&2
                exit 1
            fi
            shift 2
            ;;
        *)
            echo "❌ Opción no reconocida: $1" >&2
            exit 1
            ;;
    esac
done

# Validar que todas las opciones hayan sido asignadas
if [[ "$HOST" == "undefined" ]]; then
    echo "❌ Error: Falta especificar la opción --host." >&2
    exit 1
fi

if [[ "$ORIGIN" == "undefined" ]]; then
    echo "❌ Error: Falta especificar la opción --origin." >&2
    exit 1
fi

if [[ "$YEAR" == "undefined" ]]; then
    echo "❌ Error: Falta especificar la opción --year." >&2
    exit 1
fi

# Obtener nombre de la muestra desde dos niveles arriba del directorio actual
SAMPLE_NAME=$(basename "$(dirname "$(dirname "$PWD")")")

# Archivo de salida
OUTPUT="${SAMPLE_NAME}_secuencias_renombradas.fna"
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
        subtipo_file="${dir}/S${numero_segmento}_subtipos_top5.tsv"
        
        if [ -f "$fasta_file" ]; then
            # Extraer subtipo principal
            if [ -f "$subtipo_file" ]; then
                subtipo=$(awk 'NR==2 {print $1}' "$subtipo_file")
            else
                echo "⚠️  Archivo de subtipo no encontrado: $subtipo_file" >&2
                subtipo="HaNb"
            fi

            # Extraer header original
            ## original_header=$(grep '^>' "$fasta_file")

            ## Extraer longitud y cobertura con regex
            ## longitud=$(echo "$original_header" | grep -oP 'length_\K[0-9]+')
            ## cobertura=$(echo "$original_header" | grep -oP 'cov_\K[0-9.]+')

            # Crear nuevo header
            ## nuevo_header=">${SAMPLE_NAME}|Segment_${numero_segmento}_${protein_name}|A/${HOST}/${ORIGIN}/${YEAR}(${subtipo})|Length=${longitud}|Coverage=${cobertura}"
            nuevo_header=">${SAMPLE_NAME}|Segment_${numero_segmento}_${protein_name}|A/${HOST}/${ORIGIN}/${YEAR}(${subtipo})"

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

echo "✅ Archivo '$OUTPUT' generado correctamente."