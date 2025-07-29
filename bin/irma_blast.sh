#!/usr/bin/env bash
set -uo pipefail

# ============================================
# DECLARAR VARIABLES
# ============================================
DB="/backup/DATABASES/NCBI/NT_VIRUSES/nt_viruses"
OUTPUT="blast_best_hit_hsp.tsv"
OUTPUT_SORT="blast_best_hit_hsp_sorted.tsv"

# ============================================
# PROCESAR CADA ARCHIVO A_*.fasta
# ============================================

# Crear archivo con nombres de columnas
echo "SEGMENT QUERY BITSCORE SUBJECT ALIG_LEN Q_LEN S_LEN PIDENT Q_START Q_END S_START S_END S_TITLE " | tr ' ' '\t' >"$OUTPUT"

# Para cada fasta producido por IRMA
for fasta in A_*.fasta
do
    # Extraer nombre del archivo sin la extensión
    NAME=$(basename "$fasta" .fasta)

    # Crear nombre para archivo temporal
    TMP="${NAME}_raw.tsv"

    echo "Procesando: $fasta"

    # Buscar patrón en $NAME para obtener el número de segmento
    case "$NAME" in
        *PB2*) SEGMENT="1" ;;
        *PB1*) SEGMENT="2" ;;
        *PA*)  SEGMENT="3" ;;
        *HA*)  SEGMENT="4" ;;
        *NP*)  SEGMENT="5" ;;
        *NA*)  SEGMENT="6" ;;
        *MP*)  SEGMENT="7" ;;
        *NS*)  SEGMENT="8" ;;
        *)     echo "⚠️  Segmento no identificado en: $NAME"; SEGMENT="?" ;;
    esac

    # Ejecutar BLASTN, con NT_virus, sobre cadena plus. 
    blastn -query "$fasta" \
    -db "$DB" \
    -strand 'plus' \
    -out "$TMP" \
    -outfmt "6 qseqid bitscore saccver length qlen slen pident qstart qend sstart send stitle " \
    -num_threads 26 # Modificar si es requerido

    # Obtener mejor subject en función del mayor bitscore
    best_subject=$( sort -k2,2nr "$TMP" | head -n1 | cut -f3 )

    # Extraer todos los HSPs de ese best subject. Guardar en archivo que ya tienen los nombres de columnas
    awk -v best="$best_subject" -v seg="$SEGMENT" '$3 == best { print seg "\t" $0 }' "$TMP" >> "$OUTPUT"

    # Eliminar archivo temporal
    rm "$TMP"
done

# ============================================
# ORDENAR POR SEGMENTO Y POR POSICIÓN DE INICIO EN SUBJECT
# ============================================

# Guardar encabezado de la tabla
head -n 1 "$OUTPUT" > header.tmp

# Seleccionar columnas, excepto la de bitscore.
# Ordenar por columna 1 y luego por columna 10 numéricamente. Respetar líneas en blanco.
sed -n '2,$p' "$OUTPUT" | cut -f1,2,4-13 | LC_ALL=C sort -t $'\t' -k1,1 -k10,10n | awk -F'\t' '
BEGIN { OFS = FS; prev = "" }
{
    if (prev != "" && $1 != prev) {
        print ""  # Línea en blanco entre segmentos
    }
    print
    prev = $1
} ' > body.tmp

# Unir encabezado a la tabla ordenada
cat header.tmp body.tmp > "$OUTPUT_SORT"

# Eliminar archivo sin ordenar y temporales
rm "$OUTPUT" header.tmp body.tmp

