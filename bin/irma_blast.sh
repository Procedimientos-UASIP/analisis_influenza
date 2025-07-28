#!/usr/bin/env bash
set -uo pipefail

# ============================================
# Configuración
# ============================================
DB="/backup/DATABASES/NCBI/NT_VIRUSES/nt_viruses"
OUTPUT="blast_best_hit_hsp.tsv"
OUTPUT_SORT="blast_best_hit_hsp_sorted.tsv"

# ============================================
# Procesar cada archivo A_*.fasta
# ============================================

# Crear archivo con nombres de columnas
echo "SEGMENT QUERY BITSCORE SUBJECT ALIG_LEN Q_LEN S_LEN PIDENT Q_START Q_END S_START S_END S_TITLE " | tr ' ' '\t' >"$OUTPUT"

for fasta in A_*.fasta
do
    # Extraer nombre del archivo sin el .fasta
    base=$(basename "$fasta" .fasta)

    # Crear archivo temporal
    tmp_file="${base}_raw.tsv"

    echo "Procesando: $fasta"

    # Buscar patrón en $base para obtener el número de segmento
    case "$base" in
        *PB2*) segment="1" ;;
        *PB1*) segment="2" ;;
        *PA*)  segment="3" ;;
        *HA*)  segment="4" ;;
        *NP*)  segment="5" ;;
        *NA*)  segment="6" ;;
        *MP*)  segment="7" ;;
        *NS*)  segment="8" ;;
        *)     echo "⚠️  Segmento no identificado en: $base"; segment="?" ;;
    esac

    # Ejecutar BLASTN 
    blastn -query "$fasta" \
	-db "$DB" \
	-strand 'plus' \
        -outfmt "6 qseqid bitscore saccver length qlen slen pident qstart qend sstart send stitle " \
        -out "$tmp_file" \
	-num_threads 26

    # Obtener mejor subject con mejor mayor bitscore
    best_subject=$( sort -k2,2nr "$tmp_file" | head -n1 | cut -f3 )

    # Extraer todos los HSPs de ese best subject
    awk -v best="$best_subject" -v seg="$segment" '$3 == best { print seg "\t" $0 }' "$tmp_file" >> "$OUTPUT"

    # Eliminar archivo temporal
    rm "$tmp_file"

done

# ============================================
# Ordenar por segmento y por posicion de inicio en subject
# ============================================

# Seleccionar columnas menos bitscore. 
# Ordenar por columna 1 y luego por columna 10 numéricamente. Respetar líneas en blanco.
cut -f1,2,4-13 "$OUTPUT" | LC_ALL=C sort -t $'\t' -k1,1 -k10,10n | awk -F'\t' '
BEGIN { OFS = FS; prev = "" }
{
    if (prev != "" && $1 != prev) {
        print ""  # Línea en blanco entre segmentos
    }
    print
    prev = $1
} ' > "$OUTPUT_SORT"

# Eliminar archivo sin ordenar
rm "$OUTPUT"

