#!/bin/bash
set -euo pipefail

# Validar entrada
if [ "$#" -ne 1 ]; then
    echo "❌ Uso: $0 <archivo_fasta>" >&2
    exit 1
fi

FASTA="$1"

if [ ! -f "$FASTA" ]; then
    echo "❌ Error: El archivo '$FASTA' no existe." >&2
    exit 1
fi

# Extraer nombre de la carpeta actual (e.g., S1, S2, ...)
CARPETA_TMP=$(mktemp -d)
ARCHIVO_SALIDA="Analisis_subtipos_top5.tsv"

# Archivo temporal
BLAST_OUT=$(mktemp)
TMP_TABLE=$(mktemp)
DB="/backup/DATABASES/NCBI/NT_VIRUSES/nt_viruses"

# Separar secuencias individuales del FASTA
echo "📦 Separando secuencias individuales..."
csplit -s -z -f "$CARPETA_TMP/seq_" "$FASTA" '/^>/' '{*}'

# Inicializar archivo final de salida
echo -e "ID\tSubtipo\tConteo\tPorcentaje" > "$ARCHIVO_SALIDA"

# Iterar sobre cada secuencia
for SEQ_FILE in "$CARPETA_TMP"/seq_*; do
    ID=$(grep '^>' "$SEQ_FILE" | sed 's/^>//; s/ .*//')
    echo "🚀 Procesando secuencia: $ID"

    BLAST_OUT=$(mktemp)
    TMP_TABLE=$(mktemp)

    blastn -db "$DB" -query "$SEQ_FILE" -out "$BLAST_OUT" -outfmt "6 stitle" -max_hsps 1 -num_threads 20

    if [ ! -s "$BLAST_OUT" ]; then
        echo "⚠️  No se encontraron hits para $ID." >&2
        continue
    fi

    SUBTIPOS=$(grep -oE 'H[0-9]{1,2}N[0-9]{1,2}' "$BLAST_OUT")
    if [ -z "$SUBTIPOS" ]; then
        echo "⚠️  No se encontraron subtipos HxNy para $ID." >&2
        continue
    fi

    TOTAL=$(echo "$SUBTIPOS" | wc -l)

    echo "$SUBTIPOS" | sort | uniq -c | while read -r COUNT SUBTIPO; do
        PORC=$(awk -v c="$COUNT" -v t="$TOTAL" 'BEGIN { printf "%.4f", (c/t)*100 }')
        printf "%s\t%d\t%s\n" "$SUBTIPO" "$COUNT" "$PORC"
    done | sort -k3,3nr > "$TMP_TABLE"

    # Guardar top 5 + "Otros"
    head -n 5 "$TMP_TABLE" | while IFS=$'\t' read -r SUB COUNT PCT; do
        printf "%s\t%s\t%s\t%.2f\n" "$ID" "$SUB" "$COUNT" "$PCT" >> "$ARCHIVO_SALIDA"
    done

    tail -n +6 "$TMP_TABLE" | awk -v id="$ID" -F'\t' '
    BEGIN { total_count=0; total_pct=0.0 }
    {
        total_count += $2
        total_pct += $3
    }
    END {
        if (total_count > 0)
            printf "%s\t%s\t%d\t%.2f\n\n", id, "Otros", total_count, total_pct
    }' >> "$ARCHIVO_SALIDA"

    rm -f "$BLAST_OUT" "$TMP_TABLE"
done

# Limpieza final
rm -rf "$CARPETA_TMP"

echo "✅ Análisis completo."
