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
CARPETA_ACTUAL=$(basename "$PWD")
CARPETA_TMP=$(mktemp -d)
ARCHIVO_SALIDA="${CARPETA_ACTUAL}_subtipos_top5.tsv"


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
            printf "%s\t%s\t%d\t%.2f\n", id, "Otros", total_count, total_pct
    }' >> "$ARCHIVO_SALIDA"

    rm -f "$BLAST_OUT" "$TMP_TABLE"
done

# Limpieza final
rm -rf "$CARPETA_TMP"

echo "✅ Análisis completo. Resultados guardados en: $ARCHIVO_SALIDA"

exit


# Archivo temporal
BLAST_OUT=$(mktemp)
TMP_TABLE=$(mktemp)
DB="/backup/DATABASES/NCBI/NT_VIRUSES/nt_viruses"

# Ejecutar BLASTN
echo "🚀 Ejecutando BLASTN..."
blastn -db "$DB" -query "$FASTA" -out "$BLAST_OUT" -outfmt "6 stitle" -max_hsps 1 -num_threads 20

# Verificar que haya resultados
if [ ! -s "$BLAST_OUT" ]; then
    echo "⚠️  No se encontraron hits en BLAST." >&2
    exit 1
fi

# Extraer subtipos HxNy (respetando uno o dos dígitos)
echo "🔍 Extrayendo subtipos HxNy..."
SUBTIPOS=$(grep -oE 'H[0-9]{1,2}N[0-9]{1,2}' "$BLAST_OUT")

if [ -z "$SUBTIPOS" ]; then
    echo "⚠️  No se encontraron subtipos HxNy en los resultados de BLAST." >&2
    exit 1
fi

# Conteo total
TOTAL=$(echo "$SUBTIPOS" | wc -l)

# Preparar archivo intermedio
TMP_TABLE=$(mktemp)

# Calcular frecuencia por subtipo
echo "$SUBTIPOS" | sort | uniq -c | while read -r COUNT SUBTIPO; do
    PORC=$(awk -v c="$COUNT" -v t="$TOTAL" 'BEGIN { printf "%.4f", (c/t)*100 }')
    printf "%s\t%d\t%s\n" "$SUBTIPO" "$COUNT" "$PORC"
done | sort -k3,3nr > "$TMP_TABLE"

# Mostrar encabezado
printf "%-10s %-10s %-10s\n" "Subtipo" "Conteo" "Porcentaje"

# Mostrar top 5
head -n 5 "$TMP_TABLE" | while IFS=$'\t' read -r SUB COUNT PCT; do
    printf "%-10s %-10s %-10.2f\n" "$SUB" "$COUNT" "$PCT"
done

# Calcular "Otros"
tail -n +6 "$TMP_TABLE" | awk -F'\t' '
BEGIN { total_count=0; total_pct=0.0 }
{
    total_count += $2
    total_pct += $3
}
END {
    if (total_count > 0)
        printf "%-10s %-10d %-10.2f\n", "Otros", total_count, total_pct
}
'

# Guardar en archivo de salida
echo -e "Subtipo\tConteo\tPorcentaje" > "$ARCHIVO_SALIDA"
(head -n 5 "$TMP_TABLE"; tail -n +6 "$TMP_TABLE" | awk -F'\t' '
BEGIN { total_count=0; total_pct=0.0 }
{
    total_count += $2
    total_pct += $3
}
END {
    if (total_count > 0)
        printf "Otros\t%d\t%.2f\n", total_count, total_pct
}') >> "$ARCHIVO_SALIDA"

# Limpieza
rm -f "$BLAST_OUT" "$TMP_TABLE"