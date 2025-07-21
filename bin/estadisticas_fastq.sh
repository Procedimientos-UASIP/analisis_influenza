#!/bin/bash

# Verificar argumento
if [ $# -ne 1 ]; then
    echo "Uso: $0 archivo.fastq.gz"
    exit 1
fi

FASTQ="$1"
OUTPUT="estadisticas.tsv"

# Verificar que el archivo existe
if [ ! -f "$FASTQ" ]; then
    echo "❌ El archivo '$FASTQ' no existe."
    exit 1
fi

# Ejecutar análisis
zcat "$FASTQ" | awk '
BEGIN {
    total_reads = 0
    total_bases = 0
    total_qual = 0
    min_len = -1
    max_len = -1
    gc_count = 0
}
# Cada lectura tiene 4 líneas: línea 2 = secuencia, línea 4 = calidad
NR % 4 == 2 {
    seq = $0
    len = length(seq)
    total_reads++
    total_bases += len

    # GC count
    gc_count += gsub(/[GCgc]/, "", seq)

    if (min_len == -1 || len < min_len) min_len = len
    if (max_len == -1 || len > max_len) max_len = len
}
NR % 4 == 0 {
    qual = $0
    for (i = 1; i <= length(qual); i++) {
        q = substr(qual, i, 1)
        total_qual += (ord(q) - 33)
    }
}
# Función ord() para obtener el código ASCII
function ord(c) {
    return index(" !\"#$%&'\''()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~", c) + 31
}
END {
    avg_len = (total_reads > 0) ? total_bases / total_reads : 0
    avg_qual = (total_bases > 0) ? total_qual / total_bases : 0
    gc_pct = (total_bases > 0) ? (gc_count / total_bases) * 100 : 0

    printf "No. de lecturas\t%d\n", total_reads
    printf "No. de bases\t%d\n", total_bases
    printf "Tamaño promedio de lecturas\t%.2f\n", avg_len
    printf "Calidad promedio de lecturas\t%.2f\n", avg_qual
    printf "Rango de tamaños de lecturas\t%d - %d\n", min_len, max_len
    printf "%% GC\t%.2f%%\n", gc_pct
}
' > "$OUTPUT"

echo "✅ Estadísticas guardadas en $OUTPUT"