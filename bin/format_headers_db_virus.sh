#!/bin/bash

# Verificar que se haya proporcionado un archivo
if [ $# -ne 1 ]; then
    echo "Uso: $0 <input.fasta[.gz]>"
    exit 1
fi

INPUT="$1"

# Elegir el comando correcto: zcat o cat
if [[ "$INPUT" == *.gz ]]; then
    READER="zcat"
else
    READER="cat"
fi

# Procesar con awk
$READER "$INPUT" | awk '
/^>/ {
    line = substr($0, 2)
    n = split(line, parts, "\\|")

    # Limpiar espacios en cada parte
    for (i = 1; i <= n; i++) {
        gsub(/^ +| +$/, "", parts[i])
    }

    match(parts[2], /segment[ ]*[0-9]+/, seg)
    gsub("segment[ ]*", "", seg[0])

    printf(">segmento_%s|%s|%s|%s|%s|%s|%s\n", seg[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6])
    next
}
{ print }
'