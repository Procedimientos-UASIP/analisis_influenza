#!/bin/bash
set -euo pipefail

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
    # Seleccionar desde el > hacia adelante
    line = substr($0, 2) 

    # Separa el header por | y almacena an un arreglo llamado parts
    n = split(line, parts, "\\|")

    # Limpiar espacios en cada parte
    for (i = 1; i <= n; i++) {
        gsub(/^ +| +$/, "", parts[i])
    }

    # Aislar la descipcion del fasta. 
    desc = parts[2]

    # Eliminar todo el contenido entre paréntesis, incluyendo "A/..." con subtipo (H3N2), etc. para evitar buscar ahí dentro.
    gsub(/\(.*\)/, "", desc)

    # Buscar entonces "segment" fuera del identificador del virus. Si no lo encuentra, lo extrae de otro apartado del encabezado.
    if (match(desc, /segment[^0-9]*([0-9]+)/, seg)) {
        segment_num = seg[1]
    } else {
        segment_num = parts[5]
    }

    # Producir el header final
    printf(">segmento_%s|%s|%s|%s|%s|%s|%s\n", segment_num, parts[1], parts[2], parts[3], parts[4], parts[5], parts[6])
    next
}
{ print }
'