#!/usr/bin/env bash

# Orden definido segun los segmentos de influenza A
SEGMENTS=(A_PB2 A_PB1 A_PA A_HA A_NP A_NA A_MP A_NS)

# Definir nombre de archivo de salida
OUTPUT="prof.tsv"

# Crear directorio temporal
TMP_DIR=$(mktemp -d)

# Iniciar variable para llevar el número máximo de lineas de los archivos de profundidad
MAX_LINES=0

# Iniciar lista de encabezados, basados en el nombre del archivo sin sufijo.
HEADERS=()

# Para cada segmento según el array de nombres
for SEG in "${SEGMENTS[@]}"
do
    # Extraer nombre de ubicación de archivo de profundidades a analizar. 
    # Tiene un comodin para HA y NA
    FILE=$(ls "tables/${SEG}"*-coverage.txt)

    # Extraer el nombre del archivo de la dirección
    FILENAME=$(basename "$FILE")

    # Quitar sujido del nombre del archivo
    PREFIX=${FILENAME%-coverage.txt}

    # Añadir el prefijo al array de encabezados
    HEADERS+=("$PREFIX")

    # Contar líneas del archivo (sin encabezado)
    LINES=$(($(wc -l < "$FILE") - 1))
    if (( LINES > MAX_LINES ))
    then
        MAX_LINES=$LINES
    fi

    # Extraer tercera columna sin encabezado del archivo de profundidades y guardar.
    awk 'NR > 1 {print $3}' "$FILE" > "$TMP_DIR/${PREFIX}.col"
done

# Generar la columna de posiciones y guardarla de forma temporal
seq 1 $((MAX_LINES)) > "$TMP_DIR/position.col"

# Construir la matriz de salida con encabezados
# Primero la linea de encabezados separados por tabs
# Luego la columa de posiciones, y luego todos los archivos de profundidades sin encabezados
{
    printf "Posicion"
    for header in "${HEADERS[@]}"
    do
        printf "\t%s" "$header"
    done
    echo

    paste "$TMP_DIR/position.col" $(printf "$TMP_DIR/%s.col " "${HEADERS[@]}")
} > "$OUTPUT"

# Limpiar archivos temporales
rm -r "$TMP_DIR"