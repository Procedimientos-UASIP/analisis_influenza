#!/usr/bin/env bash
set -euo pipefail

# ===============================
# PARSEAR ARGUMENTOS
# ===============================

while [[ $# -gt 0 ]]
do
  case $1 in
    --table_blast) TABLE_BLAST="$2"; shift 2 ;;
    *) echo "❌ Argumento desconocido: $1"; exit 1 ;;
  esac
done

# ============================================
# PREPARAR FASTA FINAL
# ============================================

# Iniciar archivo de salida
true >final.fasta

for SEG in {1..8}
do
    # Obtener las lineas que empiezan con el segmento n, de la tabla de blast de los resultados de IRMA
    SEG_LINES=$(grep -P "^$SEG" "$TABLE_BLAST")

    # Obtener la identidad del query (A_XXX_XX)
    QUERY=$(printf "%s\n" "$SEG_LINES" | cut -f2 | uniq)
    
    # Contar las lineas del filtrado
    COUNT=$(printf "%s\n" "$SEG_LINES" | wc -l)

    # Obtener el nombre del fasta que almacena la secuencia a procesar
    FASTA="${QUERY}.fasta"

    # Si sólo hay un match por subject, imprimir el header y la secuencia.
    if [ "$COUNT" -eq 1 ]
    then
        echo ">$QUERY"
        sed -n '2p' "$FASTA"
    # Si hay mpas de un match por subject, se debe agregar N's. Para eso:
    elif [ "$COUNT" -gt 1 ]
    then
        # Iniciar una varialble con la posicion final (0 si es la primera iteración)
        prev_S_END=0
        
        # Bandera
        FLAG_HEADER=0

        # Obtener el nombre del fasta que almacena la secuencia a procesar
        FASTA="${QUERY}.fasta"
        
        # Obtener la secuencia integra a procesar
        SEQUENCE=$(sed -n '2p' "$FASTA")
        
        # Imprimir todas las lineas de la variable del filtrado (como si leyera un archivo)
        # Se guarda el valor de cada columna de cada línea en una variable
        # Se silencian variables que inician con _
        printf "%s\n" "$SEG_LINES" | \
        while IFS=$'\t' read -r SEG QUERY _SUBJECT _ALIG_LEN _Q_LEN _S_LEN _PIDENT Q_START Q_END S_START S_END _S_TITLE
        do #Por cada linea (match)
            # Imprimir encabezado FASTA solo una vez
            if [ "$FLAG_HEADER" -eq 0 ]; then
                echo ">$QUERY"
                FLAG_HEADER=1
            fi
            
            # Calcular índice de inicio del match
            START_INDEX=$(( Q_START - 1))

            # Calcular longitud del match
            LENGTH=$(( Q_END - Q_START + 1 ))

            # Extraer secuencia del match
            SUBSEQ=${SEQUENCE:$START_INDEX:$LENGTH}

            # Calcular gap si no es la primera iteración
            if [ "$prev_S_END" -ne 0 ]
            then
                # Tamaño del gap es la diferencia corregida entre el inicio de un registro y el final del anterior
                GAP_SIZE=$(( S_START - prev_S_END - 1 ))
                if [ "$GAP_SIZE" -gt 0 ]; then
                    # Obtener el string con la cantidad de N según la cantidad de argumentos se pasen a printf, que corresponde a la cantidad de N's
                    GAP_NS=$(printf 'N%.0s' $(seq 1 $GAP_SIZE))
                    # Imprimir las N's
                    printf "%s" "$GAP_NS"
                fi
            fi
            
            # Imprimir la subsecuencia
            printf "%s" "$SUBSEQ"
            
            # Actualizar prev_S_END para la siguiente iteración
            prev_S_END=$S_END
        done
        # Imprimir un salto de línea al final para cuadrar el formato del fasta
        printf "\n"
    fi
done
