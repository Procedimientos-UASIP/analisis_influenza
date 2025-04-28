#!/bin/bash
set -euo pipefail

# Archivo de salida
SALIDA="uso_de_lecturas.txt"

# Título informativo
echo "📋 Calculando porcentaje de lecturas usadas en ensamblaje de novo para cada segmento (S1-S8)..."

# Escribir encabezado formateado en el archivo
printf "%-10s\t%-20s\t%-20s\t%-20s\n" "Segmento" "Lecturas_alineadas" "Lecturas_ensamblaje" "Porcentaje_usadas" > "$SALIDA"

# También mostrar encabezado en pantalla
printf "%-10s\t%-20s\t%-20s\t%-20s\n" "Segmento" "Lecturas_alineadas" "Lecturas_ensamblaje" "Porcentaje_usadas"


# Bucle para cada carpeta S1 a S8
for i in {1..8}; do
    SEGMENTO="S${i}"

    # Verificar existencia de archivos necesarios
    READS_USADAS_ENSAMBLE="${SEGMENTO}/COBERTURA/${SEGMENTO}_all_sorted.bam"
    PATH_READS_POR_SEGMENTO="../../ALINEAMIENTO/BWA/${SEGMENTO}"

    if [ ! -f "$READS_USADAS_ENSAMBLE" ]; then
        echo "⚠️  No se encontró el archivo BAM para $SEGMENTO, omitiendo."
        continue
    fi

    if [ ! -d "$PATH_READS_POR_SEGMENTO" ]; then
        echo "⚠️  No se encontró la carpeta de lecturas para $SEGMENTO, omitiendo."
        continue
    fi

    # Obtener cantidad de lecturas alineadas
    LECTURAS_ENSAMBLADAS=$(samtools view -F 260 "$READS_USADAS_ENSAMBLE" | cut -f1 | sort | uniq | wc -l)

    # Obtener cantidad de lecturas usadas para ensamblaje de novo
    LECTURAS_ALINEADAS=$(zgrep '^@' "${PATH_READS_POR_SEGMENTO}/s${i}_reads_r1.fq.gz" "${PATH_READS_POR_SEGMENTO}/s${i}_reads_u1.fq.gz" "${PATH_READS_POR_SEGMENTO}/s${i}_reads_u2.fq.gz" | cut -d' ' -f1 | sed 's/^@//' | sort | uniq | wc -l)

    # Calcular porcentaje
    if [ "$LECTURAS_ALINEADAS" -gt 0 ]; then
        PORCENTAJE=$(awk "BEGIN { printf \"%.2f\", ($LECTURAS_ENSAMBLADAS / $LECTURAS_ALINEADAS) * 100 }")
    else
        PORCENTAJE="NA"
    fi

    # Mostrar resultados formateados en pantalla
    printf "%-10s\t%-20s\t%-20s\t%-20s\n" "$SEGMENTO" "$LECTURAS_ALINEADAS" "$LECTURAS_ENSAMBLADAS" "$PORCENTAJE"

    # Guardar resultados formateados en el archivo
    printf "%-10s\t%-20s\t%-20s\t%-20s\n" "$SEGMENTO" "$LECTURAS_ALINEADAS" "$LECTURAS_ENSAMBLADAS" "$PORCENTAJE" >> "$SALIDA"
done

echo "✅ Resultados guardados en '$SALIDA'."