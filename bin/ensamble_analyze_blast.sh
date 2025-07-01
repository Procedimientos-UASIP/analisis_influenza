#!/bin/bash
set -euo pipefail

# Verificar que se proporcione un argumento
if [ "$#" -ne 1 ]; then
    echo "❌ Uso: $0 <archivo_blastn-careful.txt>" >&2
    exit 1
fi

# Guardar argumento en variable y verificar que exista el archivo
INPUT_FILE="$1"
if [ ! -f "$INPUT_FILE" ]; then
    echo "❌ Error: El archivo '$INPUT_FILE' no existe." >&2
    exit 1
fi

# Extraer número de segmento del nombre del directorio actual
BASENAME=$(basename "$INPUT_FILE")
RESULT_FILE="best_result_${BASENAME}"
SEGMENT_NUMBER=$(basename "$PWD" | sed 's/^S//')

# Inicializar variables del mejor hit
best_kmer="" best_node="" best_blast_reference="" best_strand=""
best_length=0 best_slen=0 best_percentage=0 max_cov=0

# Redirigir el archivo al principio del ciclo while
exec 3< "$INPUT_FILE"
    
# Leer el archivo línea por línea
while IFS= read -r line <&3; do

if [[ $line == K-mer* ]]; then
    CURRENT_KMER=$(awk '{print $2}' <<< "$line")

    elif [[ $line == *"segmento_${SEGMENT_NUMBER}"* ]]; then
        # Extraer valores clave
        node=$(awk '{print $1}' <<< "$line")
        cov=$(awk -F"_cov_" '{print $2}' <<< "$line" | awk '{print $1}')
        #length=$(awk '{print $4}' <<< "$line")
        length=$(awk -F'_length_' '{print $2}' <<< "$node" | awk -F'_cov_' '{print $1}')
        slen=$(awk '{print $5}' <<< "$line")
        blast_reference=$(awk -F'|' '{print $2}' <<< "$line")
        #similarity_percentage=$(awk "BEGIN { printf \"%.4f\", ($length / $slen) * 100 }")
        similarity_percentage=$(awk '{print $3}' <<< "$line")
        start=$(awk '{print $9}' <<< "$line")
        end=$(awk '{print $10}' <<< "$line")

        # Determinar la orientación
        if [ "${start}" -lt "${end}" ]; then
            strand="+"
        else
            strand="-"
        fi

        # Seleccionar mejor hit basado en porcentaje y cobertura
        #if awk "BEGIN {exit !($similarity_percentage > $best_percentage || \
        #    ($similarity_percentage == $best_percentage && $cov > $max_cov))}"; then
        if awk -v l="$length" -v bl="$best_length" \
        -v c="$cov" -v bc="$max_cov" \
        -v s="$similarity_percentage" -v bs="$best_percentage" \
            'BEGIN {
                exit !(l > bl || (l == bl && c > bc) || (l == bl && c == bc && s > bs))
            }'; then

            best_kmer="$CURRENT_KMER"
            best_node="$node"
            best_blast_reference="$blast_reference"
            best_length="$length"
            best_slen="$slen"
            best_percentage="$similarity_percentage"
            max_cov="$cov"
            best_strand="$strand"
        fi
    fi
done 

# Extraer número de nodo de mejor match
NODE_NUM=$(echo "$best_node" | grep -oP 'NODE_\K[^_]+(?=_length_)')

# Obtener la descripcion de la referencia de blast
BLAST_DESCRIPTION=$(blastdbcmd \
  -db /backup/DATABASES/NCBI/NT_VIRUSES/nt_viruses \
  -dbtype nucl \
  -entry "${best_blast_reference}" \
  -outfmt "%t" | head -n 1)

# Cerrar el descriptor de archivo cuando ya no se necesite
exec 3<&-

#Guardar resultados en archivo con encabezado y tabulaciones
{
echo -e "SEGMENT\
\tKMER\
\tNODE\
\tLENGTH\
\tSLEN\
\tPERCENTAGE\
\tCOV\
\tSUBJECT_STRAND\
\tBLAST_REFERENCE\
\tBLAST_DESCRIPTION"

echo -e "${SEGMENT_NUMBER}\
\t${KMER}\
\t${NODE_NUM}\
\t${best_length}\
\t${best_slen}\
\t${best_percentage}\
\t${max_cov}\
\t${best_strand}\
\t${best_blast_reference}\
\t${BLAST_DESCRIPTION}"

#echo -e "KMER\tSEGMENT\tLENGTH\tSLEN\tPERCENTAGE\tCOV\tNODE\tBLAST_REFERENCE\tSUBJECT_STRAND"
#echo -e "${best_kmer}\t${SEGMENT_NUMBER}\t${best_length}\t${best_slen}\t${best_percentage}\t${max_cov}\t${best_node}\t${best_blast_reference}\t${best_strand}"
} > "$RESULT_FILE"
echo "✅ Resultados guardados en: $RESULT_FILE"

# Buscar archivo de secuencia correspondiente en KMERS
fasta_file=$(find . -type f -path "*/KMERS/cak${best_kmer}.fasta" | head -n 1)

if [ -z "$fasta_file" ]; then
    echo "❌ No se encontró el archivo FASTA: cak${best_kmer}.fasta en ninguna carpeta KMERS/" >&2
    exit 1
fi

# Extraer la secuencia del mejor node usando seqkit
output_fasta="best_result_S${SEGMENT_NUMBER}_sequence.fna"
seqkit grep -n -p "$best_node" "$fasta_file" > "$output_fasta"

echo "✅ Secuencia extraída en: $output_fasta"

# Corregir la orientación según la variable best_strand
corrected_output="best_result_S${SEGMENT_NUMBER}_sequence_corrected.fna"

if [ "$best_strand" == "-" ]; then
    echo "🔄 Invirtiendo la secuencia (reverso complementario) por orientación '-'..."
    seqkit seq --reverse --complement -t dna -v "$output_fasta" > "$corrected_output"
else
    echo "✅ Manteniendo la secuencia original (orientación '+')..."
    cp "$output_fasta" "$corrected_output"
fi

echo "✅ Secuencia corregida guardada en: $corrected_output"

echo "➤ Calculando profundidad de lecturas de S${SEGMENT_NUMBER} en mejor secuencia ensamblada."

# conda activate alineamiento

# Crear un archivo para calcular la profundidad en cada posición del mejor ensamble obtenido
mkdir PROFUNDIDAD && cd PROFUNDIDAD

# Crear link simbolico al fasta con el mejor ensamble para no copiar archivos 
ln -s "../best_result_S${SEGMENT_NUMBER}_sequence.fna" .

# Crear indice
echo -e "\t➤ Preparando indice"
bwa-mem2 index "best_result_S${SEGMENT_NUMBER}_sequence.fna" > /dev/null 2>&1

# Hacer alineamiento 
echo -e "\t➤ Alineando lecturas pareadas de S${SEGMENT_NUMBER} en mejor secuencia ensamblada."
bwa-mem2 mem -t 20 "best_result_S${SEGMENT_NUMBER}_sequence.fna" \
"../../../../ALINEAMIENTO/BWA/S${SEGMENT_NUMBER}/s${SEGMENT_NUMBER}_reads_r1.fq.gz" \
"../../../../ALINEAMIENTO/BWA/S${SEGMENT_NUMBER}/s${SEGMENT_NUMBER}_reads_r1.fq.gz" \
>"S${SEGMENT_NUMBER}_paired.sam" 2>/dev/null

echo -e "\t➤ Alineando lecturas no pareadas de R1 de S${SEGMENT_NUMBER} en mejor secuencia ensamblada."
bwa-mem2 mem -t 20 "best_result_S${SEGMENT_NUMBER}_sequence.fna" \
"../../../../ALINEAMIENTO/BWA/S${SEGMENT_NUMBER}/s${SEGMENT_NUMBER}_reads_u1.fq.gz" \
>"S${SEGMENT_NUMBER}_unpaired_r1.sam" 2>/dev/null

echo -e "\t➤ Alineando lecturas no pareadas de R2 de S${SEGMENT_NUMBER} en mejor secuencia ensamblada."
bwa-mem2 mem -t 20 "best_result_S${SEGMENT_NUMBER}_sequence.fna" \
"../../../../ALINEAMIENTO/BWA/S${SEGMENT_NUMBER}/s${SEGMENT_NUMBER}_reads_u2.fq.gz" \
>"S${SEGMENT_NUMBER}_unpaired_r2.sam" 2>/dev/null

samtools merge -u - "S${SEGMENT_NUMBER}_paired.sam" "S${SEGMENT_NUMBER}_unpaired_r1.sam" "S${SEGMENT_NUMBER}_unpaired_r2.sam" | samtools sort -o "S${SEGMENT_NUMBER}_all_sorted.bam" && rm "S${SEGMENT_NUMBER}_paired.sam" "S${SEGMENT_NUMBER}_unpaired_r1.sam" "S${SEGMENT_NUMBER}_unpaired_r2.sam"

echo -e "\t➤ Calculando profundidades"
samtools depth -o "S${SEGMENT_NUMBER}_profundidad" -H -a "S${SEGMENT_NUMBER}_all_sorted.bam"

# Regresar al directorio anterior
cd ..
