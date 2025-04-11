#!/bin/bash
set -euo pipefail

# Registrar tiempo de inicio del script
START=$SECONDS

############ USAGE #################
usage() {
    >&2 cat << EOF
Usage: $0
   [ -h | --help ] Show this help message
   [ -a | --r1pair ] R1 reads file
   [ -b | --r2pair ] R2 reads file
   [ -x | --unpair ] Unpaired reads file 1
   [ -y | --unpair2] Unpaired reads file 2
   [ -t | --threads ] Threads to use
   [ -o | --outdir ] Output directory
   [ --kini ] Shortest K-mer to use 
   [ --kfin ] Longest K-mer to use
   [ -c | --complete ] Perform full analysis without breaking. Default: FALSE
   [ --kmeroutdir ] Output directory for k-mers
   [ --blastoutput ] BLAST output filename
EOF
    exit 1
}

# Colectar parámetros
args=$(getopt -a -o ha:b:x:y:t:o:c --long help,r1pair:,r2pair:,unpair1:,unpair2:,threads:,outdir:,kini:,kfin:,kmeroutdir:,complete,blastoutput: -- "$@")

# Si no hay parámetros, mostrar ayuda
[[ $? -gt 0 ]] && usage && exit 1

eval set -- "${args}"

while :; do
    case $1 in
        -h | --help)     usage ;;
        -a | --r1pair)   READS_R1=$2   ; shift 2 ;;
        -b | --r2pair)   READS_R2=$2   ; shift 2 ;;
        -x | --unpair1)   READS_U1=$2    ; shift 2 ;;
        -y | --unpair2)   READS_U2=$2    ; shift 2 ;;
        -t | --threads)  THREADS=$2    ; shift 2 ;;
        -o | --outdir)   OUTDIR=$2     ; shift 2 ;;
        --kmeroutdir)    KMER_OUTDIR=$2; shift 2 ;;
        -c | --complete) COMPLETE="TRUE"   ; shift ;;
        --kini)          KINI=$2       ; shift 2 ;;
        --kfin)          KFIN=$2       ; shift 2 ;;
        --blastoutput)   BLAST_OUTPUT=$2; shift 2 ;;
        --) shift; break ;;
        *) >&2 echo "Opción incorrecta: $1"
           usage ;;
    esac
done

# Asignar valores por defecto si no se proporcionan
READS_R1=${READS_R1:-"undefined"}
READS_R2=${READS_R2:-"undefined"}
READS_U1=${READS_U1:-"undefined"}
READS_U2=${READS_U2:-"undefined"}
THREADS=${THREADS:-8}
KINI=${KINI:-11}
KFIN=${KFIN:-127}
OUTDIR=${OUTDIR:-"OUTPUT"}
KMER_OUTDIR=${KMER_OUTDIR:-"${OUTDIR}/KMERS"}
BLAST_OUTPUT=${BLAST_OUTPUT:-"${OUTDIR}/BLAST_RESULTS"}
COMPLETE=${COMPLETE:-"FALSE"}

# Validar que los valores de k-mer sean correctos
if [ $((KINI % 2)) -eq 0 ] || [ "$KINI" -lt 11 ] || [ "$KINI" -gt 127 ]; then
    echo "K-mer de inicio inválido"
    exit 1
fi

if [ $((KFIN % 2)) -eq 0 ] || [ "$KFIN" -lt 11 ] || [ "$KFIN" -gt 127 ]; then
    echo "K-mer de fin inválido"
    exit 1
fi

print_variables() {
    printf "R1: %s\nR2: %s\nU1: %s\nU2: %s\nThreads: %s\nKINI: %s\nKFIN: %s\nOUTDIR: %s\nKMER_OUTDIR: %s\nBLAST_OUTPUT: %s\nCOMPLETE: %s\n" \
        "$READS_R1" "$READS_R2" "$READS_U1" "$READS_U2" "$THREADS" "$KINI" "$KFIN" "$OUTDIR" "$KMER_OUTDIR" "$BLAST_OUTPUT" "$COMPLETE"
}

# Llamada a la función de impresión
print_variables

# Confirmar si los parámetros son correctos antes de continuar
while true; do
    read -rp "¿Los parámetros son correctos? (S/N): " RESPUESTA
    case "${RESPUESTA,,}" in  # convierte a minúsculas para admitir s/S/n/N
        s|si)
            printf "✅ Continuando con el procesamiento...\n\n"
            break
            ;;
        n|no)
            printf "🛑 Ejecución cancelada por el usuario.\n\n"
            exit 0
            ;;
        *)
            printf "❗ Entrada no válida. Por favor responde 'S' (sí) o 'N' (no).\n\n"
            ;;
    esac
done

# Hacer directorios de salida
mkdir "$OUTDIR"
mkdir "$KMER_OUTDIR"
mkdir "$BLAST_OUTPUT"

for KMER in $(seq "$KINI" 2 "$KFIN"); do
    START_FOR=$SECONDS
    echo "➤ [$(date '+%Y-%m-%d %H:%M:%S')] Ensamblando con k-mer = $KMER..."

    # Ejecutar SPAdes de forma silenciosa
    conda run -n ensamble spades.py -o "$OUTDIR" -1 "$READS_R1" -2 "$READS_R2" -s "$READS_U1" -s "$READS_U2" \
        --careful -t "$THREADS" -k "$KMER" >/dev/null 2>&1

    # Verificar salida del ensamblado
    SCAFFOLDS="$OUTDIR/scaffolds.fasta"
    DEST_FILE="$KMER_OUTDIR/cak$KMER.fasta"

    if [[ ! -f "$SCAFFOLDS" ]]; then
        echo "⚠️ [$(date '+%Y-%m-%d %H:%M:%S')] No se encontró '$SCAFFOLDS' para k=$KMER."
    else
        mv "$SCAFFOLDS" "$DEST_FILE"
        echo "✅ [$(date '+%Y-%m-%d %H:%M:%S')] Scaffold encontrado. Se renombra y se mueve a: $DEST_FILE"

        # Realizar BLAST con los primeros 3 scaffolds
        echo "🔍 Realizando BLAST para k-mer = $KMER..."
        {
            echo "K-mer: $KMER"
            conda run -n base seqkit head -n 3 "$DEST_FILE" > temp.fa && \
            conda run -n alineamiento blastn \
                -outfmt "6 qseqid sseqid pident length slen qlen qstart qend sstart send" \
                -max_target_seqs 1 -max_hsps 1 \
                -query temp.fa \
                -db /backup/DATABASES/UASIP/BLAST/INFLUENZA/influenza_db \
                -num_threads "$THREADS" 2>> "${BLAST_OUTPUT}/blastn-warnings.log" && \
            rm temp.fa
            echo ""  # Separador
        } >> "${BLAST_OUTPUT}/blastn-careful.txt"
    fi

    # Calcular tiempo del ciclo
    END_FOR=$SECONDS
    RUN_TIME_FOR=$((END_FOR - START_FOR))

    printf "⏱️ [$(date '+%Y-%m-%d %H:%M:%S')] Proceso para k=%s finalizado en %s min, %s seg.\n\n" \
        "$KMER" "$((RUN_TIME_FOR/60))" "$((RUN_TIME_FOR%60))"
done

# Registrar tiempo de inicio del script
END=$SECONDS

# Calcular segundos que tomó correr el script
RUN_TIME=$((END-START))

printf "\nEl script tardó %s minutos, %s segundos\n\n" "$((RUN_TIME/60))" "$((RUN_TIME%60))"