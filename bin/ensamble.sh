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
   [ --kmeroutdir ] Output directory for k-mers
   [ --blastoutput ] BLAST output filename
EOF
    exit 1
}

# Colectar parámetros
args=$(getopt -o ha:b:x:y:t:o: -l help,r1pair:,r2pair:,unpair1:,unpair2:,threads:,outdir:,kini:,kfin:,kmeroutdir:,blastoutput:,trusted: -- "$@")

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
        --kini)          KINI=$2       ; shift 2 ;;
        --kfin)          KFIN=$2       ; shift 2 ;;
        --blastoutput)   BLAST_OUTPUT=$2; shift 2 ;;
        --trusted)   TRUSTED_CONTIG=$2; shift 2 ;;
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
TRUSTED_CONTIG=${TRUSTED_CONTIG:-""}

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
    printf "R1: %s\nR2: %s\nU1: %s\nU2: %s\nThreads: %s\nKINI: %s\nKFIN: %s\nOUTDIR: %s\nKMER_OUTDIR: %s\nBLAST_OUTPUT: %s\nTRUSTED_CONTIG: %s\n" \
        "$READS_R1" "$READS_R2" "$READS_U1" "$READS_U2" "$THREADS" "$KINI" "$KFIN" "$OUTDIR" "$KMER_OUTDIR" "$BLAST_OUTPUT" "$TRUSTED_CONTIG"
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

# Verificación y creación del directorio
if [ -d "$OUTDIR" ]; then
        echo "❌ Error: El directorio '$OUTDIR' ya existe. No se puede sobrescribir." >&2
        exit 1
    else
        mkdir "$OUTDIR"
        echo "✅ Directorio '$OUTDIR' creado correctamente."
fi

# Creación de directorios internos
mkdir "$KMER_OUTDIR"
mkdir "$BLAST_OUTPUT"


# Para cada uno de los kmeros impares entre el inicio y el final proporcionado
for KMER in $(seq "$KINI" 2 "$KFIN"); do
    # Registrar tiempo de inicio del loop actual
    START_FOR=$SECONDS

    echo "➤  [$(date '+%Y-%m-%d %H:%M:%S')] Ensamblando con k-mer = $KMER..."

    # Si NO se suministro un contig de confianza
    if [[ -z "$TRUSTED_CONTIG" ]]; then
        # TRUSTED_CONTIG está vacía, ejecuta SPAdes sin contigs de confianza
        conda run -n ensamble spades.py -o "$OUTDIR" -1 "$READS_R1" -2 "$READS_R2" \
        -s "$READS_U1" -s "$READS_U2" --careful -t "$THREADS" -k "$KMER" >/dev/null 2>&1
    else
    
    # Si SÍ se suministro un contig de confianza
        if [[ -f "$TRUSTED_CONTIG" ]]; then
            # TRUSTED_CONTIG no está vacía y el archivo existe, usarlo como input confiable
            conda run -n ensamble spades.py -o "$OUTDIR" -1 "$READS_R1" -2 "$READS_R2" \
            -s "$READS_U1" -s "$READS_U2" --careful -t "$THREADS" -k "$KMER" \
            --trusted-contigs "$TRUSTED_CONTIG" >/dev/null 2>&1
        else
            echo "ERROR: El archivo especificado en TRUSTED_CONTIG no existe: $TRUSTED_CONTIG" >&2
            exit 1
        fi
    fi

    # Construir archivos de salida
    SCAFFOLDS="$OUTDIR/scaffolds.fasta"
    DEST_FILE="$KMER_OUTDIR/cak$KMER.fasta"

    # Si no se encuentra un scaffold para el kmero de la iteración
    if [[ ! -f "$SCAFFOLDS" ]]; then
        echo "⚠️ [$(date '+%Y-%m-%d %H:%M:%S')] No se encontró '$SCAFFOLDS' para k=$KMER."
    else
    #Si se produjo un archivo de scaffold para el kmero de la iteracion
        mv "$SCAFFOLDS" "$DEST_FILE"
        echo "✅ [$(date '+%Y-%m-%d %H:%M:%S')] Scaffold encontrado. Se renombra y se mueve a: $DEST_FILE"

        # Realizar BLAST con los primeros 3 scaffolds
        echo "🔍 [$(date '+%Y-%m-%d %H:%M:%S')] Realizando BLAST para k-mer = $KMER..."
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
        } >> "${BLAST_OUTPUT}/blastn-careful.txt"
    fi

    # Registrar tiempo de inicio del loop actual
    END_FOR=$SECONDS

    # Calcular tiempo del ciclo
    RUN_TIME_FOR=$((END_FOR - START_FOR))

    # Imprimir tiempo que tardo el loop actual
    printf "⏱️  [$(date '+%Y-%m-%d %H:%M:%S')] Proceso para k=%s finalizado en %s min, %s seg.\n\n" \
        "$KMER" "$((RUN_TIME_FOR/60))" "$((RUN_TIME_FOR%60))"
done

# Registrar tiempo de inicio del script
END=$SECONDS

# Calcular segundos que tomó correr el script
RUN_TIME=$((END-START))

# Imprimir tiempo el tiempo que tomó correr el script
printf "\nEl script tardó %s minutos, %s segundos\n\n" "$((RUN_TIME/60))" "$((RUN_TIME%60))"