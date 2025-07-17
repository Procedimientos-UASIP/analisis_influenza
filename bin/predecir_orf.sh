#!/usr/bin/env bash
set -euo pipefail

# ===============================
# FUNCIONES
# ===============================

print_usage() {
    echo -e "\n Uso:\n\
    predecir_orf.sh\n\
    --s1_min [2180] --s2_min [2180] --s3_min [2050] --s4_min [1600]\n\
    --s5_min [1400] --s6_min [1300] --s7_min [660] --s8_min [600]\n\
    --s1_max [N] --s2_max [N] --s3_max [N] --s4_max [N]\n\
    --s5_max [N] --s6_max [N] --s7_max [N] --s8_max [N]\n\
    --n_inicio [15] --n_final [15]\n\
    --input archivo.fna\n\
    \n\
    Los parámetros min son opcionales, y se indica su valor default. Los parámetros max, e input, son obligatorios.\n"
    exit 1
}

# ===============================
# PARSEAR ARGUMENTOS
# ===============================

while [[ $# -gt 0 ]]; do
  case $1 in
    --s1_min) S1_MIN="$2"; shift 2 ;;
    --s2_min) S2_MIN="$2"; shift 2 ;;
    --s3_min) S3_MIN="$2"; shift 2 ;;
    --s4_min) S4_MIN="$2"; shift 2 ;;
    --s5_min) S5_MIN="$2"; shift 2 ;;
    --s6_min) S6_MIN="$2"; shift 2 ;;
    --s7_min) S7_MIN="$2"; shift 2 ;;
    --s8_min) S8_MIN="$2"; shift 2 ;;

    --s1_max) S1_MAX="$2"; shift 2 ;;
    --s2_max) S2_MAX="$2"; shift 2 ;;
    --s3_max) S3_MAX="$2"; shift 2 ;;
    --s4_max) S4_MAX="$2"; shift 2 ;;
    --s5_max) S5_MAX="$2"; shift 2 ;;
    --s6_max) S6_MAX="$2"; shift 2 ;;
    --s7_max) S7_MAX="$2"; shift 2 ;;
    --s8_max) S8_MAX="$2"; shift 2 ;;

    --n_inicio) N_INICIO="$2"; shift 2 ;;
    --n_final) N_FINAL="$2"; shift 2 ;;
    --input) INPUT="$2"; shift 2 ;;
    *) echo "❌ Argumento desconocido: $1"; exit 1 ;;
  esac
done

# ===============================
# DEFAULTS PARA MIN Y N, SI NO SE DECLARAN
# ===============================

S1_MIN="${S1_MIN:-2180}"
S2_MIN="${S2_MIN:-2180}"
S3_MIN="${S3_MIN:-2050}"
S4_MIN="${S4_MIN:-1600}"
S5_MIN="${S5_MIN:-1400}"
S6_MIN="${S6_MIN:-1300}"
S7_MIN="${S7_MIN:-660}"
S8_MIN="${S8_MIN:-600}"

N_INICIO="${N_INICIO:-15}"
N_FINAL="${N_FINAL:-15}"

# ===============================
# VALIDACIÓN DE VARIABLES OBLIGATORIAS
# ===============================

if [[ -z "${S1_MAX:-}" || -z "${S2_MAX:-}" || -z "${S3_MAX:-}" || -z "${S4_MAX:-}" || \
      -z "${S5_MAX:-}" || -z "${S6_MAX:-}" || -z "${S7_MAX:-}" || -z "${S8_MAX:-}" || \
      -z "${INPUT:-}" ]]; then
  echo "❌ Error: faltan argumentos obligatorios."
  print_usage
fi

# ===============================
# MENSAJE DE CONFIRMACIÓN
# ===============================

echo -e "✔️  Argumentos a usar:\n\
  Segmento 1: min=$S1_MIN  max=$S1_MAX\n\
  Segmento 2: min=$S2_MIN  max=$S2_MAX\n\
  Segmento 3: min=$S3_MIN  max=$S3_MAX\n\
  Segmento 4: min=$S4_MIN  max=$S4_MAX\n\
  Segmento 5: min=$S5_MIN  max=$S5_MAX\n\
  Segmento 6: min=$S6_MIN  max=$S6_MAX\n\
  Segmento 7: min=$S7_MIN  max=$S7_MAX\n\
  Segmento 8: min=$S8_MIN  max=$S8_MAX\n\
  Bases al inicio: $N_INICIO  Bases al final: $N_FINAL\n\
  Archivo input: $INPUT\n"

# ===============================
# MAPEAR VALORES A ARRAYS
# ===============================

declare -A SEGMENT_MIN_SIZES
declare -A SEGMENT_MAX_SIZES

SEGMENT_MIN_SIZES["Segment_1"]=$S1_MIN
SEGMENT_MAX_SIZES["Segment_1"]=$S1_MAX

SEGMENT_MIN_SIZES["Segment_2"]=$S2_MIN
SEGMENT_MAX_SIZES["Segment_2"]=$S2_MAX

SEGMENT_MIN_SIZES["Segment_3"]=$S3_MIN
SEGMENT_MAX_SIZES["Segment_3"]=$S3_MAX

SEGMENT_MIN_SIZES["Segment_4"]=$S4_MIN
SEGMENT_MAX_SIZES["Segment_4"]=$S4_MAX

SEGMENT_MIN_SIZES["Segment_5"]=$S5_MIN
SEGMENT_MAX_SIZES["Segment_5"]=$S5_MAX

SEGMENT_MIN_SIZES["Segment_6"]=$S6_MIN
SEGMENT_MAX_SIZES["Segment_6"]=$S6_MAX

SEGMENT_MIN_SIZES["Segment_7"]=$S7_MIN
SEGMENT_MAX_SIZES["Segment_7"]=$S7_MAX

SEGMENT_MIN_SIZES["Segment_8"]=$S8_MIN
SEGMENT_MAX_SIZES["Segment_8"]=$S8_MAX

# ===============================
# DIVIDIR FASTA EN FASTAS DE UNA SOLA SECUENCIA
# ===============================

WORKDIR="ORFS_PER_SEGMENT"

mkdir -p "$WORKDIR"

cd "$WORKDIR"

# Divide el multifasta en archivos individuales
csplit -z -f seg_ "../$INPUT" '/^>/' '{*}' >/dev/null

# Renombra cada archivo por orden
count=1
for file in seg_*; do
  mv "$file" "Segment_${count}.fna"
  ((count++))
done

cd ..

# ===============================
# EJECUTAR getorf POR SEGMENTO
# ===============================

echo "⏳ Ejecutando getorf por segmento (solo ORF más largo)..."
echo -e "Original_Header\tStart\tEnd\tNucleotide_Sequence\tProtein_Sequence" > ./tbl_longest_orf_range.tsv

for segfile in "$WORKDIR"/Segment_*.fna
do
  segname=$(basename "$segfile" .fna)
  minsize="${SEGMENT_MIN_SIZES[$segname]}"
  maxsize="${SEGMENT_MAX_SIZES[$segname]}"

  echo "  ➜ $segname | minsize: $minsize | maxsize: $maxsize"

  OUT_ORFS="out_orfs_$segname"

  # 1) Ejecutar getorf
  conda run -n emboss getorf -sequence "$segfile" -outseq "${WORKDIR}/${OUT_ORFS}" -auto -find 1 -noreverse -minsize "$minsize" -maxsize "$maxsize"

  # 2) Filtrar solo el ORF más largo
  seqkit sort -l -r "${WORKDIR}/${OUT_ORFS}" 2>/dev/null \
  | seqkit head -n 1 \
  | seqkit seq -w 0 \
  > "${WORKDIR}/longest_orf_$segname"

  # ===============================
  # CONTROL: SI NO SE ENCUENTRA ORF, USAR TODA LA SECUENCIA
  # ===============================

  # Si el archivo no existe o está vacio
  if [[ ! -s "${WORKDIR}/longest_orf_$segname" ]]; then
    echo "⚠️  No se encontró ORF en $segname. Se tomará toda la secuencia."

    original_header=$(head -n 1 "$segfile" | sed 's/^>//')

    nuc_seq=$(grep -v '^>' "$segfile" | tr -d '\n')

    seq_length=${#nuc_seq}
    start=1
    end=$seq_length
    prot_seq="NA"

#Si el archivo existe y contienen datos:
  else
    # Extraer header
    original_header=$(head -n 1 "$segfile" | sed 's/^>//')

    # Obtener coordenadas del ORF más largo, que getorfs depositó en el header
    new_header=$(grep '^>' "${WORKDIR}/longest_orf_$segname")
    start=$(echo "$new_header" | sed -E 's/.*\[([0-9]+)[[:space:]]*-[[:space:]]*([0-9]+)\].*/\1/')
    end=$(echo "$new_header" | sed -E 's/.*\[([0-9]+)[[:space:]]*-[[:space:]]*([0-9]+)\].*/\2/')

    # Obtener secuencia de nucleótidos del ensamble
    nuc_seq=$(grep -v '^>' "$segfile" | tr -d '\n')

    # Obtener la secuencia de aminoácidos del ORF
    prot_seq=$(grep -v '^>' "${WORKDIR}/longest_orf_$segname" | tr -d '\n')
  fi

  # Guardar datos en tabla
  echo -e "${original_header}\t${start}\t${end}\t${nuc_seq}\t${prot_seq}" >> ./tbl_longest_orf_range.tsv

done

