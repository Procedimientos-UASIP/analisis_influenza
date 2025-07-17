#!/usr/bin/env bash
set -euo pipefail

# ===============================
# PARSEAR ARGUMENTOS
# ===============================

while [[ $# -gt 0 ]]; do
  case $1 in
    --n_inicio) N_INICIO="$2"; shift 2 ;;
    --n_final) N_FINAL="$2"; shift 2 ;;
    --input) INPUT="$2"; shift 2 ;;
    *) echo "❌ Argumento desconocido: $1"; exit 1 ;;
  esac
done

# ===============================
# DEFAULTS PARA MIN Y N, SI NO SE DECLARAN
# ===============================

N_INICIO="${N_INICIO:-15}"
N_FINAL="${N_FINAL:-15}"

OUTFILE="fasta_cortado.fna"

# ===============================
# PRODUCIR FASTA CORTADO
# ===============================

# Crear archivo de salida o limpiarlo si ya existe
> "$OUTFILE"

sed -n '2,$p' $INPUT | while IFS=$'\t' read -r HEADER START END NUC_SEQ PROT_SEQ 
do
  # Limpiar header
  HEADER_CLEAN=$(echo "$HEADER" | cut -d '|' -f1-3)

  # Calcular longitud
  LENGTH=$(( END - START + 1 ))

  # Calcular inicio con extensión
  NEW_START=$(( START - $N_INICIO ))

  NEW_LENGTH=$(( LENGTH + $N_INICIO + $N_FINAL ))

  # Extraer subsecuencia
  SUBSEQ=$(echo "$NUC_SEQ" | awk -v s="$NEW_START" -v l="$NEW_LENGTH" '{print substr($0, s, l)}')

  # Guardar en archivo único
  echo ">$HEADER_CLEAN" >> "$OUTFILE"
  echo "$SUBSEQ" >> "$OUTFILE"

done