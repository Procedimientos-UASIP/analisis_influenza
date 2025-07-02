#!/usr/bin/env bash
set -euo pipefail

#############################################
# FUNCIONES AUXILIARES
#############################################

print_usage() {
    echo ""
    echo "Uso:"
    echo "  predecir_orf.sh --s1 2000 --s2 2000 --s3 1700 --s4 1500 --s5 1500 --s6 1300 --s7 900 --s8 600 --input archivo.fna"
    echo ""
    exit 1
}

#############################################
# PARSEAR ARGUMENTOS
#############################################

# Inicializa valores vacíos
S1=""; S2=""; S3=""; S4=""; S5=""; S6=""; S7=""; S8=""; INPUT=""

while [[ "$#" -gt 0 ]]; do
  case $1 in
    --s1) S1="$2"; shift ;;
    --s2) S2="$2"; shift ;;
    --s3) S3="$2"; shift ;;
    --s4) S4="$2"; shift ;;
    --s5) S5="$2"; shift ;;
    --s6) S6="$2"; shift ;;
    --s7) S7="$2"; shift ;;
    --s8) S8="$2"; shift ;;
    --input) INPUT="$2"; shift ;;
    *) echo "Error: argumento desconocido: $1"; print_usage ;;
  esac
  shift
done

# Validación básica
if [[ -z "$S1" || -z "$S2" || -z "$S3" || -z "$S4" || -z "$S5" || -z "$S6" || -z "$S7" || -z "$S8" || -z "$INPUT" ]]; then
  echo "Error: faltan argumentos obligatorios."
  print_usage
fi

echo "✔️  Argumentos recibidos:"
echo "  Segmento 1: $S1"
echo "  Segmento 2: $S2"
echo "  Segmento 3: $S3"
echo "  Segmento 4: $S4"
echo "  Segmento 5: $S5"
echo "  Segmento 6: $S6"
echo "  Segmento 7: $S7"
echo "  Segmento 8: $S8"
echo "  Archivo input: $INPUT"
echo ""