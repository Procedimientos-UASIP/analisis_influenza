#!/bin/bash
set -euo pipefail

declare -A muestras
valid_keys=("S1" "S2" "S3" "S4" "S5" "S6" "S7" "S8")

# Verificar que se haya proporcionado al menos un argumento
if [[ $# -eq 0 ]]; then
  echo "Uso: $0 [--S1 archivo] [--S2 archivo] ... [--S8 archivo]"
  echo "Debe proporcionar al menos un archivo de cobertura."
  exit 1
fi

# Procesa los argumentos nombrados
while [[ $# -gt 0 ]]; do
  case "$1" in
    --S[1-8])
      key="${1#--}"
      shift
      if [[ $# -eq 0 ]]; then
        echo "Error: No se proporcionó archivo para la opción --$key"
        exit 1
      fi
      file="$1"
      if [[ ! -f "$file" ]]; then
        echo "Error: El archivo '$file' no existe."
        exit 1
      fi
      muestras["$key"]="$file"
      shift
      ;;
    *)
      echo "Error: Opción desconocida '$1'"
      exit 1
      ;;
  esac
done

if [ ${#muestras[@]} -eq 0 ]; then
  echo "Error: No se proporcionó ningún archivo de cobertura."
  exit 1
fi

# Construcción dinámica de la cabecera
header="POS"
orden=()
for key in "${valid_keys[@]}"; do
  if [[ -n "${muestras[$key]:-}" ]]; then
    header+="\t$key"
    orden+=("$key")
  fi
done

# Guardar header en archivo final
echo -e "$header" > profundidades_finales.tsv

# Construcción del comando paste dinámico
awk_commands=()
for key in "${orden[@]}"; do
  archivo="${muestras[$key]}"
  awk_commands+=("<(awk 'NR > 1 { print \$3 }' \"$archivo\")")
done

# Ejecuta el paste y añade columna POS.
eval paste "${awk_commands[*]}" | awk '{ print NR "\t" $0 }' >> prof.tsv