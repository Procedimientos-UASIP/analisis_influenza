#!/bin/bash
set -euo pipefail

# Validar entrada
if [ "$#" -ne 1 ]; then
    echo "❌ Uso: $0 <archivo_fasta>" >&2
    exit 1
fi

FASTA="$1"

if [ ! -f "$FASTA" ]; then
    echo "❌ Error: El archivo '$FASTA' no existe." >&2
    exit 1
fi

# Función para encontrar el motivo más corto
buscar_motivo_mas_corto() {
    local protein="$1"
    local best=""
    local len=${#protein}

    for ((i=0; i<=len-4; i++)); do
        if [[ "${protein:$i:1}" == "P" ]]; then
            for ((j=i+1; j<=len-3; j++)); do
                if [[ "${protein:$j:3}" == "GLF" ]]; then
                    local candidate="${protein:$i:$((j+3-i))}"
                    if [[ -z "$best" || ${#candidate} -lt ${#best} ]]; then
                        best="$candidate"
                    fi
                    break
                fi
            done
        fi
    done

    echo "$best"
}

# Procesar cada secuencia
while read -r header && read -r seq; do
    ID="${header#>}"
    echo "🔍 Análisis de: $ID"

    # Secuencia en archivo temporal
    tmp_fasta=$(mktemp)
    echo -e ">$ID\n$seq" > "$tmp_fasta"

    # Procesar strand +
    for FRAME in 0 1 2; do
        SEQ_FRAME="${seq:$FRAME}"
        LEN=${#SEQ_FRAME}
        MOD=$((3 - LEN % 3))
        if [ "$MOD" -ne 3 ]; then
            SEQ_FRAME="${SEQ_FRAME}$(printf 'N%.0s' $(seq 1 $MOD))"
        fi

        tmp_fasta_frame=$(mktemp)
        echo -e ">$ID\n$SEQ_FRAME" > "$tmp_fasta_frame"

        TRANS=$(seqkit translate -F "$tmp_fasta_frame" | tail -n +2)

        MATCH=$(buscar_motivo_mas_corto "$TRANS")

        if [ -n "$MATCH" ]; then
            POS_COUNT=$(echo "$MATCH" | grep -o '[RHK]' | wc -l || true)
            POS_COUNT=${POS_COUNT:-0}

            if (( POS_COUNT <= 4 )); then
                PRED="LPAI"
            else
                PRED="HPAI"
            fi
            echo -e "Frame +$((FRAME+1))\t$MATCH\t$PRED"
        else
            echo -e "Frame +$((FRAME+1))\t\t"
        fi

        rm -f "$tmp_fasta_frame"
    done

    # Procesar strand -
    REVSEQ=$(seqkit seq --reverse --complement -t dna -v "$tmp_fasta" | tail -n +2)

    for FRAME in 0 1 2; do
        SEQ_FRAME="${REVSEQ:$FRAME}"
        LEN=${#SEQ_FRAME}
        MOD=$((3 - LEN % 3))
        if [ "$MOD" -ne 3 ]; then
            SEQ_FRAME="${SEQ_FRAME}$(printf 'N%.0s' $(seq 1 $MOD))"
        fi

        tmp_fasta_frame=$(mktemp)
        echo -e ">$ID\n$SEQ_FRAME" > "$tmp_fasta_frame"

        TRANS=$(seqkit translate -F "$tmp_fasta_frame" | tail -n +2)

        MATCH=$(buscar_motivo_mas_corto "$TRANS")

        if [ -n "$MATCH" ]; then
            POS_COUNT=$(echo "$MATCH" | grep -o '[RHK]' | wc -l || true)
            POS_COUNT=${POS_COUNT:-0}

            if (( POS_COUNT <= 4 )); then
                PRED="LPAI"
            else
                PRED="HPAI"
            fi
            echo -e "Frame -$((FRAME+1))\t$MATCH\t$PRED"
        else
            echo -e "Frame -$((FRAME+1))\t\t"
        fi

        rm -f "$tmp_fasta_frame"
    done

    rm -f "$tmp_fasta"
done < <(seqkit seq -w 0 "$FASTA")