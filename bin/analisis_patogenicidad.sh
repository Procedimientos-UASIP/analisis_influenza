#!/bin/bash
set -euo pipefail

# Verificar entrada
if [ "$#" -ne 1 ]; then
    echo "❌ Uso: $0 <archivo_fasta>" >&2
    exit 1
fi

FASTA="$1"

if [ ! -f "$FASTA" ]; then
    echo "❌ Error: El archivo '$FASTA' no existe." >&2
    exit 1
fi

# Procesar cada secuencia individualmente
while read -r header && read -r seq; do
    ID="${header#>}"
    echo "🔍 Análisis de: $ID"

    # Crear archivo temporal para la secuencia actual
    tmp_fasta=$(mktemp)
    echo -e ">$ID\n$seq" > "$tmp_fasta"

    # Procesar frames directos
    for FRAME in 0 1 2; do
        # Recortar frame
        SEQ_FRAME="${seq:$FRAME}"
        LEN=${#SEQ_FRAME}
        MOD=$((3 - LEN % 3))
        if [ "$MOD" -ne 3 ]; then
            SEQ_FRAME="${SEQ_FRAME}$(printf 'N%.0s' $(seq 1 $MOD))"
        fi

        # Guardar secuencia recortada temporalmente
        tmp_fasta_frame=$(mktemp)
        echo -e ">$ID\n$SEQ_FRAME" > "$tmp_fasta_frame"

        # Traducir
        TRANS=$(seqkit translate -F "$tmp_fasta_frame" | tail -n +2)

        # Buscar motivo
        MATCH=$(echo "$TRANS" | grep -oP 'P.*?GLF' || true)
        if [ -n "$MATCH" ]; then
            SHORTEST=$(echo "$MATCH" | awk '{ print length, $0 }' | sort -n | head -n1 | cut -d' ' -f2-)
            POS_COUNT=$(echo "$SHORTEST" | grep -o '[RHK]' | wc -l || true)
            if [ -z "$POS_COUNT" ]; then POS_COUNT=0; fi

            if (( POS_COUNT <= 4 )); then
                PRED="LPAI"
            else
                PRED="HPAI"
            fi
            echo -e "Frame +$((FRAME+1))\t$SHORTEST\t$PRED"
        else
            echo -e "Frame +$((FRAME+1))\t\t"
        fi

        rm -f "$tmp_fasta_frame"
    done

    # Reverso complementario
    REVSEQ=$(seqkit seq --reverse --complement -t dna -v "$tmp_fasta" | tail -n +2)

    for FRAME in 0 1 2; do
        SEQ_FRAME="${REVSEQ:$FRAME}"
        LEN=${#SEQ_FRAME}
        MOD=$((3 - LEN % 3))
        if [ "$MOD" -ne 3 ]; then
            SEQ_FRAME="${SEQ_FRAME}$(printf 'N%.0s' $(seq 1 $MOD))"
        fi

        # Guardar secuencia recortada temporalmente
        tmp_fasta_frame=$(mktemp)
        echo -e ">$ID\n$SEQ_FRAME" > "$tmp_fasta_frame"

        # Traducir
        TRANS=$(seqkit translate -F "$tmp_fasta_frame" | tail -n +2)

        # Buscar motivo
        MATCH=$(echo "$TRANS" | grep -oP 'P.*?GLF' || true)
        if [ -n "$MATCH" ]; then
            SHORTEST=$(echo "$MATCH" | awk '{ print length, $0 }' | sort -n | head -n1 | cut -d' ' -f2-)
            POS_COUNT=$(echo "$SHORTEST" | grep -o '[RHK]' | wc -l || true)
            if [ -z "$POS_COUNT" ]; then POS_COUNT=0; fi

            if (( POS_COUNT <= 4 )); then
                PRED="LPAI"
            else
                PRED="HPAI"
            fi
            echo -e "Frame -$((FRAME+1))\t$SHORTEST\t$PRED"
        else
            echo -e "Frame -$((FRAME+1))\t\t"
        fi

        rm -f "$tmp_fasta_frame"
    done

    rm -f "$tmp_fasta"
done < <(seqkit seq -w 0 "$FASTA")