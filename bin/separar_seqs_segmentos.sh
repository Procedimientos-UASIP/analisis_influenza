#!/bin/bash

for i in {1..8}
do
    echo "Procesando segmento $i"

    samtools view ./ALINEAMIENTO/BWA/ALL_SEGMENTS_MAPPING/all_sorted.bam | 
    awk -v segmento="segmento_$i" '
        $0 !~ /SA:/ && tolower($0) ~ segmento {
            read_length = length($10);
            cigar = $6;
            sum = 0;
            while (match(cigar, /[0-9]+M/)) {
                sum += substr(cigar, RSTART, RLENGTH - 1);
                cigar = substr(cigar, RSTART + RLENGTH);
            }
            if (sum >= read_length * 0.7) {
                print $1;
            }
        }' |
    sort |
    uniq | sed 's/ //' >"./ALINEAMIENTO/BWA/S${i}/s${i}_lecturas.txt"

    echo "Extrayendo lecturas R1 del segmento $i"

    # Extraer las secuencias R1
    seqkit grep -f "./ALINEAMIENTO/BWA/S${i}/s${i}_lecturas.txt" reads_r1_tr.fq.gz | gzip >"./ALINEAMIENTO/BWA/S${i}/s${i}_reads_r1.fq.gz"

    echo "Extrayendo lecturas R2 del segmento $i"

    # Extraer las secuencias R2
    seqkit grep -f "./ALINEAMIENTO/BWA/S${i}/s${i}_lecturas.txt" reads_r2_tr.fq.gz | gzip >"./ALINEAMIENTO/BWA/S${i}/s${i}_reads_r2.fq.gz"

    echo "Extrayendo lecturas no pareadas de R1 del segmento $i"

    # Extraer secuencias no pareadas de secuencias R1
    seqkit grep -f "./ALINEAMIENTO/BWA/S${i}/s${i}_lecturas.txt" reads_u1_tr.fq.gz | gzip >"./ALINEAMIENTO/BWA/S${i}/s${i}_reads_u1.fq.gz"

    echo "Extrayendo lecturas no pareadas de R2 del segmento $i"

    # Extraer secuencias no pareadas de secuencias R2
    seqkit grep -f "./ALINEAMIENTO/BWA/S${i}/s${i}_lecturas.txt" reads_u2_tr.fq.gz | gzip >"./ALINEAMIENTO/BWA/S${i}/s${i}_reads_u2.fq.gz"

done
