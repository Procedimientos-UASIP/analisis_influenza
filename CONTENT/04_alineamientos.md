# Alineamiento de lecturas contra base de datos


```bash
mkdir ALINEAMIENTO && cd ALINEAMIENTO 
```

```bash
ln -s ../TRIMMING/[CLEAN_PAIRED_READS_1] reads_r1tr.fq.gz 
ln -s ../TRIMMING/[CLEAN_PAIRED_READS_1] reads_r2tr.fq.gz 
ln -s ../TRIMMING/[CLEAN_UNPAIRED_READS_1] reads_s1tr.fq.gz 
ln -s ../TRIMMING/[CLEAN_UNPAIRED_READS_1] reads_s2tr.fq.gz
```

```bash
mkdir –p BWA/ALL_SEGMENTS_MAPPING
```

```bash
#Alineamiento de secuencias pareadas 
bwa-mem2 mem -t 14 ~/DATABASES/BWAMEM2/INFLUENZA/influenza_nucl_2024-11.fna reads_r1tr.fq.gz reads_r2tr.fq.gz >BWA/ALL_SEGMENTS_MAPPING/paired.sam 

#Alineamiento de secuencias no pareadas 
bwa-mem2 mem -t 14 ~/DATABASES/BWAMEM2/INFLUENZA/influenza_nucl_2024-11.fna reads_rstr.fq.gz >BWA/ALL_SEGMENTS_MAPPING/unpaired.sam
```

```bash
cd BWA/ALL_SEGMENTS_MAPPING 
```

```bash
samtools merge all.sam paired.sam unpaired.sam && samtools view -S -b all.sam >all.bam && samtools sort all.sam >all_sorted.bam && rm all.sam paired.sam unpaired.sam all.bam
```




```bash
for i in {1..8}; do
  echo "Procesando segmento $i"
  samtools view ALL_SEGMENTS_MAPPING/all_sorted.bam |
  grep -v "SA:" | 
  grep "segmento_${i}" | 
  awk '{ print $1"\t"length($10)"\t"$6 }' | 
  awk '{ gsub(/[0-9]M/, "&", $3) }1' | 
  awk '{ gsub(/[0-9][A-LN-Z]/, "", $3) }1'| 
  awk '{ gsub(/M/, "\t", $3) }1' |
  tr " " "\t" | 
  awk '{ for(i=3;i<=NF;i++) t+=$i; { print $1, "\t", $2,"\t", t }; t=0 }' |
  awk '{ if ( $3 >= $2*0.7) print $0 }' |
  cut -f1 |
  sort |
  uniq | 
  sed 's/ //' >S${i}/s${i}_lecturas.txt  
done 
```

