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
samtools merge all.sam paired.sam unpaired.sam && samtools view -S -b all.sam >all.bam && samtools sort all.sam >all_sorted.bam && rm all.sam paired.sam unpaired.sam 
```
