# Alineamiento de lecturas contra base de datos y ensamble *de novo*
El ensamble genómico es un procedimiento en biología computacional que permite la reconstrucción de secuencias genómicas, a partir de fragmentos más pequeños (llamados "lecturas" o "reads") obtenidos por secuenciación masiva, por ejemplo, utilizando el equipo MySeq de Ilumina. Dichas lecturas están en formato "fastq" que contienen información de la secuencia y la calidad de detección de cada base de cada lectura. Una vez que las lecturas se revisan y se limpian, se alinean entonces contra una base de datos de referencia, se separan las lecturas por segmento viral con el cual hicieron *match*, y se ensamblen los fragmentos virales *de novo*.

Los códigos depositados a continuación son una resumen práctico del "PROCEDIMIENTO PARA EL ANALISIS DE DATOS DE INFLUENZA AVIAR POR SECUENCIACIÓN MASIVA MEIDANTE LA TÉCNICA DE ILLUMINA (MAPEO, ENSAMBLE GENÓMICO Y METRICAS DEL ENSAMBLE)", con código CPA-PD-#####, desarrollado e implementado por la UASIP-CENAPA-SENASICA.

## Código
Se deberá estar posicionado en la carpeta de trabajo de la muestra, ej. */home/CPA-12345-24/*. Los nombres de los archivos en las siguientes instrucciones deben ajustarse a los requerimientos de los archivos a usar. Se indica con corchetes los argumentos que deben/pueden modificarse. A partir de esta estructura de archivos, se procesan los archivos de la siguiente manera:

### Alineamiento
1. Hacer directorio para los alineamientos y ubicarse en dicho directorio
```bash
mkdir ALINEAMIENTO && cd ALINEAMIENTO 
```

2. Crear ligas simbolicas a los archivos a utilizar
```bash
ln -s ../TRIMMING/[CLEAN_PAIRED_READS_1] reads_r1tr.fq.gz 
ln -s ../TRIMMING/[CLEAN_PAIRED_READS_1] reads_r2tr.fq.gz 
ln -s ../TRIMMING/[CLEAN_UNPAIRED_READS_1] reads_s1tr.fq.gz 
ln -s ../TRIMMING/[CLEAN_UNPAIRED_READS_1] reads_s2tr.fq.gz
ln -s ../TRIMMING/reads_rstr.fq.gz
```

3. Hacer directorio para ensamblar TODAS LAS LECTURAS contra la base de datos viral de referencia.
```bash
mkdir –p BWA/ALL_SEGMENTS_MAPPING
```

4. Hacer los alineamientos de TODAS las secuencias pareadas y no pareadas.
```bash
#Alineamiento de secuencias pareadas 
bwa-mem2 mem -t 14 ~/DATABASES/BWAMEM2/INFLUENZA/influenza_nucl_2024-11.fna reads_r1tr.fq.gz reads_r2tr.fq.gz >BWA/ALL_SEGMENTS_MAPPING/paired.sam 

#Alineamiento de secuencias no pareadas 
bwa-mem2 mem -t 14 ~/DATABASES/BWAMEM2/INFLUENZA/influenza_nucl_2024-11.fna reads_rstr.fq.gz >BWA/ALL_SEGMENTS_MAPPING/unpaired.sam
```

4. Entrar al directorio donde se depositaron los resultados del alineamiento.
```bash
cd BWA/ALL_SEGMENTS_MAPPING 
```

5. Formatear los archivos sam resultantes en un archivo bam ordenado.
```bash
samtools merge all.sam paired.sam unpaired.sam && samtools view -S -b all.sam >all.bam && samtools sort all.bam >all_sorted.bam && rm all.sam paired.sam unpaired.sam all.bam
```
6. Regresar al directorio BWA (Revisar ubicación con pwd)
```bash
cd .. 
```
7. Crear ocho directorios para depositar las lecturas por el segmento al cual alinearon.
```bash
mkdir S{1..8}
```

8. Se filtran las lecturas por el segmento al que alinearon. Se incluye un filtro para sólo considerar lecturas con más del 70% de alineamiento de sus bases. Se depositan los resultados en la carpeta correspondiente.

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

9. Regresar al directorio donde están las ligas simbólicas (Revisar ubicación con pwd)
```bash
cd .. 
```

9. Extraer secuencias pareadas y no pareadas para cada segmento.
```bash
for i in {1..8}  
do
  echo "Extrayendo secuencias alineadas del segmento $i"
  #Secuencias forward
  seqkit grep -f ./BWA/S${i}/s${i}_lecturas.txt ./reads_r1tr.fq.gz | gzip >BWA/S${i}/s${i}_reads_r1.fq.gz

  #Secuencias reverse
  seqkit grep -f ./BWA/S${i}/s${i}_lecturas.txt ./reads_r2tr.fq.gz | gzip >BWA/S${i}/s${i}_reads_r2.fq.gz  

  #Secuencias no pareadas
  seqkit grep -f ./BWA/S${i}/s${i}_lecturas.txt ./reads_rstr.fq.gz | gzip >BWA/S${i}/s${i}_reads_rs.fq.gz 
done
```

9. Regresar al directorio nativo del proyecto, donde estan los archivos de lecturas crudas (Revisar ubicación con pwd).
```bash
cd .. 
```

### Ensamble
1. Hacer directorio para los ensambles, y dentro un directorio para los resultados de spades, y ubicarse en directorio de ensamble
```bash
mkdir -p ENSAMBLE/SPADES && cd ENSAMBLE 
```

```bash
ln -s ../TRIMMING/CPA-13120-24_S2_L001_R1_001_val_1.fq.gz reads_r1tr.fq.gz 
ln -s ../TRIMMING/CPA-13120-24_S2_L001_R2_001_val_2.fq.gz reads_r2tr.fq.gz 
ln -s ../TRIMMING/reads_rstr.fq.gz 
```


