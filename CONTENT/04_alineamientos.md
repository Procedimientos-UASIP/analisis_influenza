# ALINEAMIENTO DE LECTURAS CONTRA BASE DE DATOS DE INFLUENZA

El ensamble genómico es un procedimiento en biología computacional que permite la reconstrucción de secuencias genómicas, a partir de fragmentos más pequeños (llamados "lecturas" o "reads") obtenidos por secuenciación masiva, por ejemplo, utilizando el equipo MySeq de Ilumina. Dichas lecturas están en formato "fastq" que contienen información de la secuencia y la calidad de detección de cada base de cada lectura. Una vez que las lecturas se revisan y se limpian, se alinean entonces contra una base de datos de referencia, se separan las lecturas por segmento viral con el cual hicieron *match*, y se ensamblen los fragmentos virales *de novo*.

Los códigos depositados a continuación son una resumen práctico del "PROCEDIMIENTO PARA EL ANALISIS DE DATOS DE INFLUENZA AVIAR POR SECUENCIACIÓN MASIVA MEDIANTE LA TÉCNICA DE ILLUMINA
 (ALINEAMIENTO Y FILTRADO)", <ins>(CPA-PD-XXXXX)</ins>, desarrollado e implementado por la UASIP-CENAPA-SENASICA.

## Código
Se deberá estar posicionado en la carpeta de trabajo de la muestra, ej. */home/CPA-12345-24/* **(Los nombres de los archivos en las siguientes instrucciones deben ajustarse a los requerimientos de los archivos a usar)**. A partir de esta estructura de archivos, se procesan los archivos de la siguiente manera:

### Preparación
1. Hacer directorio para los alineamientos:
```bash
mkdir -p ALINEAMIENTO/BWA/ALL_SEGMENTS_MAPPING
```

2. Crear ligas simbolicas a los archivos a usar:
```bash
ln -s TRIMMING/[CLEAN_PAIRED_READS_1] reads_r1_tr
ln -s TRIMMING/[CLEAN_PAIRED_READS_1] reads_r2_tr
ln -s TRIMMING/[CLEAN_UNPAIRED_READS_1] reads_u1_tr
ln -s TRIMMING/[CLEAN_UNPAIRED_READS_1] reads_u2_tr
```
### Alineamiento
3. Hacer los alineamientos de las secuencias pareadas y no pareadas contra la base de datos viral preparada en la UASIP.
```bash
#Activar ambiente conda
conda activate alineamiento

#Alineamiento de secuencias pareadas 
bwa-mem2 mem -t 20 /backup/DATABASES/UASIP/BWAMEM2/INFLUENZA/influenza_db.fna reads_r1_tr reads_r2_tr >ALINEAMIENTO/BWA/ALL_SEGMENTS_MAPPING/paired.sam 

#Alineamiento de secuencias no pareadas de R1
bwa-mem2 mem -t 20 /backup/DATABASES/UASIP/BWAMEM2/INFLUENZA/influenza_db.fna reads_u1_tr >ALINEAMIENTO/BWA/ALL_SEGMENTS_MAPPING/unpaired_r1.sam

#Alineamiento de secuencias no pareadas de R2
bwa-mem2 mem -t 20 /backup/DATABASES/UASIP/BWAMEM2/INFLUENZA/influenza_db.fna reads_u2_tr >ALINEAMIENTO/BWA/ALL_SEGMENTS_MAPPING/unpaired_r2.sam
```

4. Entrar al directorio donde se depositaron los resultados del alineamiento.
```bash
cd ALINEAMIENTO/BWA/ALL_SEGMENTS_MAPPING 
```

5. Formatear los archivos sam resultantes en un archivo bam ordenado.
```bash
samtools merge -u - paired.sam unpaired_r1.sam unpaired_r2.sam | samtools sort -o all_sorted.bam && rm paired.sam unpaired_r1.sam unpaired_r2.sam
```

6. Regresar al directorio del proyecto:
```bash
cd ../../../
```

### Filtrado
7. Crear ocho directorios para depositar las lecturas por el segmento al cual alinearon.
```bash
mkdir ALINEAMIENTO/BWA/S{1..8}
```

8. Ejecutar el siguiente programa:
```bash
~/analisis_influenza/bin/separar_seqs_segmentos.sh
```

Este programa abre el archivo bam, elimina alineamientos con alineamientos secundarios, y obtiene las lecturas por cada segmento. Se analiza, por segmento, el formato CIGAR asociado de cada lectura. Se filtran aquellas secuencias que alinearon en al menos 70% de sus bases. Se guardan los nombres de las lecturas pertinentes y se extraen las secuencias de los archivos fastq limpios usando dichos nombres. Los archivos producidos se guardan en las ubicaciones pertinentes.
