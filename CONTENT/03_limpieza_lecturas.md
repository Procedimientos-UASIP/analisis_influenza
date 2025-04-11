# LIMPIEZA DE ARCHIVOS FASTQ
El ensamblaje de secuencias genómicas utilizando las secuencias obtenidas por secuenciación masiva se ve afectado por la calidad de las secuencias mismas. Por lo tanto, antes de realizar el alineamiento y ensamblaje de fragmentos genómicos, es prioridad eliminar todas aquellas lecturas, o partes de lecturas, cuyas calidades no sean óptimas.

Para cumplir lo anterior, es necesario realizar una evaluación inicial de todas las lecturas obtenidas por secuenciación masiva, en función de los valores de calidad asociado a cada posición de cada lectura. Posteriormente, se debe realizar una limpieza de las secuencias de acuerdo con la evaluación de los análisis obtenidos. Finalmente, es necesario una reevaluación de las calidades de las lecturas después de dicha limpieza para garantizar el mejoramiento de las calidades de secuencias en general.

Los códigos depositados a continuación son una resumen práctico del ***"PROCEDIMIENTO PARA EL ANÁLISIS DE CALIDAD, LIMPIEZA Y REANÁLISIS DE LAS SECUENCIAS PRODUCIDAS POR SECUENCIACION MASIVA MEDIANTE LA TÉCNICA DE ILLUMINA"***, <ins>(CPA-PD-XXXXX)</ins>, desarrollado e implementado por la UASIP-CENAPA-SENASICA.

## Código
Se recomieda colocar los archivos pareados crudos de la muestra a trabajar en un solo directorio (con el nombre de la muestra). Ej. */home/CPA-12345-24/R1.fastq.gz* y */home/CPA-12345-24/R2.fastq.gz* **(Los nombres de los archivos en las siguientes instrucciones deben ajustarse a los requerimientos de los archivos a usar)**. A partir de esta estructura de archivos, se procesan los archivos de la siguiente manera:

1. Ubicarse en el directorio de trabajo. Aquí deben estar los archivos de lecturas:
```bash
cd /home/CPA-12345-24
```

2. Crear un directorio para depositar el análisis de calidad de lecturas crudas
```bash
mkdir CALIDAD_CRUDA
```

3. Ejecutar el análisis de calidad.
```bash
fastqc R1.fastq.gz R1.fastq.gz --outdir ./CALIDAD_CRUDA --threads 10 
```

4. Revisar los archivos HTML de salida para identificar parámetros que ayuden a la limpieza. Esto se puede lograr transfiriendo los archivos HTML a una computadora local y visualizando en un navedador de internet.

5. Ejecutar el programa de limpieza. Ajustar parámetros de acuerdo a las necesidades. Los resultados se guardarán en ./TRIMMING/
```bash
trim_galore --paired --retain_unpaired --gzip --fastqc --cores 6 --clip_R1 18 --clip_R2 18 --three_prime_clip_R1 2 --three_prime_clip_R2 2 --length 180 -o TRIMMING/ R1.fastq.gz R1.fastq.gz
```

6. Revisar los archivos HTML de salida para verificar los resultados de la limpieza. Queda a consideración de la experiencia del usuario el quedar conforme con la limpieza realizada. 