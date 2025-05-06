# ENSAMBLE *DE NOVO* DE GENOMAS DE INFLUENZA 

El ensamble genómico es un procedimiento en biología computacional que permite la reconstrucción de secuencias genómicas largas, a partir de fragmentos más pequeños (llamados "lecturas" o "reads") obtenidos por secuenciación masiva, por ejemplo, con Ilumina. 

Los códigos depositados a continuación son una resumen práctico del "PROCEDIMIENTO PARA EL ANALISIS DE DATOS DE INFLUENZA AVIAR POR SECUENCIACIÓN MASIVA MEDIANTE LA TÉCNICA DE ILLUMINA  
 (ENSAMBLE)", <ins>(CPA-PD-XXXXX)</ins>, desarrollado e implementado por la UASIP-CENAPA-SENASICA.

 ## Código
Se deberá estar posicionado en la carpeta de trabajo de la muestra, ej. */home/CPA-12345-24/* **(Los nombres de los archivos en las siguientes instrucciones deben ajustarse a los requerimientos de los archivos a usar)**. A partir de esta estructura de archivos, se procesan los archivos de la siguiente manera:

### Preparacion
1. Hacer directorio para ensambles y ensamblador a usar. Incluir subdirectorios para ensambles de cada segmento:
```bash
mkdir -p ENSAMBLE/SPADES/S{1..8}
```

2. Entrar al directorio del segmento a ensamblar
```bash
cd ENSAMBLE/SPADES/S1
```

3. Activar ambiente conda de los ensambladores
```bash
conda activate ensamble
```

### Ensamble
4. Hacer el ensamblaje en dicha ubicación. Este paso implica especificar todos los kmeros a usar y, cuando se obtengan ensambles a nivel de scaffolds, realizar un blast de los resultados obtenidos para revisar el kmero que realizó el mejor ensamblade. Como esto implica mucho procesamiento en bucle, se encapsuló el proceso en un sólo script. **Cuidar de modificar archivos de entrada, kmeros, output y log**. Un ejemplo de su implementación es la siguiente:

```bash
echo "S" | nohup ~/analisis_influenza/bin/ensamble.sh -a ../../../ALINEAMIENTO/BWA/S1/s1_reads_r1.fq.gz -b ../../../ALINEAMIENTO/BWA/S1/s1_reads_r2.fq.gz -x ../../../ALINEAMIENTO/BWA/S1/s1_reads_u1.fq.gz -y ../../../ALINEAMIENTO/BWA/S1/s1_reads_u2.fq.gz -t 26 --kini 11 --kfin 127 -o OUT_11_127 >log_11_127 2>&1 &
```

Realizar el análisis hasta los kmeros que se consideren pertinentes (min. 11, max. 127).

### Mejor ensamble y cálculo de profundidad
5. Ejecutar el script para unir todos los BLAST resultantes depositatos en las carpetas llamadas OUT_*. Ejecutar aunque sólo haya un directorio de salida.
```bash
~/analisis_influenza/bin/merge_blast_results.sh
```

6. Revisar el archivo resultante para asegurar que existan resultados de blast con algunos kmeros.

7. Procesar el archivo con el siguiente script para identificar la el mejor alineamiento, más extenso y con la mejor covertura. Se generará además la carpeta de PROFUNDIDAD con la profundidad obtenida para el segmento. Ejecutar en ambiente alineamiento. Modificar el nombre de archivo correspondiente:
```bash
conda activate alineamiento
~/analisis_influenza/bin/ensamble_analyze_blast.sh S1_blastn-careful_merged.txt
```
NOTA: En caso de observar el mensaje "\[W::bam_merge_core2\] No @HD tag found.", ignorarlo, pues es sólo una advertencia de samtools.

8. Repetir desde el paso 2, modificando cada carpeta correspondiente a cada segmento.

### Unir profundidades y calcular porcentaje de lecturas usadas
9. Para unir todos las profundidades de todos los segmentos, posicionarse en la carpeta con todas las subcarpetas de segmentos (ENSAMBLE/SPADES/) y ejecutar el script indicando la ubicación de cada archivo a unir:
```bash
~/analisis_influenza/bin/unir_profundidades.sh --S1 S1/COBERTURA/S1_cobertura --S2 S2/COBERTURA/S2_cobertura --S3 S3/COBERTURA/S3_cobertura --S4 S4/COBERTURA/S4_cobertura --S5 S5/COBERTURA/S5_cobertura --S6 S6/COBERTURA/S6_cobertura --S7 S7/COBERTURA/S7_cobertura --S8 S8/COBERTURA/S8_cobertura
```

10. Ejecutar este script para calcular el porcentaje de lecturas ensambladas de acuerdo a las que se filtraron por el alineamiento
```bash
~/analisis_influenza/bin/ensamble_cal_porc_lect_mapeadas.sh
```

### Graficar profundidades y uso de lecturas
11. Activar ambiente conda para graficar en R
```bash
conda activate R
```

12. Graficar profundidades. Indicar archivo donde se unieron las profundidades, así como el nombre de la muestra:
```bash
~/analisis_influenza/bin/graficar_profundidad.R --input_file coberturas_finales.tsv --muestra CPA-00245-25-P1-1
```

13.  Graficar porcentaje de lecturs usadas. Indicar archivo donde se unieron los datos de las lecturas usadas, así como el nombre de la muestra:
```bash
~/analisis_influenza/bin/graficar_lecturas.R --input_file uso_de_lecturas.tsv --muestra CPA-00245-25-P1-1
```