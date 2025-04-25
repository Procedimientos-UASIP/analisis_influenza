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

### Ensamble
3. Hacer el ensamblaje en dicha ubicación. Este paso implica especificar todos los kmeros a usar y, cuando se obtengan ensambles a nivel de scaffolds, realizar un blast de los resultados obtenidos para revisar el kmero que realizó el mejor ensamblade. Como esto implica mucho procesamiento en bucle, se encapsuló el proceso en un sólo script. **Cuidar de modificar archivos de entrada, kmeros, output y log**. Un ejemplo de su implementación es la siguiente:

```bash
echo "S" | nohup ~/analisis_influenza/bin/ensamble.sh -a ../../../ALINEAMIENTO/BWA/S1/s1_reads_r1.fq.gz -b ../../../ALINEAMIENTO/BWA/S1/s1_reads_r2.fq.gz -x ../../../ALINEAMIENTO/BWA/S1/s1_reads_u1.fq.gz -y ../../../ALINEAMIENTO/BWA/S1/s1_reads_u2.fq.gz -t 26 --kini 11 --kfin 127 -o OUT_11_127 >log_11_127 2>&1 &
```

Realizar el análisis hasta los kmeros que se consideren pertinentes (min. 11, max. 127).

4. Ejecutar el script para unir todos los BLAST resultantes depositatos en las carpetas llamadas OUT_*. Ejecutar aunque sólo haya un directorio de salida.
```bash
~/analisis_influenza/bin/merge_blast_results.sh
```

5. Revisar el archivo resultante para asegurar que existan resultados de blast con algunos kmeros.

6. Procesar el archivo con el siguiente script para identificar la el mejor alineamiento, más extenso y con la mejor covertura. Se generará además la carpeta de PROFUNDIDAD con la profundidad obtenida para el segmento. Modificar el nombre de archivo correspondiente:
```bash
~/analisis_influenza/bin/ensamble_analyze_blast.sh S1_blastn-careful_merged.txt
```

7. Repetir desde el paso 2, modificando cada carpeta correspondiente a cada segmento.

8. Para unir todos las profundidades de todos los segmentos, posicionarse en la carpeta con todas las subcarpetas de segmentos (ENSAMBLE/SPADES/) y ejecutar:
```bash
~/analisis_influenza/bin/unir_profundidades.sh
```