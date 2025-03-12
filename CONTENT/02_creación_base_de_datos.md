# LIMPIEZA DE ARCHIVOS FASTQ

Una base de datos es una colección organizada de datos, los cuales se almacenan y gestionan de manera estructurada, permitiendo el acceso a ellos, su modificación, y su recuperación. El principal propósito de una base de datos es facilitar el almacenamiento de grandes volúmenes de datos y permitir que los usuarios puedan trabajar con ellos de manera eficiente. 

Una base de datos sobre secuencias de virus de influenza aviar es una implementación específica que permite a un grupo de usuarios gestionar, analizar y extraer información sobre las secuencias de dicho virus. Este tipo de base de datos es esencial para la investigación genética del virus, la detección de variantes y el seguimiento de mutaciones clave que pueden influir en la transmisión y patogenicidad del virus. 

Particularmente, la construcción de una base de datos para análisis de BLAST implica la organización de miles de secuencias que serán utilizadas por los algoritmos de BLAST para encontrar secuencias de referencias, usando alineamientos locales, que permitan identificar la naturaleza de las secuencias que sean analizadas. 

Los códigos depositados a continuación son una resumen práctico del ***"PROCEDIMIENTO PARA LA CREACIÓN DE UNA BASE DE DATOS CON SECUENCIAS GENÓMICAS DE VIRUS DE INFLUENZA AVIAR PARA EJECUTAR BLAST Y BWA-MEM2"***, desarrollado e implementado por la UASIP-CENAPA-SENASICA.

## Código
 
### Obtención de datos que conformarán la base de datos

1. Acceder a la base de datos de NCBI Virus (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/). 
2. Seleccionar el icono "Search by virus" y escribir la palabra “Influenza A”. 
3. Seleccionar el taxón identificado como "taxid:11320", ya que es el identificador taxonómico del virus de influenza aviar y sus subtipos. 
4. Al desplegarse los resultados, identificar la sección Refine Results en la columna izquierda y buscar la opción identificada como "Segment". Seleccionar el segmento el cual se desea descargar (en nuestro caso del 1 al 8, uno a la vez).  
5. Una vez seleccionado el segmento, se descarga seleccionando la opción "Download".  
6. Seleccionar "Nucleotide" y siguiente.
7. Seleccionar "Download All Records" y siguiente.
8. Seleccionar "Build custom" y agregar los campos "Accesion”, “GenBank Title”, “Host, Country”, y “Genotype”. Dichos campos aparecerán en el encabezado de cada secuencia.
9. Dar clic en en Download y guardar como segment_1_nucl. 
10. El proceso de descarga se realizará para cada uno de los ocho segmentos, modificando el nombre del archivo final de acuerdo con el segmento que se haya descargado, teniendo en total 8 archivos fasta correspondientes a los ocho segmentos del virus de influenza. 

### LIMPIEZA DE INFORMACIÓN PARA CADA SEGMENTO

1. Crear el sistema de archivos en el servidor para la base de datos
```bash
mkdir –p ~/DATABASES/RAW/NCBI_VIRUS/INFLUENZA 
```

2. Crear el sistema de archivos donde se armará y depositará las bases de datos indexadas que serán utilizadas por BLAST y BWA-MEM2 
```bash
mkdir –p ~/DATABASES/BLAST/INFLUENZA 
mkdir –p ~/DATABASES/BWAMEM2/INFLUENZA 
```

3. Subir los arvhivos al servidor
```bash
scp  segment_1_nucl.fasta [USUARIO@DIRECCION_IP]:/home/DATABASES/RAW/NCBI_VIRUS/INFLUENZA 
```

4. Corregir cada archivo fasta
```bash
cat segment_1_nucl.fasta | \
  sed 's/ /_/g' | \
  sed 's/))[^|]*|/))/g' | \
  sed 's/Influenza_A_virus_//g'| \
  sed 's/^>\([^|]*\)/>segmento_1|\1/g' \
  >s1.fasta && rm segment_1_nucl.fasta 
```
