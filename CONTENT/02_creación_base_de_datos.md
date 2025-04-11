# CREACION DE BASE DE DATOS DE SECUENCIAS DE INFLUENZA

Una base de datos es una colección organizada de datos, los cuales se almacenan y gestionan de manera estructurada, permitiendo el acceso y/o modificación de los datos contenidos. El principal propósito de una base de datos es facilitar el almacenamiento de grandes volúmenes de datos y permitir que los usuarios puedan trabajar con ellos de manera eficiente.

Una base de datos sobre secuencias de virus de influenza aviar sirve para gestionar, analizar y extraer información sobre las secuencias de dicho virus. Este tipo de base de datos es esencial para la investigación genética del virus, la detección de variantes y el seguimiento de mutaciones clave que pueden influir en la transmisión y patogenicidad del virus.

La construcción de una base de datos para realizar análisis con programas bioinformáticos implica que las secuencias que conformarán la base de datos sean formateadas e indexadas para eficientizar la consulta de los datos. Diferentes programas formatean e indexan un mismo conjunto de secuencias de forma diferente, por lo que es importante preparar cada base de datos para cada programa a usar.

Los códigos depositados a continuación son una resumen práctico del ***"PROCEDIMIENTO PARA LA CREACIÓN DE BASES DE DATOS CON SECUENCIAS GENÓMICAS DE VIRUS DE INFLUENZA AVIAR PARA EJECUTAR BLAST Y BWA-MEM2"*** <ins>(CPA-PD-XXXXX)</ins>, desarrollado e implementado por la UASIP-CENAPA-SENASICA.

Las direcciones que apuntan a los archivos están diseñados para trabajar en el servidor cenapa01. **Modificar conforme sea necesario**

## Código

### Obtención de datos que conformarán la base de datos

1.	Acceder a la base de datos de NCBI Virus (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/):
2.	Seleccionar el icono "Search by virus".
3.	Escribir “Alphainfluenzavirus influenzae”:
4.	Seleccionar el taxón "2955291". 
5.	Buscar "Segment", en la columna "Refine Results", y seleccionar el segmento a descargar (1 al 8, uno a la vez). 
6.	Seleccionar "Download".
    1. Seleccionar "Nucleotide".
    2.	Seleccionar "Download All Records".
    3.	Seleccionar "Build custom". Además de los campos "Accesion” y “GenBank Title”, agregar “Host”, “Country”, “Segment” y “Genotype”.
    4.	Dar clic en “Download”. Cuidar que las descargas concluyan con éxito.
    5.	Renombrar el archivo a seg_1_nucl.fna
7.	Repetir desde el paso 5 para todos los segmentos virales (1-8). 
8.	Comprimir los archivos con pigz. Ejemplo de compresión en Windows con Powershell:
```powershell
& "C:\Users\david.rendon.i\Downloads\pigz-win32\pigz.exe" -9 -p 4 seg_1_nucl.fna
```

10.	En el servidor, crear el sistema de archivos donde se depositarán los archivos descargados. Modificar las fechas al año y mes en que se actualiza la base de datos:
```bash
mkdir -p /backup/RAW/UASIP_DB_ALPHAINFLUENZAVIRUS_2025_04
```

11.	Ubicarse en dicho directorio:
```bash
cd /backup/RAW/UASIP_DB_ALPHAINFLUENZAVIRUS_2025_04
```

13.	En la computadora personal, usar scp para cargar los archivos comprimidos al servidor. Un ejemplo de transferencia es:
```bash
scp seg_1_nucl.fna.gz admcenapa@10.24.34.18:/backup/RAW/UASIP_DB_ALPHAINFLUENZAVIRUS_2025_04
```

### Limpieza de información en cada segmente

12.	En el servidor, corregir los encabezados de cada fasta. Para eso, se diseñó un script con los comandos necesarias. Sólo funcionará si se han seguido las instrucciones de este procedimiento al pie de la letra. Hacer un loop for para operar sobre todos los archivos:
```bash
for i in {1..8}; do ~/bin/format_headers_db_virus.sh seg_"$i"_nucl.fasta.gz | pigz -9 -p20 -c >clean_seg_"$i"_nucl.fna.gz; done
```
14.	Crear el sistema de archivos donde se construirán las bases de datos indexadas que serán utilizadas por BLAST y BWA-MEM2:
```bash
mkdir –p /backup/DATABASES/UASIP/BLAST/INFLUENZA
mkdir –p /backup/DATABASES/UASIP/BWAMEM2/INFLUENZA
```

15.	Ingresar al directorio donde se construirá  la base de datos de influenza para BLAST:
```bash
cd /backup/DATABASES/UASIP/BLAST/INFLUENZA
```
16.	Crear un archivo temporal con todas las secuencias que conformarán la base de datos:
```bash
cat /backup/RAW/UASIP_DB_ALPHAINFLUENZAVIRUS_2025_04/clean_seg_* > influenza_db.fna.gz
```

17.	Construir la base de datos. Se realiza una descompresión previa en tubería para no generar un archivo intermediario:
```bash
zcat influenza_db.fna.gz | makeblastdb -in - -dbtype nucl -out influenza_db -title influenza_db
```

18.	Mover el archive temporal a la carpeta donde se construirá la base de datos de BWAMEM2:
```bash
mv influenza_db.fna.gz /backup/DATABASES/UASIP/BWAMEM2/INFLUENZA
```

19.	Ingresar al directorio donde se construirá la base de datos de influenza para BWAMEM2:
```bash
cd /backup/DATABASES/UASIP/BWAMEM2/INFLUENZA
```

20.	BWAMEM2 no puede construir la base de datos en una tubería, por lo que hay que descomprimir el archivo fasta:
```bash
pigz -k -d -p 20 influenza_db.fna.gz
```

21.	Armar la base de datos para BWAMEM2. Usar ambiente conda alineamiento. Puede tardar unos 15 minutos:
```bash
conda activate alineamiento

bwa-mem2 index influenza_db.fna
```

 


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

5. Unir archivos fasta. Procurar seguir la siguiente convención para el nombre del archivo resultante:

&emsp; [nombre del virus]\_[tipo de secuencias]_[fecha de creación(Año-Mes)].fasta

 ```bash
cat s*.fasta >influenza_nucl_2024-oct.fasta 
```

6. Crear link simbólico en los lugares donde se indexará la base de datos.
```bash
cd ~/DATABASES/BLAST/INFLUENZA 
ln -s ~/DATABASES/RAW/NCBI_VIRUS/INFLUENZA/influenza_nucl_2024-11.fna 

cd ~/DATABASES/BWAMEM2/INFLUENZA 
ln -s ~/DATABASES/RAW/NCBI_VIRUS/INFLUENZA/influenza_nucl_2024-11.fna 
```

7. Indexar bases de datos. Revisar estár en el ambiente conda donde se encuentren los programas para indexar las bases de datos.
```bash
cd ~/DATABASES/BLAST/INFLUENZA 
makeblastdb –in influenza_nucl.fasta –dbtype nucl

cd ~/DATABASES/BWAMEM2/INFLUENZA 
bwa-mem2 index influenza_nucl.fasta 
```

