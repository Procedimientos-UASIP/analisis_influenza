# CREAR BASE DE DATOS DE VIRUS (NCBI Virus)

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
14. En el servidor, revisar que los archivos se hayan transferido:
```bash
ls -la
```

### Limpieza de la información en encabezados de secuencias de cada segmento

15.	Corregir los encabezados de cada fasta. Para eso, se diseñó un script con los comandos necesarias. Sólo funcionará si se han seguido las instrucciones de este procedimiento al pie de la letra hasta este punto. Se incluye en un loop *for* para operar sobre todos los archivos transferidos:
```bash
#Revisar primero con esta instruccion. Sólo deben haber una cantidad por archivo fasta de segmento procesado.
#for i in {1..8}; do ~/analisis_influenza/bin/format_headers_db_virus.sh seg_"$i"_nucl.fasta.gz | grep ">" | cut -f1 -d"|" | sort | uniq -c; done

# Si lo genera discrepancia, reportar para revisar el script de limpieza. Posteriormente:
for i in {1..8}; do ~/analisis_influenza/bin/format_headers_db_virus.sh seg_"$i"_nucl.fasta.gz | pigz -9 -p20 -c >clean_seg_"$i"_nucl.fna.gz; done
```

16.	Crear el sistema de archivos donde se construirán las bases de datos indexadas que serán utilizadas por BLAST y BWA-MEM2:
```bash
mkdir –p /backup/DATABASES/UASIP/BLAST/INFLUENZA
mkdir –p /backup/DATABASES/UASIP/BWAMEM2/INFLUENZA
```

17.	Ingresar al directorio donde se construirá  la base de datos de influenza para BLAST:
```bash
cd /backup/DATABASES/UASIP/BLAST/INFLUENZA
```

18.	Crear un archivo temporal con todas las secuencias comprimidas que conformarán la base de datos:
```bash
cat /backup/RAW/UASIP_DB_ALPHAINFLUENZAVIRUS_2025_04/clean_seg_* > influenza_db.fna.gz
```

19.	Construir la base de datos. Se realiza una descompresión previa en tubería para no generar un archivo intermediario:
```bash
zcat influenza_db.fna.gz | makeblastdb -in - -dbtype nucl -out influenza_db -title influenza_db
```

20.	Mover el archivo temporal a la carpeta donde se construirá la base de datos de BWAMEM2:
```bash
mv influenza_db.fna.gz /backup/DATABASES/UASIP/BWAMEM2/INFLUENZA
```

21.	Ingresar al directorio donde se construirá la base de datos de influenza para BWAMEM2:
```bash
cd /backup/DATABASES/UASIP/BWAMEM2/INFLUENZA
```

22.	BWAMEM2 no puede construir la base de datos en una tubería, por lo que hay que descomprimir el archivo fasta:
```bash
pigz -k -d -p 20 influenza_db.fna.gz
```

21.	Armar la base de datos para BWAMEM2. Usar **ambiente conda alineamiento**. Puede tardar unos 15 minutos:
```bash
conda activate alineamiento

bwa-mem2 index influenza_db.fna
```
