# OPCIONAL: DESCARGAR BASES DE DATOS PRECOMPILADAS (NCBI Y KRAKEN2)

Ademas de descargar las secuencias correspondientes a cada segmento del genoma del virus de influenza aviar, existe la posibilidad de descargar bases de datos precompiladas y construidas con secuencias de numerosos virus, bacterias, archeas, y otros organismos.

Este procedimiento pretende establecer los pasos para descargar la base de datos de NCBI denominada "viral", y la base de datos de KRAKEN2, denominada "Standar"

## Esquema del filesystem final para bases de datos preformateadas

Se recomienda que la estructura final de los archivos tenga la siguiente estructura, en una partición independiente (/backup) dedicada a depositar datos que no serán modificados con frecuencia.

```
/backup/
        └── DATABASES
            ├── KRAKEN2
            │   ├── STANDARD
            │   └── VIRAL
            │       ├── library
            │       └── taxonomy
            └── NCBI
                ├── CORE_NT
                ├── NT_VIRUSES
                └── TAXDB
```

## Código

### Descarga de base de datos NCBI viral (nt_viruses)

1. Durante la instalación de blast, se instala un programa que nos permite descargar de forma local bases de datos depositadas en NCBI. Para descargar toda la base de datos correspondiente a virus (aprox. 64 Gb), ejecutar el siguiente comando. Se recomienda ejecutarlo en /backup/DATABASES/NCBI/NT_VIRUSES/.
```bash
update_blastdb.pl --source "aws" --num_threads 20 nt_viruses
```

### Descarga de base de datos NCBI principal condensada (core_nt)
2. Para descargar toda la base de datos correspondiente al núcleo central de secuencias de NCBI (aprox. 237 Gb), ejecutar el siguiente comando. Se recomienda ejecutarlo en /backup/DATABASES/NCBI/CORE_NT/.
```bash
update_blastdb.pl --source "aws" --num_threads 20 core_nt
```

3. Es probable que esta base de datos pida los datos de la base de datos taxonómicos (aprox. 261 Mb). Descargarlos en una carpeta independiente. Se recomienda ejecutarlo en /backup/DATABASES/NCBI/TAXDB/.
```bash
update_blastdb.pl --source "aws" --num_threads 20 taxdb
```
4. Los archivos de la base de datos de taxonomia a copiar son taxdb.btd, taxdb.bti y taxonomy4blast.sqlite3. Se deben copiar en /backup/DATABASES/NCBI/CORE_NT/

### Descarga de base de datos KRAKEN2 Viral

5. Al instalar kraken2, se instalan programas para facilitar bases de datos de ciertos organismos. Sin embargo, existen otras bases de datos que cuentan con información de varios organismos, por lo que no es necesario descargarlas individualmente y luego unirlas. Para eso, se copiará la dirección del enlace a descargar. Visitar https://benlangmead.github.io/aws-indexes/k2 e identificar el link **.tar.gz** de la colección **"Viral" (0.5 Gb)**. Luego, usar un descargador de alta velocidad para realizar la descarga. Por ejemplo:

```bash
aria2c -c -x 16 https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20250402.tar.gz
```
6. Descomprimir el archivo resultante usando varios núcleos para acelerar la descompresión:

```bash
pigz -p 20 -dc k2_viral_20250402.tar.gz | tar -xf -
```

### Descarga de base de datos KRAKEN2 Standard

7. Los pasos son similares a los anteriores salvo por la liga de la base de datos, que en este caso será la colección **"Standard" (70 Gb)**:

```bash
aria2c -c -x 16 https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20250402.tar.gz

pigz -p 20 -dc k2_standard_20250402.tar.gz | tar -xf -
```
