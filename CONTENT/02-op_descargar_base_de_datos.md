# OPCIONAL: DESCARGAR BASES DE DATOS PRECOMPILADAS (NCBI Y KRAKEN2)

Ademas de descargar las secuencias correspondientes a cada segmento del genoma del virus de influenza aviar, existe la posibilidad de descargar bases de datos precompiladas y construidas con secuencias de numerosos virus, bacterias, archeas, y otros organismos.

Este procedimiento pretende establecer los pasos para descargar la base de datos de NCBI denominada "viral", y la base de datos de KRAKEN2, denominada "Standar"

## Código

### Descarga de base de datos NCBI viral (nt_viruses)

1. Durante la instalación de blast, se instala un pequeño programa que nos permite descargar de forma local bases de datos depositadas en NCBI. Para descargar toda la base de datos correspondiente a virus, realizar:
```bash
update_blastdb.pl --source "aws" --num_threads 20 nt_viruses
```

### Descarga de base de datos NCBI principal condensada (core_nt)
2. Para descargar toda la base de datos correspondiente a virus, realizar:
```bash
update_blastdb.pl --source "aws" --num_threads 20 core_nt
```

3. Es probable que esta base de datos pida los datos de la base de datos taxonomicos. Descargarlos en una carpeta independiente y copiar archivos importantes con:
```bash
update_blastdb.pl --source "aws" --num_threads 20 core_nt
update_blastdb.pl --source "aws" --num_threads 20 taxdb
```
Los archivos de la base de datos de taxonomia a copiar son taxdb.btd, taxdb.bti y taxonomy4blast.sqlite3


### Descarga de base de datos KRAKEN2 Standar


