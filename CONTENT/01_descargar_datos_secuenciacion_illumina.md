# DESCARGAR DATOS DE SECUENCIACION MASIVA DE LA NUBE DE ILLUMINA

Los datos obtenidos del equipo de secuenciación masiva son inmediatamente depositados en una base de datos en linea de illumina, conocida como **basespace**. Este procedimiento indica la forma de explorar dicha nube y descargar los archivos pertinentes necesarios para hacer el ensamblaje de segmentos de virus de influenza aviar.

## Código

### Descarga de datos de secuenciacion masiva.

1. Se debe contar con el programa "bs" de ilummina descargado y ubicado. Para facilitar la operación, se incluye en la carpeta bin de este repositorio.

2. Authorizar el programa para que acceda a la nube
```bash
bs auth
```
Esto proporcionará una URL donde se deberán suministrar las credenciales apropiadas. Preguntar por ellas. Una vez configurada la autorización, el programa puede utilizarse.

3. Ubicar el ID del projecto a descargar con:
```bash
~/analisis_influenza/bin/bs project list
```

4. Armar la descarga, por ejemplo:
```bash
~/analisis_influenza/bin/bs project download --id=402166304 --output=PROYECTO --extension=fastq.gz
```

5. Tambien se puede descargar muestras individuales si se conoce su código de biomuestra:
```bash
~/analisis_influenza/bin/bs biosample list
```
6. Armar la descarga, por ejemplo:
```bash
~/analisis_influenza/bin/bs biosample download --id=745389546 --output=TEST --extension=fastq.gz
```

