# Análisis de secuencias de virus influenza
Este repositorio contiene los scripts necesarios para el análisis de la secuenciación masiva de muestras de virus de influenza, que incluye la limpieza, alineamiento y ensamblaje de las lecturas obtenidas. Su objetivo principal es encapsular, en la medida de lo posible, los procedimientos necesarios para el procesamiento bioinformático de muestras de secuenciación masiva, documentar los mismos, y garantizar la reproducibilidad de los análisis.  

Los scripts informáticos aquí depositados son utilizados por la Unidad de Atención Sanitaria Integral en Polinizadores (UASIP) del Centro Nacional de Referencia en Parasitología Animal y Tecnología Analítica (CENAPA) del Servicio Nacional de Sanidad Inocuidad y Calidad Agroalimentaria (SENASICA). Dichas organizaciones son propietarias de la información contenida en este repositorio.

# Secciones

## CONTENT
1. [Descargar lecturas de la nube de Illumina](https://github.com/Procedimientos-UASIP/analisis_influenza/blob/main/CONTENT/01_descargar_datos_secuenciacion_illumina.md)
2. [Crear base de datos de virus (NCBI Virus) para BLAST y BWA-MEM2](https://github.com/Procedimientos-UASIP/analisis_influenza/blob/main/CONTENT/02_creaci%C3%B3n_base_de_datos.md)
3. [OPCIONAL: Descargar base de datos precompiladas (NCBI y KRAKEN2)](https://github.com/Procedimientos-UASIP/analisis_influenza/blob/main/CONTENT/02-op_descargar_base_de_datos.md)
4. [Limpiar lecturas](https://github.com/Procedimientos-UASIP/analisis_influenza/blob/main/CONTENT/03_limpieza_lecturas.md)
5. [Alinear lecturas](https://github.com/Procedimientos-UASIP/analisis_influenza/blob/main/CONTENT/04_alineamientos.md)
6. [Ensamblar lecturas](https://github.com/Procedimientos-UASIP/analisis_influenza/blob/main/CONTENT/05_ensambles.md)
7. [Análisis de patogenicidad](https://github.com/Procedimientos-UASIP/analisis_influenza/blob/main/CONTENT/06_subtipos_y_patogenicidad.md)
8. Análisis filogenético

## NOTEBOOKS

Se incluyen archivos de jupyter notebook con kernel de bash para tener una plantilla de trabajo que condense cada una de las operaciones necesarias contenidas en CONTENT.