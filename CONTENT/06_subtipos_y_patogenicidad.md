# Análisis de subtipificación y patogenicidad

La Influenza aviar es una enfermedad viral causada por el virus de la gripe aviar (Virus de la Influenza Tipo A), que afecta principalmente a las aves, pero tiene el potencial de infectar a humanos y otros animales, por lo que es necesario la subtipificación y la evaluación de la patogenicidad de las cepas del virus, ya que son elementos críticos para la vigilancia, control y prevención de brotes de esta enfermedad en aves y en la salud pública. Evitando una transmisión rápida y mitigando con la prevención de vacunas y estrategias de bioseguridad. 

Los códigos depositados a continuación son una resumen práctico del "PROCEDIMIENTO PARA EL ANÁLISIS DE DATOS DE INFLUENZA AVIAR POR SECUENCIACIÓN MASIVA MEDIANTE LA TÉCNICA DE ILLUMINA (SUBTIPIFICACIÓN Y PATOGENICIDAD)", con código CPA-PD-#####, desarrollado e implementado por la UASIP-CENAPA-SENASICA.

## Código

### Subtipificación

Si bien la identificación de un genoma viral en su totalidad se da en función de la identidad de sus segmentos 4 y 6 (que codifican para las proteinas hemaglutinina y neuraminidasa, respectivamente), es posible darles una identidad a los otros segmentos en función de su presencia en distintos contextos HxNy, lo cual a su vez es una evidencia indirecta de eventos de recombinación genética entre segmentos de distintos virus. Para determinar el subtipo viral, se hará un blast de la secuencia ensamblada de cada segmento contra la base de datos de virus de NCBI, se obtendrán todos los subjects con los que se hagan match, y se obtendrá un porcentaje que representará la cantidad de subjects que hicieron match con determinado segmento HxNy.

Para realizar lo anterios, se diseño un script con todas las instrucciones para realizar lo anterior:

1. Ubicarse en la carpeta correspondiente a cada segmento ensamblado, donde se encuentra el archivo fasta con el mejor ensamble corregido.

2. Ejecutar el script indicando el archivo que va a utilizarse como query:
```bash
~/analisis_influenza/bin/analisis_subtipificacion.sh best_result_S1_sequence_corrected.fna
```

3. Se genera un archivo Sx_subtipos_top5.tsv, que contiene una tabla con los principales subtipos contra los que se hizo match. 

### Patogenicidad

Se deberá estar posicionado en la carpeta correspondiente al segmento 4. Además, se debe suministrar al script el archivo que tiene el mejor ensamble con su direccionalidad corregida.

4. Hacer una traduccion de los 3 marcos de lectura para identificar el marco paros en la seccion media (correspondiente al ORF).
```bash
seqkit translate [BEST_RESULT_S4_SEQUENCE_CORRECTED] -f1,2,3 -F 
```

5. Ejecutar el script para determinar la presencia del motivo caracteristico de alta patogenicidad (PXnGLF, donde Xn son residuos positivos (K, H, R), principalmente)
```bash
~/analisis_influenza/bin/analisis_patogenicidad.sh  [BEST_RESULT_S4_SEQUENCE_CORRECTED]
```

6. El analisis se hace en todos los marcos de lectura de la secuencia original y su reverso complementario, para cubrir todas las opciones. Identificar el marco donde apareció el ORF en la dirección sentido (+). Debería mostrarse el resultado del analisis (HPAI o LPAI) y la secuencia encontrada.