# Análisis de subtipificación y patogenicidad

La Influenza aviar es una enfermedad viral causada por el virus de la gripe aviar (Virus de la Influenza Tipo A), que afecta principalmente a las aves, pero tiene el potencial de infectar a humanos y otros animales, por lo que es necesario la subtipificación y la evaluación de la patogenicidad de las cepas del virus, ya que son elementos críticos para la vigilancia, control y prevención de brotes de esta enfermedad en aves y en la salud pública. Evitando una transmisión rápida y mitigando con la prevención de vacunas y estrategias de bioseguridad. 

Los códigos depositados a continuación son una resumen práctico del "PROCEDIMIENTO PARA EL ANÁLISIS DE DATOS DE INFLUENZA AVIAR POR SECUENCIACIÓN MASIVA MEDIANTE LA TÉCNICA DE ILLUMINA (SUBTIPIFICACIÓN Y PATOGENICIDAD)", con código CPA-PD-#####, desarrollado e implementado por la UASIP-CENAPA-SENASICA.

## Código

Se deberá estar posicionado en la carpeta correspondiente al segmento 4. Se debe suministrar al script el archivo que tiene el mejor ensamble con su direccionalidad corregida.

1. Hacer una traduccion de los 3 marcos de lectura para identificar el marco paros en la seccion media (correspondiente al ORF).
``bash
seqkit translate [BEST_RESULT_S4_SEQUENCE_CORRECTED] -f1,2,3 -F 
```

2. Ejecutar el script para determinar la presencia del motivo caracteristico de alta patogenicidad (PXnGLF, donde Xn son residuos positivos (K, H, R), principalmente)
``bash
~/analisis_influenza/bin/analisis_patogenicidad.sh  [BEST_RESULT_S4_SEQUENCE_CORRECTED]
```

3. El analisis se hace en todos los marcos de lectura de la secuencia original y su reverso complementario, para cubrir todas las opciones. Identificar el marco donde apareció el ORF en la dirección sentido (+). Debería mostrarse el resultado del analisis (HPAI o LPAI) y la secuencia encontrada.
