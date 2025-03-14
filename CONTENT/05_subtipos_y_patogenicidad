# Análisis de subtipificación y patogenicidad

La Influenza aviar es una enfermedad viral causada por el virus de la gripe aviar (Virus de la Influenza Tipo A), que afecta principalmente a las aves, pero tiene el potencial de infectar a humanos y otros animales, por lo que es necesario la subtipificación y la evaluación de la patogenicidad de las cepas del virus, ya que son elementos críticos para la vigilancia, control y prevención de brotes de esta enfermedad en aves y en la salud pública. Evitando una transmisión rápida y mitigando con la prevención de vacunas y estrategias de bioseguridad. 

Los códigos depositados a continuación son una resumen práctico del "PROCEDIMIENTO PARA EL ANÁLISIS DE DATOS DE INFLUENZA AVIAR POR SECUENCIACIÓN MASIVA MEDIANTE LA TÉCNICA DE ILLUMINA (SUBTIPIFICACIÓN Y PATOGENICIDAD)", con código CPA-PD-#####, desarrollado e implementado por la UASIP-CENAPA-SENASICA.

## Código
Se deberá estar posicionado en la carpeta de trabajo de la muestra, ej. */home/CPA-12345-24/*.  Los nombres de los archivos en las siguientes instrucciones deben ajustarse a los requerimientos de los archivos a usar. Se indica con corchetes los argumentos que deben/pueden modificarse. A partir de esta estructura de archivos, se procesan los archivos de la siguiente manera:

### Alineamiento

1. Hacer carpeta de resultado finales
``bash
mkdir resultados_enviar 
```
En esa carpeta se colocarán todos los archivos fasta y graficas que se hicieron para las muestras o muestra

2. Crear carpeta de estadisticas en TRIMMED
``bash
mkdir stats_all  
```

3. Copiar ensambles de cada segmento  a la carpeta stats_all
```bash
cat ../s{1..8}/s{1..8}_$(pwd | awk -F'/' '{print $(NF-2)}').fasta 2>/dev/null | \
  awk '/^>/ { if (NR != 1) print ""; print $0; next } {printf "%s", $0 }' > $(pwd | awk -F'/' '{print $(NF-2)}').fasta 




