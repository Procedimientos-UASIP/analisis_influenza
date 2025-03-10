# Limpieza de archivos fastqc

# Código
```bash
mkdir CALIDAD_CRUDA 
```

```bash
fastqc [INPUT_1] [INPUT_1] --outdir ./CALIDAD_CRUDA --threads 10 
```

```bash
trim_galore --paired --retain_unpaired --gzip --fastqc --cores 6 --clip_R1 18 --clip_R2 18 --three_prime_clip_R1 2 --three_prime_clip_R2 2 --length 180 --o TRIMMING/ [INPUT_1] [INPUT_2]
```
