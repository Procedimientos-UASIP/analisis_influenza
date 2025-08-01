#!/usr/bin/env Rscript

# ============================================
# CARGAR LIBRERIAS
# ============================================

# Cargar paquetes necesarios
suppressPackageStartupMessages({
  library(argparse)
  library(Biostrings)
  library(readr)
  library(dplyr)
})

# ============================================
# LEER ARGUMENTOS
# ============================================

# Definir parser de argumentos
parser <- ArgumentParser(description = "Script para analizar FASTA y tabla de cobertura")

parser$add_argument("--fasta", required = TRUE,
                    help = "Archivo FASTA con las secuencias (e.g. CPA-03878-25.fna)")
parser$add_argument("--prof", required = TRUE,
                    help = "Archivo de cobertura por posición (e.g. CPA-03878-25_prof.tsv)")

args <- parser$parse_args()

# ============================================
# LEER ARCHIVOS
# ============================================

# Leer archivo FASTA
secuencias <- readDNAStringSet(args$fasta)

# Leer tabla de cobertura
tabla_cobertura <- read_tsv(args$prof, col_types = cols()) |> select(-1)

# Convertir la tabla de cobertura a lista (cada columna es un vector)
lista_coberturas <- as.list(tabla_cobertura)

# ============================================
# ANALISIS
# ============================================

for (i in seq_along(secuencias)) {
  sec <- secuencias[[i]]
  letras <- strsplit(as.character(sec), "")[[1]]

  if (any(letras == "N")) {
    n_pos <- which(letras == "N")
    n_runs <- IRanges::IRanges(start = n_pos)

    cobertura <- lista_coberturas[[i]]
    offset <- 0

    for (j in seq_along(n_runs)) {
      inicio <- start(n_runs)[j] + offset
      largo  <- width(n_runs)[j]
      cobertura <- append(cobertura, rep(0, largo), after = inicio - 1)
      offset <- offset + largo
    }

    lista_coberturas[[i]] <- cobertura
  }
}

# Igualar longitudes
max_len <- max(sapply(lista_coberturas, length))
lista_coberturas <- lapply(lista_coberturas, function(x) {
  length(x) <- max_len
  x
})

# Volver a tabla
tabla_ajustada <- as_tibble(lista_coberturas)

# Agregar columna Posicion
tabla_ajustada <- tibble(Posicion = 1:max_len) |> bind_cols(tabla_ajustada)

write_tsv(tabla_ajustada, "./table_prof_corr.tsv", na = "")

