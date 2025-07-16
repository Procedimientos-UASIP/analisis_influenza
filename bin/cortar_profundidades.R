#!/usr/bin/env Rscript

# ==========================================
# cortar_profundidades.R
# ==========================================
# Uso:
# Rscript cortar_profundidades.R \
#   --profundidades archivo.tsv \
#   --rangos rangos.tsv \
#   --n_inicio 15 \
#   --n_final 15
# ==========================================

# ---------------------------
# 1. Leer argumentos
# ---------------------------

args <- commandArgs(trailingOnly = TRUE)

# ---------------------------
# 2. Función auxiliar para extraer valor
# ---------------------------

get_arg <- function(name, default = NULL) {
  if (name %in% args) {
    pos <- match(name, args)
    if (length(args) >= pos + 1) {
      return(args[pos + 1])
    } else {
      stop(paste("Falta valor para argumento:", name))
    }
  } else if (!is.null(default)) {
    return(default)
  } else {
    stop(paste("Argumento requerido:", name))
  }
}

# ---------------------------
# 3. Extraer cada argumento nombrado
# ---------------------------
archivo_profundidades <- get_arg("--profundidades")
archivo_rangos <- get_arg("--rangos")
start_offset <- as.numeric(get_arg("--n_inicio"))
end_offset <- as.numeric(get_arg("--n_final"))

# ---------------------------
# 4- Ejecutar Script
# ---------------------------

profundidades <- read.delim(archivo_profundidades)
coord_cds <- read.delim(archivo_rangos)

segmentos_cortados <- list()

for (i in 1:8) {
  start <- coord_cds$Start[i] - 15
  end   <- coord_cds$End[i] + 15
  
  # Filtrar las filas con POS dentro del rango
  idx <- which(profundidades$POS >= start & profundidades$POS <= end)
  
  # Extraer la columna S correspondiente y guardar en la lista
  segmentos_cortados[[ paste0("S", i) ]] <- profundidades[idx, c(paste0("S", i))]
}

#Calcular la longitud máxima
max_len <- max(sapply(segmentos_cortados, length))

# 2. Rellenar cada vector para que tenga la misma longitud
for (i in 1:length(segmentos_cortados)) {
  length(segmentos_cortados[[i]]) <- max_len
}

# Crear data frame final y guardar
profundidades_cortadas <- data.frame(POS = 1:max_len, segmentos_cortados)
write.table(profundidades_cortadas, "profundidades_cortadas.tsv", 
  quote = FALSE, sep = "\t", row.names = FALSE, eol = "\r", na = "")
