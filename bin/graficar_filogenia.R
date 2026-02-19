#!/usr/bin/env Rscript

# ---------------------- #
# 1. VERIFICACIÓN DE LIBRERIAS
# ---------------------- #

# Lista de paquetes requeridos
required_packages <- c("Biostrings", "tidyverse", "ggtree", "optparse", "phangorn")

# Función para verificar los paquetes
check_required_packages <- function(packages) {

  # Verifica cuáles no están instalados
  not_installed <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
  
  # Si hay paquetes no instalados, detiene la ejecución
  if (length(not_installed) > 0) {
    stop(
      paste( "Los siguientes paquetes no están instalados:", paste(not_installed, collapse = ", ")), call. = FALSE
      )
  }
}

# Ejecutar la verificación
check_required_packages(required_packages)

# Cargar paquetes requeridos. No se cargan todos para evitar enmascaramientos
# No se carga phangorn para no enmascarar rotate.
for (pkg in c("Biostrings", "tidyverse", "ggtree", "optparse")) {
  suppressMessages(library(pkg, character.only = TRUE))
}

# ---------------------- #
# 2. DEFINIR OPCIONES CON OPTPARSE
# ---------------------- #

option_list <- list(
  make_option(c("--tree_HA"), dest = "tree_HA", type = "character", help = "Archivo de filogenia con bootstrap (ej: CPA-XXXXX-XX_HA_tree.contree)"),
  make_option(c("--tree_NA"), dest = "tree_NA", type = "character", help = "Archivo de filogenia con bootstrap (ej: CPA-XXXXX-XX_NA_tree.contree)"),
  make_option(c("--aln_HA"), dest = "aln_HA", type = "character", help = "Archivo de alineamiento (ej: NCBI_virus_HA.aln)"),
  make_option(c("--aln_NA"), dest = "aln_NA", type = "character", help = "Archivo de alineamiento (ej: NCBI_virus_NA.aln)"),
  make_option(c("--muestra"), dest = "muestra", type = "character", help = "Nombre de la muestra a procesar")
)

# Parsear los argumentos
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Verificar que ambos argumentos fueron suministrados
if (is.null(opt$tree_HA) || is.null(opt$tree_NA) || is.null(opt$aln_HA) || is.null(opt$aln_NA) || is.null(opt$muestra)) {
  stop("Faltan argumentos. Debes proporcionar los archivos de filogenia con bootstrap, los alineamientos con la muestra problema, y el nombre de la muestra. 
  Ejemplo: ./graficar_filogenia.R 
  --tree_HA CPA-XXXXX-XX_HA.contree --tree_NA CPA-XXXXX-XX_NA.contree 
  --aln_HA NCBI_virus_HA.aln --aln_NA NCBI_virus_NA.aln
  --muestra CPA-XXXX-XX", call. = FALSE)
}

# Asignar a variables
tree_HA <- opt$tree_HA
tree_NA <- opt$tree_NA
aln_HA <- opt$aln_HA
aln_NA <- opt$aln_NA

sample_name <- opt$muestra

# ---------------------- #
# 3. REVISAR ARCHIVO
# ---------------------- #

if (!file.exists(tree_HA)) {
  stop(paste("El archivo especificado no existe:", tree_HA), call. = FALSE)
}

if (!file.exists(tree_NA)) {
  stop(paste("El archivo especificado no existe:", tree_NA), call. = FALSE)
}

if (!file.exists(aln_HA)) {
  stop(paste("El archivo especificado no existe:", aln_HA), call. = FALSE)
}

if (!file.exists(aln_NA)) {
  stop(paste("El archivo especificado no existe:", aln_NA), call. = FALSE)
}

# ---------------------- #
# 4. SCRIPT 
# ---------------------- #

# ---------------------- #
# 4.1 PARA HA 
# ---------------------- #

# Leer el archivo con coverturas
aln <- readDNAStringSet(aln_HA)
tree <- read.tree(tree_HA)
tree_mid <- phangorn::midpoint(tree)

# Confirmación de lectura y datos
message("✅ Archivo de filogenias de HA leído correctamente: ", tree_HA)
message("✅ Archivo de alineamiento de HA leído correctamente: ", aln_HA)
message("- Muestra procesada: ", sample_name)
message("💻 Dibujando árbol filogenlético de HA con muestra problema")

# Extraer nombre de headers
seq_headers <- names(aln)

# Construir df de metadados. Muestra problema queda como "NA" en HA_subtype 
metadata <- data.frame(seq_headers = seq_headers) %>% 
  separate(seq_headers, 
  into = c("accesion", "HA_subtype", "host", "country"),
  sep = "\\|", 
  remove = FALSE) %>%
  select(seq_headers, HA_subtype) %>% 
  mutate(HA_subtype = ifelse(is.na(HA_subtype), "NA", as.character(HA_subtype)))

  # Paleta de colores segun el subtipo de HA. Muestra problema en negro.
ha_colors <- c(
  H1 = "#e41a1c",
  H2 = "#377eb8",
  H3 = "#4daf4a",
  H4 = "#984ea3",
  H5 = "#ff7f00",
  H6 = "#ffff33",
  H7 = "#a65628",
  H8 = "#f781bf",
  H9 = "#999999",
  H10 = "#66c2a5",
  H11 = "#fc8d62",
  H12 = "#8da0cb",
  H13 = "#e78ac3",
  "NA" = "#000000"
)

# Dibujar el arbol
p <- ggtree(tree_mid, layout = "fan", open.angle = 180) %<+% metadata +
    geom_tippoint(aes(fill = HA_subtype), size = 5, shape = 21, color = "black") +
    geom_tree(size = 1) +
    geom_tiplab2(
      aes(color = ifelse(HA_subtype == "NA", "NA", "Other")),
      size = 2.8,
      offset = 0.05,
      align = TRUE,
      linetype = "dashed",
      linesize = 0.4
    ) +
    scale_fill_manual(values = ha_colors) + 
    scale_color_manual(
      values = c(
        "NA" = "red",
        "Other" = "black"
      ),
      guide = "none"
    ) +
    theme(plot.margin = margin(t = 120),
          legend.position = "bottom",
          legend.key.size = unit(1, "cm"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 20),
          legend.box.margin = margin(t = -350))

ggsave(paste0(sample_name, "_HA_tree.pdf"), plot = p, device = "pdf", width = 12, height = 8, units = "in")


# ---------------------- #
# 4.2 PARA NA 
# ---------------------- #

# Leer el archivo con coverturas
aln <- readDNAStringSet(aln_NA)
tree <- read.tree(tree_NA)
tree_mid <- phangorn::midpoint(tree)

# Confirmación de lectura y datos
message("✅ Archivo de filogenias de NA leído correctamente: ", tree_NA)
message("✅ Archivo de alineamiento de NA leído correctamente: ", aln_NA)
message("- Muestra procesada: ", sample_name)
message("💻 Dibujando árbol filogenlético de NA con muestra problema")

# Extraer nombre de headers
seq_headers <- names(aln)

# Construir df de metadados. Muestra problema queda como "NA" en HA_subtype 
metadata <- data.frame(seq_headers = seq_headers) %>% 
  separate(seq_headers, 
           into = c("accesion", "NA_subtype", "host", "country"),
           sep = "\\|", 
           remove = FALSE) %>%
  select(seq_headers, NA_subtype) %>% 
  mutate(NA_subtype = ifelse(is.na(NA_subtype), "NA", as.character(NA_subtype)))


# Paleta de colores segun el subtipo de NA. Muestra problema en negro.
na_colors <- c(
  N1 = "#e41a1c",
  N2 = "#377eb8",
  N3 = "#4daf4a",
  N4 = "#984ea3",
  N5 = "#ff7f00",
  N6 = "#ffff33",
  N7 = "#a65628",
  N8 = "#f781bf",
  N9 = "#999999",
  "NA" = "#000000"
)

# Dibujar el arbol
p <- ggtree(tree_mid, layout = "fan", open.angle = 180) %<+% metadata +
  geom_tippoint(aes(fill = NA_subtype), size = 5, shape = 21, color = "black") +
  geom_tree(size = 1) +
  geom_tiplab2(
    aes(color = ifelse(NA_subtype == "NA", "NA", "Other")),
    size = 2.8,
    offset = 0.05,
    align = TRUE,
    linetype = "dashed",
    linesize = 0.4
  ) +
  scale_fill_manual(values = na_colors) + 
  scale_color_manual(
    values = c(
      "NA" = "red",
      "Other" = "black"
    ),
    guide = "none"
  ) +
  theme(plot.margin = margin(t = 120),
        legend.position = "bottom",
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 20),
        legend.box.margin = margin(t = -350))

ggsave(paste0(sample_name, "_NA_tree.pdf"), plot = p, device = "pdf", width = 13, height = 8, units = "in")
