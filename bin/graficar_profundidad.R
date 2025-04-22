#!/usr/bin/env Rscript

# ---------------------- #
# 1. Verificación de paquetes
# ---------------------- #

# Lista de paquetes requeridos
required_packages <- c("tidyverse", "scales", "ggrepel", "pals", "optparse")

# Función para verificar los paquetes
check_required_packages <- function(packages) {

  # Verifica cuáles no están instalados
  not_installed <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
  
  # Si hay paquetes no instalados, detiene la ejecución
  if (length(not_installed) > 0) {
    stop(
      paste(
        "Los siguientes paquetes no están instalados:",
        paste(not_installed, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

# Ejecutar la verificación
check_required_packages(required_packages)

# Cargar paquetes requeridos. No se cargan todos para evitar enmascaramientos
for (pkg in c("tidyverse", "ggrepel", "optparse")) {
  suppressMessages(library(pkg, character.only = TRUE))
}

# ---------------------- #
# 2. Definición de opciones con optparse
# ---------------------- #

option_list <- list(
  make_option(c("--input_file"), dest = "input_file", type = "character", help = "Archivo de datos (ej: coberturas_finales.tsv)"),
  make_option(c("--muestra"), dest = "muestra", type = "character", help = "Nombre de la muestra a procesar")
)

# Parsear los argumentos
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Verificar que ambos argumentos fueron suministrados
if (is.null(opt$input_file) || is.null(opt$muestra)) {
  stop("Faltan argumentos. Debes proporcionar el archivo de datos y el nombre de la muestra. Ejemplo: ./mi_script.R --in coberturas_finales.tsv --muestra muestra_1", call. = FALSE)
}

# Asignar a variables
input_file <- opt$input_file
sample_name <- opt$muestra

# ---------------------- #
# 3. Revisar el archivo
# ---------------------- #

if (!file.exists(input_file)) {
  stop(paste("El archivo especificado no existe:", input_file), call. = FALSE)
}

# ---------------------- #
# 4. SCRIPT
# ---------------------- #

# Puedes usar read_tsv de readr (parte de tidyverse)
df_prof <- suppressMessages(read_tsv(input_file, show_col_types = FALSE) %>%
  rename("Posición" = "POS",
         "S1(PB2)" = "S1",
         "S2(PB1)" = "S2",
         "S3(PA)" = "S3",
         "S4(HA)" = "S4",
         "S5(NP)" = "S5",
         "S6(NA)" = "S6",
         "S7(M1, M2)" = "S7",
         "S8(Ns1, Ns2)" = "S8"))

# Confirmación de lectura y datos
message("Archivo leído correctamente: ", input_file)
message("Número de filas: ", nrow(df_prof))
message("Número de columnas: ", ncol(df_prof))
message("Muestra procesada: ", sample_name)

# Transformar, ordenar y filtrar
df_prof_long <- df_prof %>% 
  pivot_longer(names_to = "Segmento", values_to = "Profundidad", cols = 2:9) %>% 
  arrange(Segmento) %>% 
  filter(!is.na(Profundidad))

# Graficar
plot_prof <- df_prof_long %>% 
  ggplot(aes(x = Posición, y = Profundidad, color = Segmento))+
  facet_wrap(~ Segmento, ncol = 1, strip.position = "right", scales = "free_y") +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = pals::brewer.set1(length(unique(df_prof_long$Segmento)))) +
  labs(title = "Profundidad de secuenciación de los segmentos\nde Influenza tipo A",
       subtitle = paste0( "Muestra: ",sample_name)) +
  scale_x_continuous(breaks = seq (0,max(df_prof_long$Posición), 100), labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold"), 
        plot.subtitle = element_text(size = 16),
        axis.text.x = element_text(color = "black", size = 10, angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(color = "black", size = 10),
        legend.position = "none", 
        axis.title = element_text(size = 20),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"))

# Guardar plot
ggsave(filename = paste0(sample_name, "_profundidades.png"), plot = plot_prof, device = "png", units = "cm", width = 20, height = 26)

message("Gráfica producidad. Finalizando")