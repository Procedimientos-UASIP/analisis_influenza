# ---------------------- #
# 1. VERIFICACIÓN DE LIBRERIAS
# ---------------------- #

# Lista de paquetes requeridos
required_packages <- c("tidyverse", "scales", "pals", "optparse")

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
for (pkg in c("tidyverse", "optparse")) {
  suppressMessages(library(pkg, character.only = TRUE))
}

# ---------------------- #
# 2. DEFINIR OPCIONES CON OPTPARSE
# ---------------------- #

option_list <- list(
  make_option(c("--input_file"), dest = "input_file", type = "character", help = "Archivo de datos (ej: profundidades_finales.tsv)"),
  make_option(c("--muestra"), dest = "muestra", type = "character", help = "Nombre de la muestra a procesar")
)

# Parsear los argumentos
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Verificar que ambos argumentos fueron suministrados
if (is.null(opt$input_file) || is.null(opt$muestra)) {
  stop("Faltan argumentos. Debes proporcionar el archivo de datos y el nombre de la muestra. Ejemplo: ./mi_script.R --in profundidades_finales.tsv --muestra CPA-XXXX-2025", call. = FALSE)
}

# Asignar a variables
input_file <- opt$input_file
sample_name <- opt$muestra

# ---------------------- #
# 3. REVISAR ARCHIVO
# ---------------------- #

if (!file.exists(input_file)) {
  stop(paste("El archivo especificado no existe:", input_file), call. = FALSE)
}

# ---------------------- #
# 4. SCRIPT
# ---------------------- #

# Leer el archivo con coverturas
df_prof <- suppressMessages(read_tsv(input_file, show_col_types = FALSE))

# Confirmación de lectura y datos
message("✅ Archivo leído correctamente: ", input_file)
message("- Número de filas: ", nrow(df_prof))
message("- Número de columnas: ", ncol(df_prof))
message("- Muestra procesada: ", sample_name)
message("💻 Limpiando y procesando")

# Definir las columnas esperadas
#cols_esperadas <- paste0("S", 1:8)

# Añadir columnas faltantes con NA
#for (col in cols_esperadas) {
#  if (!col %in% names(df_prof)) {
#    df_prof[[col]] <- NA
#  }
#}

# Reordenar columnas y renombrar
df_prof <- df_prof %>% 
  rename("Posición" = "Posicion",
         "S1(PB2)" = 2,
         "S2(PB1)" = 3,
         "S3(PA)" = 4,
         "S4(HA)" = 5,
         "S5(NP)" = 6,
         "S6(NA)" = 7,
         "S7(M1, M2)" = 8,
         "S8(Ns1, Ns2)" = 9)

# Transformar, ordenar y filtrar
df_prof_long <- df_prof %>% 
  pivot_longer(names_to = "Segmento", values_to = "Profundidad", cols = 2:9) %>% 
  arrange(Segmento) %>% 
  filter(!is.na(Profundidad))

# Calcular promedio de profundidades y darles formato
df_avg <- df_prof_long %>%
  group_by(Segmento) %>%
  summarise(profundidad_promedio = mean(Profundidad, na.rm = TRUE)) %>%
  mutate(
    profundidad_formateada = scales::comma_format(accuracy = 0.1)(profundidad_promedio),
    titulo = paste0(Segmento, ". Profundidad promedio = ", profundidad_formateada))

# Sustituir etiquedas
df_prof_long_labeled <- df_prof_long %>%
  left_join(df_avg, by = "Segmento") %>%
  mutate(Segmento = factor(titulo, levels = unique(titulo)))

my_colors <- c(
  "#1b9e77", # verde oscuro
  "#d95f02", # naranja óxido
  "#7570b3", # púrpura medio
  "#e7298a", # magenta oscuro
  "#66a61e", # verde oliva
  "#e6ab02", # mostaza
  "#a6761d", # marrón dorado
  "#666666"  # gris oscuro
)

# Graficar
p1 <- df_prof_long_labeled  %>% 
  ggplot(aes(x = Posición, y = Profundidad, color = Segmento))+
  facet_wrap(~ Segmento, ncol = 2, strip.position = "top", scales = "free_y") +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = my_colors) +
  labs(title = "Profundidad de secuenciación de segmentos de IA-A",
       subtitle = paste0( "Muestra: ",sample_name)) +
  scale_x_continuous(breaks = seq (0,max(df_prof_long_labeled$Posición), 100), labels = scales::comma) +
  scale_y_continuous(labels = scales::comma, limits = c(1, NA)) +
  theme_bw() +
  theme(plot.title = element_text(size = 20), 
        plot.subtitle = element_text(size = 16),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1), 
        axis.text.y = element_text(size = 8),
        legend.position = "none", 
        axis.title = element_text(size = 20),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"))

# Guardar plot
ggsave(filename = paste0(sample_name, "_prof.png"), plot = p1, device = "png", units = "cm", width = 24, height = 14)

message("📄 Gráfica producida. Finalizando")