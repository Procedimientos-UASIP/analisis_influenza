#!/usr/bin/env Rscript

# ---------------------- #
# 1. VERIFICACIÓN DE LIBRERIAS
# ---------------------- #

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

# Lista de paquetes requeridos
required_packages <- c("tidyverse", "scales", "pals", "optparse")

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
  make_option(c("--input_file"), dest = "input_file", type = "character", help = "Archivo de datos (ej: coberturas_finales.tsv)"),
  make_option(c("--muestra"), dest = "muestra", type = "character", help = "Nombre de la muestra a procesar")
)

# Parsear los argumentos
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Verificar que ambos argumentos fueron suministrados.
if (is.null(opt$input_file) || is.null(opt$muestra)) {
  stop("Faltan argumentos. Debes proporcionar el archivo de datos y el nombre de la muestra. Ejemplo: ./mi_script.R --in coberturas_finales.tsv --muestra CPA-XXXX-2025", call. = FALSE)
}

# Asignar a variables.
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

# Leer datos.
df_reads <- suppressMessages(read_tsv(input_file, show_col_types = FALSE))

# Confirmación de lectura y datos.
message("✅ Archivo leído correctamente: ", input_file)
message("- Muestra procesada: ", sample_name)
message("💻 Limpiando y procesando")

# Calcular cuántas lecturas alineadas no se usaron en los ensambles.
Lecturas_no_usadas = sum(df_reads$Lecturas_alineadas) - sum(df_reads$Lecturas_ensamblaje) 

# Calcular lecturas no ensambladas por segmento
df_limpio <- df_reads %>%
  mutate(Lecturas_no_usadas = Lecturas_alineadas - Lecturas_ensamblaje) %>%
  pivot_longer(
    cols = c(Lecturas_ensamblaje, Lecturas_no_usadas),
    names_to = "Tipo",
    values_to = "Lecturas"
  ) %>%
  mutate(Tipo = recode(Tipo, "Lecturas_ensamblaje" = "Usadas", "Lecturas_no_usadas" = "No usadas"))

# Calcular totales por segmento para porcentajes
df_totales <- df_limpio %>%
  group_by(Segmento) %>%
  summarise(Total = sum(Lecturas), .groups = "drop")

# Unir los totales y calcular porcentaje de cada componente
df_limpio_2 <- df_limpio %>%
  left_join(df_totales, by = "Segmento") %>%
  mutate(Porcentaje = (Lecturas / Total) * 100,
         etiqueta = if_else(Tipo == "Usadas",
                       sprintf("%s\n(%.2f%%)", scales::comma(Lecturas), Porcentaje),
                       NA_character_)  # Solo etiquetar parte usada
  )

# Calcular sumatorias totales
total_alineadas <- sum(df_reads$Lecturas_alineadas)
total_ensambladas <- sum(df_reads$Lecturas_ensamblaje)
porcentaje_total <- (total_ensambladas / total_alineadas) * 100

# Crear subtítulo dinámico
subtitulo <- sprintf("Muestra: %s\n%s / %s (Alineadas / Ensambladas) (%.2f%%)",
  sample_name,
  scales::comma(total_alineadas),
  scales::comma(total_ensambladas),
  porcentaje_total)

# Graficar
p1 <- df_limpio_2 %>% 
  ggplot(aes(x = Segmento , y = Lecturas, fill = Tipo)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(data = df_limpio_2 %>% filter(!is.na(etiqueta)),
    aes(label = etiqueta, y = Lecturas_alineadas),
    vjust = -0.2,  # ajusta verticalmente
    size = 5
  ) +
  labs(
    title = "Porcentajes de uso de lecturas en ensambles",
    subtitle = subtitulo,
    x = "Segmento"
  ) +
  scale_y_continuous(
    name = "Lecturas",
    breaks = seq(0,1000000, 25000),
    labels = scales::comma,
    #sec.axis = sec_axis(~ . / factor_conversion, name = "Porcentaje (%)", breaks = breaks_secundario), 
    expand = expansion(mult = c(0,0.15))) +
  scale_fill_manual(
    values = c("Usadas" = "#4DAF4A", "No usadas" = "#E41A1C"),
    name = "Tipo de lectura"
  ) + 
  theme_bw() +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 18),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 16, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )

# Guardar plot
ggsave(filename = paste0(sample_name, "_lecturas.png"), plot = p1, device = "png", units = "cm", width = 28, height = 22)

message("📄 Gráfica producida. Finalizando")