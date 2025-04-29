#!/usr/bin/env Rscript

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
  make_option(c("--input_file"), dest = "input_file", type = "character", help = "Archivo de datos (ej: coberturas_finales.tsv)"),
  make_option(c("--muestra"), dest = "muestra", type = "character", help = "Nombre de la muestra a procesar")
)

# Parsear los argumentos
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Verificar que ambos argumentos fueron suministrados
if (is.null(opt$input_file) || is.null(opt$muestra)) {
  stop("Faltan argumentos. Debes proporcionar el archivo de datos y el nombre de la muestra. Ejemplo: ./mi_script.R --in coberturas_finales.tsv --muestra CPA-XXXX-2025", call. = FALSE)
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

# Leer datos y renombrar
df_reads <- suppressMessages(read_tsv(input_file, show_col_types = FALSE))

# Confirmación de lectura y datos
message("✅ Archivo leído correctamente: ", input_file)
message("- Número de filas: ", nrow(df_reads))
message("- Número de columnas: ", ncol(df_reads))
message("- Muestra procesada: ", sample_name)
message("💻 Limpiando y procesando")

# Calcular cuántas lecturas alineadas no se usaron en los ensambles
Lecturas_no_usadas = sum(df_reads$Lecturas_alineadas) - sum(df_reads$Lecturas_ensamblaje) 

# Limpiar df, agregar lecturas no usadas, y crear etiquetas
df_limpio <- df_reads %>% 
  select(Segmento, Lecturas_alineadas) %>% 
  rename(Lecturas = Lecturas_alineadas) %>% 
  add_row(Segmento = "No usadas", Lecturas = Lecturas_no_usadas) %>% 
  mutate(Porcentaje = Lecturas/sum(Lecturas) * 100,
         etiqueta = sprintf("%s\n(%.2f%%)", scales::comma(Lecturas), Porcentaje))

# Definir la relación entre Lecturas y Porcentaje
# Por ejemplo, suponiendo que 100% = suma total de lecturas
total_lecturas <- sum(df_limpio$Lecturas)
factor_conversion <- total_lecturas / 100

# Define los breaks del eje secundario (porcentajes deseados)
breaks_secundario <- seq(0, 30, by = 5)
breaks_principal <- breaks_secundario * factor_conversion

# Graficar
p1 <- df_limpio %>% 
  ggplot(aes(x = Segmento , y = Lecturas, fill = Segmento)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = etiqueta), vjust = -0.5, size = 5) +
    labs(title = "Lecturas usadas para el ensamblaje de\nlos segmentos de Influenza tipo A",
    subtitle = paste0( "Muestra: ",sample_name)) +
  scale_y_continuous(
    name = "Lecturas",
    labels = scales::comma,
    sec.axis = sec_axis(~ . / factor_conversion, name = "Porcentaje (%)", breaks = breaks_secundario), 
    expand = expansion(mult = c(0,0.15))) +
  scale_fill_manual(values = pals::brewer.set1(9)) + 
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold"), 
        plot.subtitle = element_text(size = 16),
        axis.title = element_text(size = 22, face = "bold"),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none")

# Guardar plot
ggsave(filename = paste0(sample_name, "_lecturas.png"), plot = p1, device = "png", units = "cm", width = 28, height = 22)

message("📄 Gráfica producida. Finalizando")