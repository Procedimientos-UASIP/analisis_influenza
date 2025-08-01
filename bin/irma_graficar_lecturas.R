suppressPackageStartupMessages(library(tidyverse))

# Leer archivo
ruta <- "tables/READ_COUNTS.txt"
df <- read_tsv(ruta, col_types = cols(.default = "c"))

# Convertir columnas numéricas
df <- df %>%
  mutate(across(Reads:PairsAndWidows, ~ as.numeric(.)))

# Filtrar las filas que comienzan con "4-"
df4 <- df %>%
  filter(str_starts(Record, "4-"))

# Agrupar por el prefijo antes del segundo guión
# Por ejemplo: "4-A_HA_H7" -> grupo = "4-A_HA"
df4 <- df4 %>%
  mutate(grupo = sub("^(4-A_[^_]+).*", "\\1", Record)) %>%
  group_by(grupo) %>%
  slice_max(order_by = Reads, n = 1, with_ties = FALSE) %>%
  ungroup()

df4 <- df4 %>%
  mutate(Record = factor(Record, levels = Record))

# Crear gráfico de barras
plot_reads <- df4 %>%
  ggplot(aes(x = Record, y = Reads)) +
   geom_col(fill = "#2a9d8f") +
   geom_text(aes(label = Reads), vjust = -0.5, size = 3.5) +
   labs(title = "Lecturas por segmento viral",
        x = "Segmento",
        y = "Número de lecturas") +
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "plot_reads.png", plot = plot_reads, device = "png", units = "in", width = 10, height = 7)