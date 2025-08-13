suppressPackageStartupMessages(library(tidyverse))

# Leer archivo
ruta <- "tables/READ_COUNTS.txt"
df <- read_tsv(ruta, col_types = cols(.default = "c"))

# Convertir columnas pertinentes a numéricas
df <- df %>%
  mutate(across(Reads:PairsAndWidows, ~ as.numeric(.)))

# Filtrar las filas que comienzan con "4-"
df4 <- df %>%
  filter(str_starts(Record, "4-")) 

# Agrupar por el prefijo antes del segundo guión
# Por ejemplo: "4-A_HA_H7" -> grupo = "4-A_HA"
# Quedarse con el máximo por cada grupo
# Extraer el nombre de segmento y ordenarlo como factor
df4_clean <- df4 %>%
  mutate(Segment = sub("^(4-A_[^_]+).*", "\\1", Record)) %>%
  group_by(Segment) %>%
  slice_max(order_by = Reads, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(Segment = str_remove(Segment, "^.*_")) %>% 
  mutate(Segment_num = case_when(
      Segment == "PB2" ~ 1,
      Segment == "PB1" ~ 2,
      Segment == "PA" ~ 3,
      Segment == "HA" ~ 4, 
      Segment == "NP" ~ 5,
      Segment == "NA" ~ 6,
      Segment == "MP" ~ 7,
      Segment == "NS" ~ 8
  )) %>% 
  mutate(Segment = fct_reorder(Segment, Segment_num)) %>% 
  arrange(Segment)
  
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

# Crear gráfico de barras
plot_reads <- df4_clean %>%
  ggplot(aes(x = Segment, y = Reads)) +
  geom_col(fill = my_colors) +
  geom_text(aes(label = scales::comma(Reads)), vjust = -0.5, size = 3.5) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0,0.1))) +
  labs(title = "Lecturas usadas en ensamblaje",
       x = "Segmento",
       y = "Número de lecturas") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
        plot.subtitle = element_text(size = 16),
        axis.text.x = element_text(size = 16, vjust = 0.5, hjust = 1, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "none", 
        axis.title = element_text(size = 20),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"))

ggsave(filename = "plot_reads.png", plot = plot_reads, device = "png", units = "in", width = 10, height = 7)