library(ggplot2)
library(dplyr)

# =========================
# INDIVIDUAL
# =========================

base_ind <- "Bep03"
ind_LC <- paste0(base_ind, "_1X")
ind_HC <- paste0(base_ind, "_30X")

# =========================
# PATHS
# =========================

base_dir <- "/media/daniel/Espacio/RefinedIBD"
path_LC <- file.path(base_dir, "LC")
path_HC <- file.path(base_dir, "HC")

output_file <- file.path(base_dir, paste0("IBD_mosaic_by_partner_", base_ind, ".tiff"))
output_table <- file.path(base_dir, paste0("IBD_mosaic_by_partner_", base_ind, ".csv"))

# =========================
# LENGTH CLASSES
# =========================

short_max_cm <- 3
long_min_cm <- 5

length_colors <- c(
  "Short (<3 cM)" = "#4E79A7",
  "Medium (3-5 cM)" = "#F28E2B",
  "Long (>=5 cM)" = "#D62728"
)

classify_length <- function(cm) {
  case_when(
    cm < short_max_cm ~ "Short (<3 cM)",
    cm < long_min_cm ~ "Medium (3-5 cM)",
    TRUE ~ "Long (>=5 cM)"
  )
}

clean_id <- function(x) {
  x %>%
    gsub("_1X", "", .) %>%
    gsub("_30X", "", .)
}

# =========================
# CHROMOSOMES AND GENOME OFFSETS
# =========================

cromosomas <- paste0("chr", 1:22)

chr_sizes <- c(
  chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
  chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
  chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
  chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
  chr17=83257441, chr18=80373285, chr19=58617616, chr20=64444167,
  chr21=46709983, chr22=50818468
)

chr_offsets <- c(0, cumsum(chr_sizes)[-length(chr_sizes)])
names(chr_offsets) <- names(chr_sizes)

chr_midpoints <- chr_offsets + chr_sizes / 2
chr_boundaries <- chr_offsets[-1]

# RefinedIBD columns:
# ID1 hap1 ID2 hap2 Chr Inicio Fin LOD cM
colnames_ibd <- c("ID1", "hap1", "ID2", "hap2", "Chr", "Inicio", "Fin", "LOD", "cM")

# =========================
# FUNCTIONS
# =========================

leer_chr <- function(path, chr) {
  archivo <- file.path(path, paste0(chr, "_ibd_merged.ibd"))
  df <- read.table(archivo, header = FALSE)
  colnames(df) <- colnames_ibd
  df
}

filtrar_individuo <- function(df, ind, dataset) {
  df %>%
    filter(ID1 == ind | ID2 == ind) %>%
    mutate(
      target = clean_id(ind),
      partner = ifelse(ID1 == ind, clean_id(ID2), clean_id(ID1)),
      dataset = dataset,
      length_class = classify_length(cM),
      x_start = chr_offsets[Chr] + Inicio,
      x_end = chr_offsets[Chr] + Fin
    )
}

leer_dataset <- function(path, ind, dataset) {
  bind_rows(lapply(cromosomas, function(chr) {
    leer_chr(path, chr) %>% filtrar_individuo(ind, dataset)
  }))
}

# =========================
# BUILD PLOT DATA
# =========================

LC <- leer_dataset(path_LC, ind_LC, "lcWGS")
HC <- leer_dataset(path_HC, ind_HC, "hcWGS")

plot_data <- bind_rows(LC, HC)

if(nrow(plot_data) == 0) {
  stop(paste("No IBD segments found for", base_ind))
}

partner_order <- plot_data %>%
  group_by(partner) %>%
  summarise(total_bp = sum(Fin - Inicio), .groups = "drop") %>%
  arrange(desc(total_bp), partner) %>%
  pull(partner)

row_levels <- as.vector(t(outer(partner_order, c("hcWGS", "lcWGS"), paste, sep = " - ")))

plot_data <- plot_data %>%
  mutate(
    row_label = paste(partner, dataset, sep = " - "),
    row_label = factor(row_label, levels = rev(row_levels)),
    length_class = factor(
      length_class,
      levels = c("Short (<3 cM)", "Medium (3-5 cM)", "Long (>=5 cM)")
    )
  )

write.csv(plot_data, output_table, row.names = FALSE)

# =========================
# PLOT
# =========================

grafico <- ggplot(plot_data) +
  geom_segment(
    aes(
      x = x_start,
      xend = x_end,
      y = row_label,
      yend = row_label,
      color = length_class
    ),
    linewidth = 1.8,
    lineend = "butt"
  ) +
  geom_vline(
    xintercept = chr_boundaries,
    color = "gray85",
    linewidth = 0.25
  ) +
  scale_color_manual(values = length_colors, drop = FALSE, name = "IBD length") +
  scale_x_continuous(
    breaks = chr_midpoints,
    labels = names(chr_midpoints),
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y = element_text(size = 6),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = paste0("IBD segments for ", base_ind),
    x = "Chromosome",
    y = "Partner - dataset"
  )

ggsave(
  output_file,
  grafico,
  width = 34,
  height = max(18, length(row_levels) * 0.28),
  units = "cm",
  dpi = 600,
  compression = "lzw"
)

cat("Segments plotted:", nrow(plot_data), "\n")
cat("Partners:", length(partner_order), "\n")
cat("Saved plot:", output_file, "\n")
cat("Saved table:", output_table, "\n")
