suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(circlize))

base_ind <- "Bep03"
partner_filter <- "LR174"
min_cm <- 0

base_dir <- "/media/daniel/Espacio/RefinedIBD"
path_LC <- file.path(base_dir, "LC")
path_HC <- file.path(base_dir, "HC")

output_file <- file.path(base_dir, paste0("IBD_circus_", base_ind, "_vs_", partner_filter, ".tiff"))
output_table <- file.path(base_dir, paste0("IBD_circus_", base_ind, "_vs_", partner_filter, ".csv"))

ind_LC <- paste0(base_ind, "_1X")
ind_HC <- paste0(base_ind, "_30X")

chromosomes <- paste0("chr", 1:22)

chr_sizes <- c(
  chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
  chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
  chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
  chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
  chr17=83257441, chr18=80373285, chr19=58617616, chr20=64444167,
  chr21=46709983, chr22=50818468
)

colnames_ibd <- c("ID1", "hap1", "ID2", "hap2", "Chr", "Inicio", "Fin", "LOD", "cM")

clean_id <- function(x) {
  x %>% gsub("_1X", "", .) %>% gsub("_30X", "", .)
}

length_class <- function(cm) {
  case_when(
    cm < 3 ~ "Short (<3 cM)",
    cm < 5 ~ "Medium (3-5 cM)",
    TRUE ~ "Long (>=5 cM)"
  )
}

read_chr <- function(path, chr) {
  archivo <- file.path(path, paste0(chr, "_ibd_merged.ibd"))
  df <- read.table(archivo, header = FALSE)
  colnames(df) <- colnames_ibd
  df
}

read_dataset <- function(path, ind, dataset) {
  bind_rows(lapply(chromosomes, function(chr) {
    df <- read_chr(path, chr) %>%
      filter(ID1 == ind | ID2 == ind) %>%
      filter(clean_id(ifelse(ID1 == ind, ID2, ID1)) == partner_filter) %>%
      filter(cM >= min_cm)

    if(nrow(df) == 0) {
      return(data.frame())
    }

    df %>%
      mutate(
        target = clean_id(ind),
        partner = ifelse(ID1 == ind, clean_id(ID2), clean_id(ID1)),
        dataset = dataset,
        class = length_class(cM)
      )
  }))
}

plot_data <- bind_rows(
  read_dataset(path_HC, ind_HC, "hcWGS"),
  read_dataset(path_LC, ind_LC, "lcWGS")
)

if(nrow(plot_data) == 0) {
  stop(paste("No segments >=", min_cm, "cM found for", base_ind))
}

plot_data <- plot_data %>%
  mutate(
    class = factor(class, levels = c("Short (<3 cM)", "Medium (3-5 cM)", "Long (>=5 cM)"))
  )

write.csv(plot_data, output_table, row.names = FALSE)

colors <- c(
  "Short (<3 cM)" = "#4E79A7",
  "Medium (3-5 cM)" = "#F28E2B",
  "Long (>=5 cM)" = "#D62728"
)

plot_track <- function(df, label) {
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = 0.12,
    bg.border = "gray85",
    panel.fun = function(x, y) {
      chr <- CELL_META$sector.index
      chr_df <- df[df$Chr == chr, , drop = FALSE]

      if(nrow(chr_df) > 0) {
        for(i in seq_len(nrow(chr_df))) {
          circos.rect(
            chr_df$Inicio[i], 0.12,
            chr_df$Fin[i], 0.88,
            col = colors[as.character(chr_df$class[i])],
            border = NA
          )
        }
      }

    }
  )
}

tiff(output_file, width = 18, height = 18, units = "cm", res = 600, compression = "lzw")

circos.clear()
circos.par(
  start.degree = 90,
  gap.after = c(rep(1.5, 21), 6),
  cell.padding = c(0, 0, 0, 0),
  track.margin = c(0.006, 0.006)
)

chr_df <- data.frame(
  chr = names(chr_sizes),
  start = 0,
  end = as.numeric(chr_sizes)
)

circos.initialize(factors = chr_df$chr, xlim = chr_df[, c("start", "end")])

circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.08,
  bg.col = "gray95",
  bg.border = "gray70",
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    circos.text(
      CELL_META$xcenter,
      0.5,
      gsub("chr", "", chr),
      facing = "clockwise",
      niceFacing = TRUE,
      cex = 0.55
    )
  }
)

plot_track(plot_data %>% filter(dataset == "hcWGS"), "hcWGS")
plot_track(plot_data %>% filter(dataset == "lcWGS"), "lcWGS")

legend(
  "bottomleft",
  legend = names(colors),
  fill = colors,
  border = NA,
  bty = "n",
  cex = 0.75,
  title = "IBD length"
)

legend(
  "bottomright",
  legend = c("outer track: hcWGS", "inner track: lcWGS"),
  lwd = 3,
  col = c("gray35", "gray35"),
  bty = "n",
  cex = 0.75,
  title = "Tracks"
)

title(paste0("IBD: ", base_ind, " vs ", partner_filter), cex.main = 1.1)

dev.off()

cat("Segments plotted:", nrow(plot_data), "\n")
cat("Partner:", partner_filter, "\n")
cat("Minimum cM:", min_cm, "\n")
cat("Saved plot:", output_file, "\n")
cat("Saved table:", output_table, "\n")
