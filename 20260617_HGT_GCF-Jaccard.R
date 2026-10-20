library(ggplot2)
library(dplyr)
library(vegan)
library(ggpubr)

setwd("/mnt/c/Users/andre/OneDrive - UW-Madison/Projects/HGT_project")

MASH_MAX   <- 0.3      
MAX_POINTS <- 2000000000 

mash_raw <- read.table("20260608_AllVsAll_mash.txt", header = FALSE, sep = "\t",
                       row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
mash_mat <- as.matrix(mash_raw)
storage.mode(mash_mat) <- "numeric"
rownames(mash_mat) <- sub("\\.fna$", "", rownames(mash_mat))
colnames(mash_mat) <- rownames(mash_mat)

jac  <- readRDS("20260611_jac_distance.rds")
jac_mat <- as.matrix(jac)

common <- intersect(rownames(jac_mat), rownames(mash_mat))
cat("Genomes in common:", length(common), "\n")

#Align matrices
common <- intersect(rownames(jac_mat), rownames(mash_mat))
cat("Genomes in common:", length(common), "\n")
stopifnot(length(common) > 2)
jac_mat  <- jac_mat[common, common]
mash_mat <- mash_mat[common, common]

#pairwise values
common <- intersect(rownames(jac_mat), rownames(mash_mat))
cat("Genomes in common:", length(common), "\n")
stopifnot(length(common) > 2)
jac_mat  <- jac_mat[common, common]
mash_mat <- mash_mat[common, common]

lower <- lower.tri(jac_mat)
df <- data.frame(
  jaccard = jac_mat[lower],
  mash    = mash_mat[lower]
)


#Apply mash threshold
df <- df %>% filter(!is.na(mash), mash <= MASH_MAX)
cat("pairs kept (mash <=", MASH_MAX, "):", nrow(df), "\n")

#Log of Mash and Jaccard
df_log <- df
#replace 0s with 0.000001
#df_log <- df_log %>% 
#  mutate(across(where(is.numeric), ~ifelse(.x == 0, 0.00001, .x)))
#log
df_log$Mash_log10    <- log10(df$mash)
df_log$Jaccard_log10 <- log10(df$jaccard)
df_log <- subset(df_log, select = -c(mash, jaccard))
df_log <- rename(df_log, mash = Mash_log10, jaccard = Jaccard_log10)


#Subsample for plotting
set.seed(42)
df_plot <- if (nrow(df) > MAX_POINTS) df[sample(nrow(df), MAX_POINTS), ] else df



#Plot 
p <- ggplot(df_plot, aes(mash, jaccard)) +
  geom_bin2d(bins = 150) +
  scale_fill_viridis_c(trans = "log10", name = "pairs") +
  #geom_smooth(method = "loess", color="black", fill="blue", se = FALSE, linewidth = 0.7) + #loess for log, lm linear
  coord_cartesian(xlim = c(0, MASH_MAX)) +
  labs(title = sprintf("GCF (Jaccard) vs genome (Mash) distance — Mash <= %.2f", MASH_MAX),
       x = "Mash distance",
       y = "Jaccard distance (GCF composition)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
p


ggsave(sprintf("jaccard_vs_mash_max%.2f_fulldata.png", MASH_MAX),
       p, width = 24, height = 18, units = "cm")
ggsave(sprintf("jaccard_vs_mash_max%.2f_fulldata.pdf", MASH_MAX),
       p, width = 24, height = 18, units = "cm")


#Subsample for plotting log
set.seed(42)
df_plot <- if (nrow(df_log) > MAX_POINTS) df_log[sample(nrow(df_log), MAX_POINTS), ] else df_log


#Plot log
p <- ggplot(df_plot, aes(mash, jaccard)) +
  geom_bin2d(bins = 150) +
  scale_fill_viridis_c(trans = "log10", name = "pairs") +
  geom_smooth(method = "lm", color="black", fill="blue", se = FALSE, linewidth = 0.7) + #loess for log, lm linear
  stat_regline_equation(label.x=-2, label.y=-1.55) +
  stat_cor(aes(label=..rr.label..), label.x=-2, label.y=-1.6) +
  #coord_cartesian(xlim = c(-5, 0)) +
  labs(title = sprintf("GCF distance vs genome distance — Mash <= %.2f", MASH_MAX),
       x = "Log10 Mash distance",
       y = "Log10 Jaccard distance") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
p

ggsave(sprintf("jaccard_vs_mash_max%.2f_fulllog.png", MASH_MAX),
       p, width = 24, height = 18, units = "cm")
ggsave(sprintf("jaccard_vs_mash_max%.2f_fulllog.pdf", MASH_MAX),
       p, width = 24, height = 18, units = "cm")



