#Load libraries
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Matrix)
library(data.table)

setwd("/mnt/c/Users/andre/OneDrive - UW-Madison/Projects/HGT_project")

#Load files
bgc <- fread("20260608_HGT_BiGSLICE_T04_parsedv2.csv",select = c("Genome", "gcf_id"), data.table=TRUE)
meta <- fread("20260610_Taxonomy_HighQActinos.csv", data.table=TRUE)

COLOR_BY    <- "Order" 
MAX_LEGEND  <- 40 
TOP_N       <- 20          # how many top categories to show
SHOW_OTHER  <- TRUE

#Make presence absence of gcfs vs Genomes
pa <- bgc %>%
  select(Genome, gcf_id) %>%
  filter(!is.na(Genome), !is.na(gcf_id)) %>%
  distinct()

# Sparse incidence matrix (rows = genomes, cols = GCFs), values = 1.
g_levels <- sort(unique(pa$Genome)) #genome IDs
f_levels <- sort(unique(pa$gcf_id)) #gcf_ids
mat <- sparseMatrix(
  i = match(pa$Genome, g_levels),
  j = match(pa$gcf_id, f_levels),
  x = 1,
  dims = c(length(g_levels), length(f_levels)),
  dimnames = list(g_levels, f_levels)
)
mat@x[mat@x > 1] <- 1   # ensure binary

cat("Matrix:", nrow(mat), "genomes x", ncol(mat), "GCFs\n")

#Jaccard Dist and PCoA
#Jaccard on presence absence
mat@x[mat@x > 1] <- 1
mat <- drop0(mat)
inter <- as.matrix(tcrossprod(mat))   # n x n, ~7.6 GB
a     <- rowSums(mat)
union <- outer(a, a, "+") - inter
sim   <- inter / union
sim[!is.finite(sim)] <- 0
diag(sim) <- 1
jac   <- as.dist(1 - sim)
rm(inter, union, sim); gc()

#saveRDS(jac, "20260611_jac_distance.rds")



#PCoA
#K is top five principal coordinate axis
pcoa <- cmdscale(jac, k = 5, eig = TRUE, add = TRUE)
#saveRDS(pcoa, "20260611_pcoa_result_jaccard.rds")





#Reload jaccard and PCoA
#Start from here if data has been calculated

jac  <- readRDS("20260611_jac_distance.rds")
pcoa <- readRDS("20260611_pcoa_result_jaccard.rds")
#write.csv(as.matrix(jac), "20260611_jac_matrix.csv")

#% Variation in each axis
eig_pos <- pcoa$eig[pcoa$eig > 0]
var_exp <- 100 * eig_pos / sum(eig_pos)
cat(sprintf("PCo1 = %.1f%%, PCo2 = %.1f%%\n", var_exp[1], var_exp[2]))

#Merge coordinates with metadata
scores <- data.frame(Genome = rownames(pcoa$points),
                     PCo1 = pcoa$points[, 1],
                     PCo2 = pcoa$points[, 2],
                     stringsAsFactors = FALSE)
scores <- left_join(scores, meta, by = "Genome")

scores[[COLOR_BY]][is.na(scores[[COLOR_BY]]) | scores[[COLOR_BY]] == ""] <- "Unclassified"
top_levels <- scores %>% count(.data[[COLOR_BY]], sort = TRUE) %>%
  slice_head(n = MAX_LEGEND) %>% pull(1)
scores$plot_group <- ifelse(scores[[COLOR_BY]] %in% top_levels,
                            scores[[COLOR_BY]], "Other")


ggplot(scores, aes(PCo1, PCo2, color = plot_group)) +
  geom_point(size = 2, alpha = 0.6) +
  labs(title = "PCoA of GCF composition (Jaccard)",
       x = sprintf("PCo1 (%.1f%%)", var_exp[1]),
       y = sprintf("PCo2 (%.1f%%)", var_exp[2]),
       color = COLOR_BY) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

ggsave("20260611_HGT_PCoA_order.png", device = "png", width = 40, height = 30, units = "cm")
ggsave("20260611_HGT_PCoA_order.pdf", device = "pdf", width = 40, height = 30, units = "cm")


#Plot only top n
COLOR_BY    <- "Class"  
TOP_N       <- 10          # how many top categories to show
SHOW_OTHER  <- TRUE

scores[[COLOR_BY]][is.na(scores[[COLOR_BY]]) | scores[[COLOR_BY]] == ""] <- "Unclassified"

#Rank
top_levels <- scores %>%
  count(.data[[COLOR_BY]], sort = TRUE) %>%   # most frequent first
  slice_head(n = TOP_N) %>%
  pull(1)

#Make non top into Other category
scores$plot_group <- ifelse(scores[[COLOR_BY]] %in% top_levels,
                            scores[[COLOR_BY]], "Other")

# Order the legend by frequency
scores$plot_group <- factor(scores$plot_group,
                            levels = c(top_levels, "Other"))
plot_df <- if (SHOW_OTHER) scores else dplyr::filter(scores, plot_group != "Other")

p <- ggplot() +
  {if (SHOW_OTHER)
    geom_point(data = dplyr::filter(plot_df, plot_group == "Other"),
               aes(PCo1, PCo2), color = "grey80", size = 1.5, alpha = 0.4)} +
  geom_point(data = dplyr::filter(plot_df, plot_group != "Other"),
             aes(PCo1, PCo2, color = plot_group), size = 2, alpha = 0.7) +
  labs(title = sprintf("PCoA of GCF composition (Jaccard) — top %d %s", TOP_N, COLOR_BY),
       x = sprintf("PCo1 (%.1f%%)", var_exp[1]),
       y = sprintf("PCo2 (%.1f%%)", var_exp[2]),
       color = COLOR_BY) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))+
 #zoom on axis
   #coord_cartesian(xlim = c(-0.105, 0.01), ylim = c(-0.275, 0.105))        
p

ggsave("20260611_HGT_PCoA_class_top10_zoom.png", device = "png", width = 40, height = 30, units = "cm")
ggsave("20260611_HGT_PCoA_class_top10_zoom.pdf", device = "pdf", width = 40, height = 30, units = "cm")





##PERMANOVA
FACTOR <- "Genus"

meta_df <- as.data.frame(meta)
grp_df <- data.frame(Genome = labels(jac)) %>%
  left_join(meta_df[, c("Genome", FACTOR)], by = "Genome")
keep <- !is.na(grp_df[[FACTOR]]) & grp_df[[FACTOR]] != ""

jac_mat <- as.matrix(jac)[keep, keep]
jac_sub <- as.dist(jac_mat)
grouping <- factor(grp_df[[FACTOR]][keep])

# PERMANOVA: difference in group centroids (composition).
set.seed(42)
perm <- adonis2(jac_sub ~ grouping, permutations = 999)
print(perm)
capture.output(print(perm), file = "20260611_HGT_permanova_Genus.txt")

# PERMDISP / betadisper: difference in group dispersion (the companion check).
disp <- betadisper(jac_sub, grouping)
print(anova(disp))
capture.output(print(perm), file = "20260611_HGT_anova_Genus.txt")

