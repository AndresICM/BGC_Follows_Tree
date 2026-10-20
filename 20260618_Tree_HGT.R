library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(data.table)
library(phyloseq)
library(dplyr)
library(ggnewscale)
library(tidyr)
#library(phytools) 
#library(ggrepel)

setwd("/mnt/c/Users/andre/OneDrive - UW-Madison/Projects/HGT_project")

#Load tree file
treefile <- "20260610_Allactinoscluster_tree_95_2.newick"
tree <- read.tree(treefile)

#Load taxonomy
tax <- fread("20260610_Taxonomy_HighQActinos.csv", data.table=TRUE)
tax <- tax %>% rename(label = Clade)

#load cluster info
cluster_info <- fread("20260616_HGT_all_clades_summary_95_withGCF.csv", data.table=TRUE)


#Get GCF info into long format
cluster_info_long <- cluster_info %>%
  mutate(
    gcf_non_exc = N_GCFs - N_GCFs_clade_exclusive
  ) %>%
  select(Clade, Strain_Count, gcf_non_exc, N_GCFs_clade_exclusive) %>%
  pivot_longer(
    cols = c(gcf_non_exc, N_GCFs_clade_exclusive),
    names_to = "GCF_type",
    values_to = "Count"
  )



# Convert ggtree object to data.frame to access labels
tree_data <- ggtree(tree)
tree_df <- tree_data$data

# Merge
merged_data <- merge(tax, tree_df, by = "label", all.x = TRUE)

#Check unique Phyla
unique_order <- unique(merged_data$Order)


#Define colors
#class really
phylum_colors <- c(
  "Actinomycetes"    = "#CC79A7", 
  "Acidimicrobiia"  = "#0072B2",  
  "Rubrobacteria"  = "#009E73", 
  "Nitriliruptoria" = "#D55E00",
  "Coriobacteriia" = "#E69F00", 
  "Thermoleophilia" = "#56B4E9", 
  "Candidatus Geothermincolia"     = "#F0E442")

#Get class abundance
#Color by most abundant genus
class_abundance <- tax[!is.na(Class), .N, by = Class][order(-N)]
#fwrite(class_abundance, "20260622_HGT_classabund.csv")

#order_abundance <- tax[!is.na(Order), .N, by = Order][order(-N)]
#fwrite(order_abundance, "20260622_HGT_orderabund.csv")

#family_abundance <- tax[!is.na(Family), .N, by = Family][order(-N)]
#fwrite(family_abundance, "20260622_HGT_familyabund.csv")

#Color by most abundant genus
genus_abundance <- tax[!is.na(Genus), .N, by = Genus][order(-N)]
fwrite(genus_abundance, "20260618_HGT_genusabund.csv")
print(genus_abundance)

#Define colors

Genus_colors <- c(
  "Streptomyces"    = "#4277b6", 
  "Mycobacterium"  = "#ebce2b",  
  "Rhodococcus"     = "#df8461",
  "Microbacterium"  = "#ba1c30", 
  "Micromonospora"    = "#5fa641",
  "Nocardia"    = "#d485b2",
  "Cutibacterium"     = "#463397",
  "Kitasatospora"    = "#91218c",
  "Arthrobacter"     = "#7e1510",
  "Micrococcus"     = "#2b3514",
  "Other"= "#7f7e80")

# Create a version with only the top genus from each cluster
cluster_info$Top_Genus <- sapply(cluster_info$Genus_List, function(x) {
  genera <- unlist(strsplit(x, ",\\s*"))
  top <- genera[1]
  if (top %in% names(Genus_colors)) top else "Other"  #Uses my list in genus_colors
})


p <- ggtree(tree, layout = "fan", open.angle=12,  branch.length = "none") %<+% cluster_info

#Add heatmap of Strain counts
p <- p +  geom_fruit(
  geom = geom_tile,
  mapping = aes(
    y = Cluster,
    x = 1,
    fill = Strain_Count_app),
  width = 2.5,
  offset=-0.15) +
  scale_fill_viridis_c(name = "Strain Count", 
                       option = "D",
                       trans = "log10") +
  #geom_tiplab(linesize = 0.001, hjust = -0.001, size = 2) +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 15)  ) +
  new_scale_fill()

#Add tiles with Genus
p <- p +  geom_fruit(
  geom = geom_tile,
  mapping = aes(
    y = Cluster,
    x = 1,
    fill = Top_Genus),
  width = 2.5,
  offset=-0.15) +
  scale_fill_manual(name = "Top 10 Genera", values = Genus_colors) +
  new_scale_fill()

# Add bar plot with GCFs
p <- p +    geom_fruit(
  data = cluster_info_long,
  geom = geom_bar,
  mapping = aes(
    y = Clade,
    x = N_GCFs,
    fill = GCF_type  
  ),
  orientation = "y",
  stat = "identity",
  width = 2,
  offset = 0.03
) +
  scale_fill_manual(
    values = c(
      gcf_non_exc = "#8a8a8a",     
      N_GCFs_clade_exclusive = "#CC79A7"    
    ),
    name = "GCF Type",
    labels = c("Total GCFs", "Clade exclusive GCFs")  
  )

p <- rotate_tree(p, -180)
p

ggsave("20260618_HGT_tree_v1.png", width = 48, height = 48, dpi = 300, plot =last_plot())
ggsave("20260618_HGT_tree_v1.pdf", width = 15, height = 15, dpi = 300, plot =last_plot())


#Highlight some clades
highlighted_clades <- c("clade_4731", "clade_0141", "clade_3448", "clade_6909",
                        "clade_5031", "clade_3755", "clade_4282", "clade_4684",
                        "clade_7385", "clade_2347", "clade_1145", "clade_5692")



p <- p +
  geom_point2(aes(subset = label %in% highlighted_clades),
              size = 1, shape = 21, fill = "red", stroke = 0.8, color = "black") +
  geom_text2(aes(subset = label %in% highlighted_clades, label = label),
             size = 5, hjust = -0.2)


p

ggsave("20260618_HGT_tree_v2.png", width = 35, height = 35, dpi = 300, plot =last_plot())
ggsave("20260618_HGT_tree_v2.pdf", width = 15, height = 15, dpi = 300, plot =last_plot())




#Backup
#Class tiles, don't run
p <- p +  geom_fruit(
  geom = geom_tile,
  mapping = aes(
    y = Cluster,
    x = 1,
    fill = Class_List),
  width = 2.5,
  offset=-0.15) +
  scale_fill_manual(name = "Class", values = phylum_colors) +
  new_scale_fill()
