library(ampvis2)
library(BiocManager)
library(ggplot2)
library(dplyr)
library(phyloseq)
library(ranacapa)
library(readxl)
library(tidyverse)
library(vegan)
library(psadd)
library(ape)

#Import files-------------------------------------------------------------------------------------------------------------
otu_mat = read.delim("D:/Postdoc/1) R_Studio/MX56/MX56_9G/Table_9G.txt",row.names=1)
tax_mat = read.delim("D:/Postdoc/1) R_Studio/MX56/MX56_9G/Clasification_9G.txt",row.names=1)
samples_df = read.delim("D:/Postdoc/1) R_Studio/MX56/MX56_9G/Metadata_9G.txt",row.names=1)

#Transform into matrixes otu and tax tables
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#Transform to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
MX56 <- phyloseq (OTU, TAX, samples)
MX56 <- prune_samples(sample_names(MX56) != "Mock", MX56) # Remove potential synthetic samples
MX56

#Distance analysis (Beta diversity)-------------------------------------------------------------------------------------------------------

d = distance(MX56, method='bray')                   #method bray_curtis
plot(hclust(d, method="ward"), xlab="bray_curtis")
#
d = distance(MX56, method='jaccard')                #method bray_curtis
plot(hclust(d, method="ward"), xlab="jaccard")
#
hc <- hclust(d, method="ward")
my_tree <- as.phylo(hc) 
write.tree(phy=my_tree, file="bray_treeV.newick")   #exporting tree as nwk for edition

#Composition analysis-------------------------------------------------------------------------------------------------------------------------------------------------
#Set colors
phylum_colors <-c("purple","red","#009E73", "#FF2F80", "red4","blue","#CBD588","#FC8E00","cyan","#652CFF","maroon4","#575329", "#00FECF", "#B05B6F",
                          "#8CD0FF","darkgrey", "#F6E442", "#04F757", "#FFAA92", "#FF90C9", "#B903AA", "#C8A1A1", "royalblue4","#CBD588",
                  "#1E6E00", "#320033", "#66E1D3", "#CFCDAC", "#4FC601", "#4A3B53",
                          "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#D16100", "#575329", 
                           "#04F757", "#C8A1A1", "#1E6E00", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                          "#3B5DFF", "#549E79","#7B4F4B", "#A1C299", "purple", "#999999", "#E69F00",  "#009E73","darkgrey","#FC4E07",  
                          "#652926", "red4",  "#009E73","#00BA38" , "#252A52" , "brown",    "#D55E00", "#5F7FC7", 
                          "#CC79A7","orange","#56B4E9","#DA5724", "#508578","#C84248", "darkorchid", "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5",
                          "darkseagreen", "#5E738F","#D1A33D","#8A7C64", "#599861", "yellow","darkgoldenrod1", "#56B4E9","darkolivegreen1","#F0E442", "#0072B2", "#D55E00")
#Barplot order
MX <- MX56%>%
  tax_glom(taxrank = "order") %>%                                 # agglomerate at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>%            # Transform to rel. abundance
  psmelt() %>%                                                    # Melt to long format
  filter(Abundance > 0.05) %>%                                    # Filter out low abundance taxa
  arrange(order)     
#
Z<- ggplot(MX, aes(x = "sample", y = Abundance, fill = order)) +
  geom_bar(stat = "identity",position = "Fill") + 
  scale_fill_manual(values = phylum_colors) +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 5%) \n") +
  ggtitle("Order_abundance_MX56")

Z <- Z+ theme(axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
Z <- Z + facet_grid(~Bray_Group, scales = "free", space='free')
Z
#For excluding Unknown orders
#MX56P <- subset_taxa(MX56, (order!="Unknown_Agaricomycetes") | is.na(order))
#MX56P <- subset_taxa(MX56P, (order!="Unknown_Ascomycota") | is.na(order))
#MX56P <- subset_taxa(MX56P, (order!="Unknown_Fungi_subkgd_incertae_sedis") | is.na(order))
#MX56P <- subset_taxa(MX56P, (order!="Unknown_Pezizomycotina_class_incertae_sedis") | is.na(order))
#MX56P <- subset_taxa(MX56P, (order!="Unknown_Rozellomycota") | is.na(order))
#MX56P

#Barplot primary lifestyle
MX <- MX56%>%
  tax_glom(taxrank = "primary_lifestyle") %>%                      # agglomerate at primary lifestyle level
  transform_sample_counts(function(x) {x/sum(x)} ) %>%             # Transform to rel. abundance
  psmelt() %>%                                                     # Melt to long format
  filter(Abundance > 0.01) %>%                                     # Filter out low abundance taxa
  arrange(primary_lifestyle)     
#
Z<- ggplot(MX, aes(x = Localidad, y = Abundance, fill = order)) +
  geom_bar(stat = "identity",position = "Fill") + 
  scale_fill_manual(values = phylum_colors) +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 1%) \n") +
  ggtitle("primary_lifestyle_abundance_MX56")

Z <- Z+ theme(axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
Z <- Z + facet_grid(~Bray_Group, scales = "free", space='free')
Z
#For merging samples by bray groups (mean)
MX56mean <- merge_samples(MX56, "Bray_Group", fun=mean)
MX56mean                                                
glom_o <- tax_glom(MX56mean, taxrank = "primary_lifestyle")
write.csv(glom_o@tax_table, file = "D:/Postdoc/1) R_Studio/MX56/MX56_9G/9Gmean_primary_lifestyle9G_tax.csv")
write.csv(glom_o@otu_table, file = "D:/Postdoc/1) R_Studio/MX56/MX56_9G/9Gmean_primary_lifestyle9G_table.csv") #repeat barplots with new table, clasification, and metadata

#Alpha diversity---------------------------------------------------------------------------------------

rich = estimate_richness(MX56)
rich

alpha_meas = c("Observed")
(p <- plot_richness(MX56, "sample", color = "vegetacion", measures=alpha_meas))
z <- p + facet_grid(~Bray_Group, scales = "free", space='free')
z
alpha_meas = c("Chao1")
(p <- plot_richness(MX56, "sample", color = "vegetacion", measures=alpha_meas))
z <- p + facet_grid(~Bray_Group, scales = "free", space='free')
z
alpha_meas = c("Shannon")
(p <- plot_richness(MX56, "sample", color = "vegetacion", measures=alpha_meas))
z <- p + facet_grid(~Bray_Group, scales = "free", space='free')
z

#Ampvis (Heatmap Top 50)--------------------------------------------------------------------------------------------------------------
#From phyloseq to ampvis
Totu_table =t(otu_table(MX56))
otu_table(MX56)=Totu_table
av2_otutable <- data.frame(OTU = rownames(t(phyloseq::otu_table(MX56)@.Data)),
                           t(phyloseq::otu_table(MX56)@.Data),
                           phyloseq::tax_table(MX56)@.Data,
                           check.names = F)
#Extract metadata from the phyloseq object:
av2_metadata <- data.frame(phyloseq::sample_data(MX56), 
                           check.names = F)
av2_metadata <- cbind(rownames(av2_metadata), av2_metadata)
#Load the data with amp_load:
av2_obj <- amp_load(av2_otutable, av2_metadata)
amp_rankabundance(av2_obj, group_by = "Localidad")
sqrt
log10

#Heatmap 
amp_heatmap(av2_obj,
            group_by = "sample",
            facet_by = "Bray_Group",
            tax_aggregate = "Species",
            tax_show = 50,
            color_vector = c("white", "blue"),
            plot_colorscale = "sqrt",
            plot_values = FALSE) +
  theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
        axis.text.y = element_text(size=12),
        legend.position="left")

#NMDS----------------------------------------------------------------------------------------------
tab = read.delim("D:/Postdoc/1) R_Studio/MX56/MX56_9G/Table_9G.txt",row.names=1)
env = read.delim("D:/Postdoc/1) R_Studio/MX56/MX56_9G/Metadata_9G.txt",row.names=1)
tab
env
#definir columnas tabla y metadata
com = tab[,1:56]
com
envi = env[,20:24]
envi
#convert com to a matrix
m_com = as.matrix(com)
##Perform the NMDS ordination,nmds code
set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds
#Graph
plot(nmds, type = "t") 
plot(envi)

