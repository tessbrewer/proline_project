#contains code for simple-ish figures
#figures 1, 2, 5, s1
#uses files: 
#supplemental_dataset_1.txt, concatenated_root.tre, figure5_data.txt

library(ggplot2)
library(RColorBrewer)
library(scales)
library(reshape2)
library(caper)
library(ape)
library(phylolm)
library(phytools)

######figure 1#######
#phyla distribution
setwd('wherever/you/saved')
mapping <- read.delim('supplemental_dataset_1.txt')
mapping$phylum_plot <- ifelse(grepl('Firmicutes', mapping$phylum_gtdbtk), 'Firmicutes', 
                              ifelse(mapping$phylum_gtdbtk == 'Proteobacteria', mapping$class_gtdbtk, mapping$phylum_gtdbtk))
countz <- table(mapping$phylum_plot)
mapping$phylum_plot <- ifelse(mapping$phylum_plot %in% names(countz)[countz >= 5], mapping$phylum_plot, 'Other')

colors <- brewer.pal(11, 'Spectral')
show_col(colors)
plot <- subset(mapping, phylum_plot != 'Other')
figure_1 <- ggplot(plot, aes(y=total_motifs, x=reorder(phylum_plot, total_motifs, median))) +
  geom_jitter(aes(fill = (total_gc_coding/1000000)), size = 5, pch = 21, alpha = 0.75) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA) +
  xlab('') +    
  ylab('Total di-prolyl motifs') +    
  scale_fill_gradientn(colours = rev(colors), name = 'Coding GC content (Mbp)') +
  theme_linedraw() +
  theme(legend.justification = c(0.05, 0.95)) +
  theme(legend.position = c(0.05, 0.95)) +
  theme(text = element_text(family='Avenir', size=20, color = '#666666')) +
  theme(title = element_text(family='Avenir', size=20, color = '#666666')) +
  theme(axis.text.y = element_text(size=20, color = '#666666')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=20, color = '#666666'))

aggregate((total_gc_coding/1000000) ~ phylum_plot, mapping, mean)
# ggsave('figure_1.svg', plot=figure_1, width = 12, height = 8)


####create phylogenetic tree####
tree <- read.tree('concatenated_root.tre')
mapping <- mapping[(mapping$taxon_oid %in% tree$tip.label), ]
table(tree$tip.label %in% mapping$taxon_oid) 
#some tree tips not in mapping file but that's ok they fell out at the genome completeness/contamination screen
rooted_tree <- root(tree, 'GCF_000009185.1_ASM918v1_genomic', resolve.root=TRUE)
is.rooted(rooted_tree)
#root tree
table(duplicated(rooted_tree$node.label)) #doesn't like duplicated node.labels
rooted_tree$node.label <- make.unique(rooted_tree$node.label)
#edges with length 0 will mess up test, just change to small number
table(rooted_tree$edge.length == 0)
rooted_tree$edge.length <- ifelse(rooted_tree$edge.length == 0, 0.00001, rooted_tree$edge.length)
row.names(mapping) <- mapping$taxon_oid


####phylogenetic anova function####
phylogenetic_anova <- function(df, groupz, valuez, tree){
  df <- df[(is.na(df[[groupz]]) == FALSE & is.na(df[[valuez]]) == FALSE), ]
  groups_to_test <- as.vector(df[[groupz]])
  names(groups_to_test) <- df$taxon_oid
  test_values <- as.vector(df[[valuez]])
  names(test_values) <- df$taxon_oid
  test_tree <- drop.tip(phy = tree, tip = tree$tip.label[!tree$tip.label %in% names(groups_to_test)])
  print(paste('Tips in tree:', length(test_tree$tip.label)))
  anova_results <- phylANOVA(test_tree, groups_to_test, test_values, p.adj = 'fdr')
  return(anova_results)
}


######figure 2#######
mapping$motifs_per_gccoding <- (mapping$total_motifs / (mapping$total_gc_coding/1000000))

#thermophiles vs. mesophiles
phylogenetic_anova(subset(mapping, temperature_class != 'psychrophilic'), 'temperature_class', 'motifs_per_gccoding', rooted_tree)
table(mapping$temperature_class)
aggregate(motifs_per_gccoding ~ temperature_class, subset(mapping, temperature_class != 'psychrophilic'), mean)

#multicellular vs. unicellular
phylogenetic_anova(mapping, 'multicellular', 'motifs_per_gccoding', rooted_tree)
table(mapping$multicellular)
aggregate(motifs_per_gccoding ~ multicellular, mapping, mean)

plot <- melt(mapping, id.vars = c('motifs_per_gccoding'), measure.vars = c('multicellular', 'temperature_class'))
plot <- subset(plot, value != 'psychrophilic')
plot$value <- factor(plot$value, levels = c('mesophilic', 'thermophilic', 'unicellular', 'multicellular'))
plot$final_variable <- ifelse(plot$variable == 'multicellular', 'Complex lifecycle', 'Growth temperature')
plot$final_variable <- factor(plot$final_variable, levels = c('Growth temperature', 'Complex lifecycle'))

figure_2 <- ggplot(plot, aes(y=motifs_per_gccoding, x= value)) +
  facet_wrap(~final_variable, scales = 'free_x') + 
  geom_jitter(aes(fill = value), pch = 21, size = 5, alpha = 0.75) +
  geom_boxplot(aes(fill = value), alpha = 0.75, outlier.shape = NA) +
  xlab('') +
  ylab('Di-prolyl motifs\nper GC coding Mbp') +  
  scale_fill_manual(values = c('#ABDDA4', '#F46D43','#E6F598', '#66C2A5'),
                    labels = c('Mesophilic', 'Thermophilic', 'Unicellular', 'Multicellular')) +
  theme_linedraw() +
  theme(strip.background =element_rect(fill='white'))+
  theme(legend.position = 'NONE') +
  theme(strip.text = element_text(family='Avenir', size=24, color = '#666666')) +
  theme(axis.text = element_text(family='Avenir', size=24, color = '#666666')) +
  theme(text = element_text(family='Avenir', size=22, color = '#666666'))

# ggsave('figure_2.svg', plot=figure_2, width = 12, height = 6)


######figure 5#######
#example serine threonine kinase disorder vs. motifs
ser_thr_example <- read.delim('figure5_data.txt')
ser_thr_example$stall <- factor(ser_thr_example$stall, levels = c('S', 'M', 'W', ''))

figure_5 <- ggplot(ser_thr_example, aes(y=IUPRED2, x= POS)) +
  geom_line(color = '#666666') +
  geom_point(data = subset(ser_thr_example, stall != ''), 
             aes(x = POS, y = IUPRED2, fill = stall), pch = 21, size = 5) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = '#5E4FA2', size = 1.1) +
  geom_rect(xmin=21, xmax=291, ymin=-0.005, ymax=0.005, fill = '#66C2A5', color = '#66C2A5') + #limits of pfam
  xlab('Position in protein') +
  ylab('IUPred disorder score') +  
  scale_fill_manual(values = c('#D53E4F', '#FEE08B', '#3288BD'), name = 'Stall strength') +
  theme_linedraw() +
  theme(legend.box.background = element_rect(colour = '#666666', size = 0.75)) +
  theme(legend.justification = c(0.05, 0.95)) +
  theme(legend.position = c(0.05, 0.95)) +
  theme(axis.text = element_text(family='Avenir', size=20, color = '#666666')) +
  theme(text = element_text(family='Avenir', size=20, color = '#666666')) + 
  theme(legend.text = element_text(family='Avenir', size=15, color = '#666666'))

# ggsave('figure_5.svg', plot=figure_5, width = 8, height = 5)


######figure s1#######
#overall correlation between measured and estimated doubling 
mapping <- mapping[order(mapping$doubling_estimated),]
unique_mapping <- mapping[!duplicated(mapping$species), ]
#because I collected measured doubling times using species name, 
#need to remove duplicate species to prevent duplicated data
unique_mapping <- subset(unique_mapping, is.na(doubling_measured)  == FALSE)
nrow(unique_mapping)
cor.test(unique_mapping$doubling_measured, unique_mapping$doubling_estimated)

#check correlation in unique thermophilic genomes
cor.test(subset(unique_mapping, temperature_class == 'thermophilic')$doubling_measured, 
         subset(unique_mapping, temperature_class == 'thermophilic')$doubling_estimated)

#check correlation in unique mesophilic genomes
nrow(subset(unique_mapping, temperature_class == 'mesophilic'))
cor.test(subset(unique_mapping, temperature_class == 'mesophilic')$doubling_measured, 
         subset(unique_mapping, temperature_class == 'mesophilic')$doubling_estimated)


figure_s1 <- ggplot(unique_mapping, aes(x=doubling_measured, y= doubling_estimated)) +
  geom_smooth(data = subset(unique_mapping, temperature_class != 'psychrophilic'), aes(color = temperature_class), method = 'lm', se = FALSE) +
  geom_point(aes(fill = temperature_class), pch = 21, size = 5, alpha = 0.75) +
  ylab('Predicted doubling times (hrs)') +
  xlab('Experimentally measured\ndoubling times (hrs)') +
  scale_y_log10() +
  scale_x_log10() +
  scale_color_manual(values = c('#ABDDA4', '#F46D43')) +
  scale_fill_manual(values = c('#ABDDA4', '#E6F598', '#F46D43')) +
  theme_linedraw() +
  theme(legend.justification = c(0.98, 0.02)) +
  theme(legend.position = c(0.98, 0.02)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(family='Avenir', size=16, color = '#666666')) +
  theme(axis.text = element_text(family='Avenir', size=20, color = '#666666')) +
  theme(text = element_text(family='Avenir', size=20, color = '#666666'))

# ggsave('figure_s1.svg', plot=figure_s1, width = 7, height = 6.75)
