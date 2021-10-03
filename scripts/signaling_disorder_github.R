#covers signaling/disorder analyses
#contains model H, figure s5, s6
#uses files: 
#supplemental_dataset_1.txt, gene_data.csv, kegg_key.txt, concatenated_root.tre
#iupred_by_gene.csv, iupred_by_stall.csv, supplemental_dataset_2.txt
library(ggalt)
library(caper)
library(car)
library(ape)
library(reshape2)
library(phytools)
library(ggplot2)


######files########
setwd('wherever/you/saved')
mapping <- read.delim('supplemental_dataset_1.txt')
enc <- read.csv('gene_data.csv')
kegg_key <- read.delim('kegg_key.txt', quote='', fill=FALSE) 
kegg_key$c_level[grep('serine protein kinase', kegg_key$description)] <- 'Protein kinases'
#unify description

####create phylogenetic tree####
tree <- read.tree('concatenated_root.tre')
mapping <- mapping[(mapping$taxon_oid %in% tree$tip.label), ]
table(tree$tip.label %in% mapping$taxon_oid) 
#some tree tips not in mapping file but that's ok they fell out at completeness screen
rooted_tree <- root(tree, 'GCF_000009185.1_ASM918v1_genomic', resolve.root=TRUE)
#Haloquadratum walsbyi as root
is.rooted(rooted_tree)

table(duplicated(rooted_tree$node.label)) #duplicated node.labels trip up later programs
rooted_tree$node.label <- make.unique(rooted_tree$node.label)
#edges with length 0 will mess up pgls, just change to small number
table(rooted_tree$edge.length == 0)
rooted_tree$edge.length <- ifelse(rooted_tree$edge.length == 0, 0.00001, rooted_tree$edge.length)
row.names(mapping) <- mapping$taxon_oid
rm(tree)

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

######motifs in signaling proteins########
signaling_keggs <- subset(kegg_key, b_level == 'Protein kinases'| b_level == 'Signal transduction')
signaling <- subset(enc, ko_id %in% signaling_keggs$ko_number)
signaling$c_level <- kegg_key$c_level[match(signaling$ko_id, kegg_key$ko_number)]
sort(table(signaling$c_level))

signal_stalls <- dcast(signaling, genome ~ c_level, value.var = 'total_stalls', fun.aggregate = sum)
signal_stalls$all_stalls <- rowSums(signal_stalls[, 2:ncol(signal_stalls)])

mapping$signaling_stalls <- signal_stalls$all_stalls[match(mapping$taxon_oid, signal_stalls$genome)]
phylogenetic_anova(mapping, 'multicellular', 'signaling_stalls', rooted_tree)


######figure s5#######
#ser-threonine kinases countz multicellular vs. unicellular
ser_thr_kinases <- as.data.frame(table(subset(signaling, ko_id == 'K08884')$genome))
mapping$ser_thr_kinases <- ifelse(mapping$taxon_oid %in% ser_thr_kinases$Var1, 
                                  ser_thr_kinases$Freq[match(mapping$taxon_oid, ser_thr_kinases$Var1)], 0)

figure_s5 <- ggplot(mapping, aes(y=ser_thr_kinases, x= multicellular)) +
  geom_jitter(aes(fill = multicellular), pch = 21, size = 5, alpha = 0.75) +
  geom_boxplot(aes(fill = multicellular), alpha = 0.75, outlier.shape = NA) +
  xlab('') +
  ylab('Serine-threonine kinases') +  
  scale_fill_manual(values = c('#66C2A5', '#E6F598')) +
  theme_linedraw() +
  theme(legend.position = 'NONE') +
  theme(axis.text = element_text(family='Avenir', size=20, color = '#666666')) +
  theme(text = element_text(family='Avenir', size=20, color = '#666666'))

phylogenetic_anova(mapping, 'multicellular', 'ser_thr_kinases', rooted_tree)
ggsave('figure_s5.svg', plot=figure_s5, width = 4.5, height = 4.5)

####PGLS tests -- no motifs from signaling proteins model H####
mapping$nosignaling_stalls <- mapping$total_motifs - mapping$signaling_stalls
pretty_data <- comparative.data(rooted_tree, mapping, taxon_oid, vcv = TRUE, na.omit = FALSE)

#model H
no_signaling_PGLS <- pgls(log(nosignaling_stalls, base = 10) ~ log(total_gc_coding, base = 10) + log(doubling_estimated, base = 10) 
                          + log(trna_count, base = 10) + log(rRNA_count_estimated, base = 10) + temperature_class + multicellular, pretty_data)
summary(no_signaling_PGLS)
coef(no_signaling_PGLS)
rm(signal_genes, signal_stalls, signaling, signaling_keggs)
rm(no_signaling_PGLS, pretty_data, ser_thr_kinases)

#un-annotated stalls
mean(mapping$kegg_unannotated)
aggregate(kegg_unannotated ~ multicellular, mapping, mean)
phylogenetic_anova(mapping, 'multicellular', 'kegg_unannotated', rooted_tree)


######disorder by motifs########
common <- read.delim('supplemental_dataset_2.txt')
common_signaling <- subset(common, b_level == 'Protein kinases' & genomes_w_motif >= 50 |
                             b_level == 'Signal transduction' & genomes_w_motif >= 50)

#look at these signaling proteins, what proportion of residues are disordered
#and do the di-prolyl motifs occur in disordered regions

iupred_per_gene <- read.csv('signaling_iupred_per_gene.csv')
#disordered residues with each protein w/ ko id in common signaling 
iupred_per_motif <- read.csv('signaling_iupred_per_motif.csv')
#disorder of each di-prolyl motif within each protein w/ ko id in common signaling 


chisq_results <- data.frame(matrix(nrow = 0, ncol = 5))
colnames(chisq_results) <- c('ko_id', 'disordered_stalls', 'ordered_stalls', 'disordered_residues', 'p_value')
count = 0
for (kegg in names(table(iupred_per_gene$ko_id))){
  print(kegg)
  count = count + 1
  disordered_stalls = nrow(subset(iupred_per_motif, ko_id == kegg & disordered == 'disorder'))
  ordered_stalls = nrow(subset(iupred_per_motif, ko_id == kegg & disordered == 'order'))
  disordered_residues_percent = sum(subset(iupred_per_gene, ko_id == kegg)$disordered_residues) / sum(subset(iupred_per_gene, ko_id == kegg)$protein_length)
  ch_test <- chisq.test(c(disordered_stalls, ordered_stalls), 
                        p = c(disordered_residues_percent, (1-disordered_residues_percent)))
  chisq_results[count, 'ko_id'] <- kegg
  chisq_results[count, 'disordered_stalls'] <- disordered_stalls
  chisq_results[count, 'ordered_stalls'] <- ordered_stalls
  chisq_results[count, 'disordered_residues'] <- disordered_residues_percent
  chisq_results[count, 'p_value'] <- ch_test$p.value
  rm(disordered_stalls, ordered_stalls, disordered_residues_percent, ch_test)
}

chisq_results$adj_pvalue <- p.adjust(chisq_results$p_value, method = 'bonferroni', n = length(chisq_results$p_value))
chisq_results$percent_disordered_stalls <- chisq_results$disordered_stalls / (chisq_results$disordered_stalls + chisq_results$ordered_stalls)
chisq_results$c_level <- kegg_key$c_level[match(chisq_results$ko_id, kegg_key$ko_number)]
chisq_results$description <- kegg_key$description[match(chisq_results$ko_id, kegg_key$ko_number)]
chisq_results <- chisq_results[order(chisq_results$adj_pvalue), ]
chisq_results$diff <- chisq_results$percent_disordered_stalls - chisq_results$disordered_residues
plot <- subset(chisq_results, diff >= 0.3 & adj_pvalue < 0.0001)

plot$description_short <- c('Two-component system\nOmpR family, sensor kinase', 'Two-component system\nchemotaxis family\nsensor kinase CheA',
                            'Serine protease Do/DegP/HtrA', 'Bacterial Ser-Thr protein kinase', 'Two-component system\nresponse regulator RegA')

figure_s6 <- ggplot(plot, aes(x=disordered_residues, xend = percent_disordered_stalls, y= reorder(description_short, diff), group = description_short)) +
  geom_dumbbell(color = '#7E2872', colour_x = '#9E0142', colour_xend = '#5E4FA2', size = 2.5, 
                size_x = 7, size_xend = 7, show.legend = TRUE) +
  geom_text(data=subset(plot, description_short == 'Two-component system\nresponse regulator RegA'),
            aes(x=disordered_residues, y=description_short, label='Disordered residues', hjust=0),
            color='#9E0142', size=5, vjust=-1.3, fontface='bold', family='Avenir') +
  geom_text(data=plot, aes(x=disordered_residues, y=description_short, label = paste(round(disordered_residues * 100), '%', sep = '')),
            color='#9E0142', size=5, vjust=2.5, family='Avenir') +
  geom_text(data=subset(plot, description_short == 'Two-component system\nresponse regulator RegA'),
            aes(x=percent_disordered_stalls, y=description_short, label='Di-proyl motifs\nin disordered region'),
            color='#5E4FA2', size=5, vjust=-0.3, hjust=1, fontface='bold', family='Avenir') +
  geom_text(data=plot, aes(x=percent_disordered_stalls, y=description_short, label = paste(round(percent_disordered_stalls * 100), '%', sep = '')),
            color='#5E4FA2', size=5, vjust=2.5, family='Avenir') +
  xlab('Proportion') +
  ylab('') +  
  theme_linedraw() +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text = element_text(family='Avenir', size=20, color = '#666666')) +
  theme(text = element_text(family='Avenir', size=20, color = '#666666'))

ggsave('figure_s6.svg', plot=figure_s6, width = 10, height = 8)

#does the enrichment of di-prolyl motifs in IDRs of ser-thr kinases differ between multicellular and unicellular genomes?
iupred_per_motif$multicellular <- mapping$multicellular[match(iupred_per_motif$genome, mapping$taxon_oid)]
table(subset(iupred_per_motif, ko_id == 'K08884')$multicellular, subset(iupred_per_motif, ko_id == 'K08884')$disordered)


