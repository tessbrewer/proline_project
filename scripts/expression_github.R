#covers expected expression rank section
#contains model G, figure 4, s3, s4
#uses files: 
#supplemental_dataset_1.txt, gene_data.csv, kegg_key.txt, concatenated_root.tre
library(phylolm)
library(phytools)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(reshape2)
library(caper)
library(car)
library(ape)

######functions########
scientific_10 <- function(x) {
  #classes up scientific notation
  parse(text=gsub('e', ' %*% 10^', scales::scientific_format()(x)))
}

######files########
setwd('wherever/you/saved')
kegg_key <- read.delim('kegg_key.txt', quote='', fill=FALSE) 
mapping <- read.delim('supplemental_dataset_1.txt')

#######phylogenetic tree#######
tree <- read.tree('concatenated_root.tre')
mapping <- mapping[(mapping$taxon_oid %in% tree$tip.label), ]
table(tree$tip.label %in% mapping$taxon_oid) 
#some tree tips not in mapping file but that's ok they fell out at completeness screen
rooted_tree <- root(tree, 'GCF_000009185.1_ASM918v1_genomic', resolve.root=TRUE)
is.rooted(rooted_tree)
#root tree
table(duplicated(rooted_tree$node.label)) #doesn't like duplicated node.labels
rooted_tree$node.label <- make.unique(rooted_tree$node.label)
#edges with length 0 will mess up test, just change to small number
table(rooted_tree$edge.length == 0)
rooted_tree$edge.length <- ifelse(rooted_tree$edge.length == 0, 0.00001, rooted_tree$edge.length)
row.names(mapping) <- mapping$taxon_oid

#######prepare MEDIAN rank file#######
enc <- read.csv('gene_data.csv')
#aggregate rank of each ko
ko_count <- aggregate(rank ~ ko_id, subset(enc, !duplicated(paste(enc$genome, enc$ko_id))), length)
#count each ko occurrence per genome once
#median rank
avg_rank_median <- aggregate(rank ~ ko_id, enc, median)
avg_rank_median <- subset(avg_rank_median, ko_id != '')
avg_rank_median$genomes_w_ko <- ko_count$rank[match(avg_rank_median$ko_id, ko_count$ko_id)]
rm(ko_count)

#genomes w/ motif / strong motif in each ko
w_stalls <- subset(enc, total_stalls > 0)
stall_count <- aggregate(total_stalls ~ ko_id, subset(w_stalls, !duplicated(paste(w_stalls$genome, w_stalls$ko_id))), length)
w_str_stalls <- subset(w_stalls, S > 0)
str_stall_count <- aggregate(S ~ ko_id, subset(w_str_stalls, !duplicated(paste(w_str_stalls$genome, w_str_stalls$ko_id))), length)

avg_rank_median$genomes_w_stall <- stall_count$total_stalls[match(avg_rank_median$ko_id, stall_count$ko_id)]
avg_rank_median$genomes_w_str_stall <- str_stall_count$S[match(avg_rank_median$ko_id, str_stall_count$ko_id)]
rm(w_stalls, w_str_stalls, stall_count, str_stall_count)
avg_rank_median$genomes_w_stall <- (avg_rank_median$genomes_w_stall / avg_rank_median$genomes_w_ko) * 100
avg_rank_median$genomes_w_str_stall <- (avg_rank_median$genomes_w_str_stall / avg_rank_median$genomes_w_ko) * 100
avg_rank_median$percent_genomes_w_ko <- (avg_rank_median$genomes_w_ko / length(table(enc$genome))) * 100

#average motifs gc controlled
avg_stalls <- aggregate(motifs_per_100aa_gc ~ ko_id, enc, mean)
avg_rank_median$avg_stalls_100aa_gc <- avg_stalls$motifs_per_100aa_gc[match(avg_rank_median$ko_id, avg_stalls$ko_id)]
#average strong motifs gc controlled
avg_stalls <- aggregate(str_motifs_per_100aa_gc ~ ko_id, enc, mean)
avg_rank_median$avg_str_stalls_100aa_gc <- avg_stalls$str_motifs_per_100aa_gc[match(avg_rank_median$ko_id, avg_stalls$ko_id)]
rm(avg_stalls)

#add annotations
avg_rank_median$a_level <- kegg_key$a_level[match(avg_rank_median$ko_id, kegg_key$ko_number)]
avg_rank_median$a_level[is.na(avg_rank_median$a_level)] <- 'Unclassified'
avg_rank_median$b_level <- kegg_key$b_level[match(avg_rank_median$ko_id, kegg_key$ko_number)]
avg_rank_median$c_level <- kegg_key$c_level[match(avg_rank_median$ko_id, kegg_key$ko_number)]
avg_rank_median$description <- kegg_key$description[match(avg_rank_median$ko_id, kegg_key$ko_number)]

common <- subset(avg_rank_median, genomes_w_ko > length(table(mapping$taxon_oid))*0.25) # has to be in 25% of genomes
common$a_level[common$a_level == 'Human Diseases'] <- 'Environmental Information Processing' #weird classification
# write.table(common, 'supplemental_dataset_2.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#######rank expression vs motifs#######
#top 5 highly expressed by median
head(common[order(common$rank), ], n = 5)[,c(2, 4, 12)] 
#rank vs. avg motifs
cor.test(common$rank, common$avg_stalls_100aa_gc, method = 'spearman', exact = FALSE)
#rank vs. avg strong motifs
cor.test(common$rank, common$avg_str_stalls_100aa_gc, method = 'spearman', exact = FALSE)

#how do percentiles compare
round(mean(subset(common, rank <= quantile(rank, 0.01)[[1]])$avg_stalls_100aa_gc), 3) #top 1%
round(mean(subset(common, rank >= quantile(rank, 0.99)[[1]])$avg_stalls_100aa_gc), 3) #bottom 1%

#figure 4
colors <- brewer.pal(length(table(common$a_level)), 'Spectral')
categories <- names(table(common$a_level))
show_col(colors)
plot <- melt(common, id.vars = c('ko_id', 'rank', 'a_level'), measure.vars = c('avg_stalls_100aa_gc', 'avg_str_stalls_100aa_gc'))
plot$a_level[is.na(plot$a_level)] <- 'Unclassified'
plot$a_level <- gsub('Environmental Information Processing', 'Environmental\nInformation Processing', plot$a_level)
plot$a_level <- gsub('Genetic Information Processing', 'Genetic\nInformation Processing', plot$a_level)

#just the total motifs
figure4 <- ggplot(subset(plot, variable == 'avg_stalls_100aa_gc'), aes(y = rank, x = value)) +
  geom_jitter(aes(fill = a_level), pch = 21, size = 5, alpha = 0.75) +
  geom_smooth(aes(group = variable), method = 'lm', color = 'purple', fill = 'dark grey') +
  ylab('Median codon bias rank') +
  xlab('Average di-prolyl motifs per 100AA\n(GC-controlled)') +
  scale_x_log10(label=scientific_10) +
  scale_fill_brewer(palette = 'Spectral') +
  theme_linedraw() +
  theme(text = element_text(family='Avenir', size=22, color = '#666666')) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(family='Avenir', size=18, color = '#666666')) +
  theme(axis.text = element_text(family='Avenir', size=22, color = '#666666'))

subset(plot, variable == 'avg_stalls_100aa_gc' & value == 0) 
#lost one ko with 0 motifs on average b/c of log transformed axis
# ggsave('figure_4.svg', plot=figure4, width = 10.25, height = 7)

#just the strong motifs
figureS3 <- ggplot(subset(plot, variable == 'avg_str_stalls_100aa_gc'), aes(y = rank, x = value)) +
  geom_jitter(aes(fill = a_level), pch = 21, size = 5, alpha = 0.75) +
  geom_smooth(aes(group = variable), method = 'lm', color = 'purple', fill = 'dark grey') +
  ylab('Median codon bias rank') +
  xlab('Average strong di-prolyl motifs per 100AA\n(GC-controlled)') +
  scale_x_log10(label=scientific_10) +
  scale_fill_brewer(palette = 'Spectral') +
  theme_linedraw() +
  theme(text = element_text(family='Avenir', size=22, color = '#666666')) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(family='Avenir', size=18, color = '#666666')) +
  theme(axis.text = element_text(family='Avenir', size=22, color = '#666666'))

nrow(subset(plot, variable == 'avg_str_stalls_100aa_gc' & value == 0))
#lost nine kos with 0 motifs on average b/c of log transformed axis
# ggsave('figure_s3.svg', plot=figureS3, width = 10.25, height = 7)


#######genes in top percentile for CUB per genome#######
percentile_1 <- aggregate(Ncp ~ genome, enc, function(x){quantile(x, 0.01)}) 
#find the top percentile for Ncp for each genome
percentile_1$stalls_in_quantile <- NA
percentile_1$str_stalls_in_quantile <- NA
percentile_1$no_genes_in_quantile <- NA

top_0.5 <- enc[FALSE, ]
for(i in 1:nrow(percentile_1)){
  #extract genes under the top percentile cutoff
  temp <- subset(enc, genome == percentile_1$genome[i] & Ncp <= percentile_1$Ncp[i])
  percentile_stalls = sum(temp$total_stalls)
  str_percentile_stalls = sum(temp$S)
  no_genes =  nrow(temp)
  percentile_1$stalls_in_quantile[i] <- percentile_stalls
  percentile_1$str_stalls_in_quantile[i] <- str_percentile_stalls
  percentile_1$no_genes_in_quantile[i] <- no_genes
  top_0.5 <- rbind(top_0.5, temp)
  rm(percentile_stalls, temp, str_percentile_stalls, no_genes)
}

top_0.5$species <- mapping$species[match(top_0.5$genome, mapping$taxon_oid)]
top_0.5$b_level <- kegg_key$b_level[match(top_0.5$ko_id, kegg_key$ko_number)]
top_0.5$c_level <- kegg_key$c_level[match(top_0.5$ko_id, kegg_key$ko_number)]

#figure S4
countz <- table(top_0.5$b_level)
top_0.5$b_level_plot <- ifelse(top_0.5$b_level %in% names(head(rev(sort(countz)), n = 9)), top_0.5$b_level, 'Other')
#only give color to top 8 most abundant categories
plot <- dcast(species ~ b_level_plot, data = top_0.5, value.var = 'species', fun.aggregate = length)
plot[,2:ncol(plot)] <- plot[,2:ncol(plot)] / rowSums(plot[,2:ncol(plot)]) #proportional
plot <- melt(plot, id.vars = 'species')
plot$doubling <- mapping$doubling_estimated[match(plot$species, mapping$species)]
plot <- subset(plot, is.na(species) == FALSE & value != 0)
plot$growth <- ifelse(plot$doubling > 5, 'slow', 'fast')

#mean percent of top CUB in translation slow vs fast
aggregate(value ~ growth, subset(plot, variable == 'Translation'), mean)

figureS4 <- ggplot(plot, aes(y = value, x = reorder(species,doubling), fill = variable, color = variable)) +
  geom_bar(position='fill', stat='identity') +
  xlab('') +
  ylab('Proportion of genes\nin top percentile for CUB') +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = 'Spectral') +
  scale_color_brewer(palette = 'Spectral') +
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank()) +
  theme(text = element_text(family='Avenir', size=20, color = '#666666'))

# ggsave('figure_s4.svg', plot=figureS4, width = 12, height = 7)


#average motifs in highly expressed proteins
percentile_1$stalls_per_gene <- percentile_1$stalls_in_quantile / percentile_1$no_genes_in_quantile
percentile_1$str_stalls_per_gene <- percentile_1$str_stalls_in_quantile / percentile_1$no_genes_in_quantile
mapping$stalls_percentile <- percentile_1$stalls_per_gene[match(mapping$taxon_oid, percentile_1$genome)]

#model G
pretty_data <- comparative.data(rooted_tree, mapping, taxon_oid, vcv=TRUE, na.omit = FALSE)
high_exp_motifs <- pgls(stalls_percentile ~ log(total_gc_coding, base = 10) + log(doubling_estimated, base = 10) 
                          + log(trna_count, base = 10) + log(rRNA_count_estimated, base = 10) + temperature_class + multicellular, pretty_data)
summary(high_exp_motifs)
round(coef(high_exp_motifs), 3)

