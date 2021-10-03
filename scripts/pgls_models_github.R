#covers main models
#contains model A-F, figure 3, s2
#uses files: 
#supplemental_dataset_1.txt, concatenated_root.tre
library(ggplot2)
library(caper)
library(ape)
library(ggpubr)

setwd('wherever/you/saved')
mapping <- read.delim('supplemental_dataset_1.txt')

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


####PGLS tests -- most models####
pretty_data <- comparative.data(rooted_tree, mapping, taxon_oid, vcv = TRUE, na.omit = FALSE)

#model A
temperature_PGLS <- pgls(log(total_motifs, base = 10) ~ log(total_gc_coding, base = 10) + temperature_class, pretty_data)
summary(temperature_PGLS)
round(coef(temperature_PGLS), 3)
rm(temperature_PGLS)

#model B
multicellular_PGLS <- pgls(log(total_motifs, base = 10) ~ log(total_gc_coding, base = 10) + temperature_class + multicellular, pretty_data)
summary(multicellular_PGLS)
round(coef(multicellular_PGLS), 3)
rm(multicellular_PGLS)

#model C
nomyxo_data <- comparative.data(rooted_tree, subset(mapping, phylum_gtdbtk != 'Myxococcota'), taxon_oid, vcv = TRUE, na.omit = FALSE)
multicellular_nomyxo_PGLS <- pgls(log(total_motifs, base = 10) ~ log(total_gc_coding, base = 10) + temperature_class + multicellular, nomyxo_data)
summary(multicellular_nomyxo_PGLS)
round(coef(multicellular_nomyxo_PGLS), 3)
rm(nomyxo_data, multicellular_nomyxo_PGLS)

#model D
pred_growth_rates_PGLS <- pgls(log(total_motifs, base = 10) ~ log(total_gc_coding, base = 10) + log(doubling_estimated, base = 10) 
                         + log(trna_count, base = 10) + log(rRNA_count_estimated, base = 10) + temperature_class + multicellular, pretty_data)
summary(pred_growth_rates_PGLS)
round(coef(pred_growth_rates_PGLS), 3)

#model E
mapping <- mapping[order(mapping$doubling_estimated),]
unique_exp <- mapping[!duplicated(mapping$species), ]
unique_exp <- subset(unique_exp, is.na(doubling_measured) == FALSE)
exp_growth_data <- comparative.data(rooted_tree, unique_exp, taxon_oid, vcv = TRUE, na.omit = FALSE)
exp_growth_rates_PGLS <- pgls(log(total_motifs, base = 10) ~ log(total_gc_coding, base = 10) + log(doubling_measured, base = 10) 
                                    + log(trna_count, base = 10) + log(rRNA_count_estimated, base = 10) + temperature_class + multicellular, exp_growth_data)
summary(exp_growth_rates_PGLS)
round(coef(exp_growth_rates_PGLS), 3)
rm(exp_growth_rates_PGLS, exp_growth_data)

#model F
pretty_data_minus_thermo <- comparative.data(rooted_tree, subset(unique_exp, temperature_class == 'mesophilic'), taxon_oid, vcv=TRUE)
exp_growth_rates_PGLS_nothermo <- pgls(log(total_motifs, base = 10) ~ log(total_gc_coding, base = 10) + log(doubling_measured, base = 10) 
                    + log(trna_count, base = 10) + log(rRNA_count_estimated, base = 10) + multicellular, pretty_data_minus_thermo)
summary(exp_growth_rates_PGLS_nothermo)
round(coef(exp_growth_rates_PGLS_nothermo), 3)
rm(unique_exp)


####Figure 3####
#use the equation from model D (pred_growth_rates_PGLS), plug in values for fast and slow example species
slow_example <- subset(mapping, species == 'Methylomagnum ishizawai')
fast_example <- subset(mapping, species == 'Propionigenium maris')

slowest_growth = log(slow_example$doubling_estimated, base = 10)*coef(pred_growth_rates_PGLS)['log(doubling_estimated, base = 10)']
fastest_growth = log(fast_example$doubling_estimated, base = 10)*coef(pred_growth_rates_PGLS)['log(doubling_estimated, base = 10)']

fast_trna = log(fast_example$trna_count, base = 10)*coef(pred_growth_rates_PGLS)['log(trna_count, base = 10)']
slow_trna = log(slow_example$trna_count, base = 10)*coef(pred_growth_rates_PGLS)['log(trna_count, base = 10)']

fast_rrna = log(fast_example$rRNA_count_estimated, base = 10)*coef(pred_growth_rates_PGLS)['log(rRNA_count_estimated, base = 10)']
slow_rrna = log(slow_example$rRNA_count_estimated, base = 10)*coef(pred_growth_rates_PGLS)['log(rRNA_count_estimated, base = 10)']

#calculations -- both unicellular so add that
m = coef(pred_growth_rates_PGLS)['log(total_gc_coding, base = 10)'][[1]]
b_slow = coef(pred_growth_rates_PGLS)[[1]] + coef(pred_growth_rates_PGLS)[['multicellularunicellular']] + slowest_growth + slow_rrna + slow_trna
b_fast = coef(pred_growth_rates_PGLS)[[1]] + coef(pred_growth_rates_PGLS)[['multicellularunicellular']] + fastest_growth + fast_rrna + fast_trna

input_gc_1 = 8.00e06
fast_motifs_1 = 10^(m*(log(input_gc_1, base = 10)) + b_fast[[1]])
slow_motifs_1 = 10^(m*(log(input_gc_1, base = 10)) + b_slow[[1]])
input_gc_2 = 3.5e05
fast_motifs_2 = 10^(m*(log(input_gc_2, base = 10)) + b_fast[[1]])
slow_motifs_2 = 10^(m*(log(input_gc_2, base = 10)) + b_slow[[1]])
(slow_motifs_1-fast_motifs_1)/fast_motifs_1
#slow growth 14% increase

#full plot
pgls_plot <- ggplot(pretty_data$data, aes(y = total_motifs, x = total_gc_coding)) +
  geom_point(pch = 21, size = 2, alpha = 0.75, color = 'white', fill = '#666666') +
  geom_abline(slope = coef(pred_growth_rates_PGLS)['log(total_gc_coding, base = 10)'],
              intercept = coef(pred_growth_rates_PGLS)[[1]] + coef(pred_growth_rates_PGLS)[['multicellularunicellular']] +
              slowest_growth + slow_trna + slow_rrna,
              color = 'white', size = 1.25) +
  geom_abline(slope = coef(pred_growth_rates_PGLS)['log(total_gc_coding, base = 10)'],
              intercept = coef(pred_growth_rates_PGLS)[[1]] + coef(pred_growth_rates_PGLS)[['multicellularunicellular']] +
              slowest_growth + slow_trna + slow_rrna,
              color = '#3288BD', size = 0.75) +
  #double abline is just to give a white outline to the final line
  geom_abline(slope = coef(pred_growth_rates_PGLS)['log(total_gc_coding, base = 10)'],
              intercept = coef(pred_growth_rates_PGLS)[[1]] + coef(pred_growth_rates_PGLS)[['multicellularunicellular']] +
              fastest_growth + fast_trna + fast_rrna,
              color = 'white', size = 1.25) +
  geom_abline(slope = coef(pred_growth_rates_PGLS)['log(total_gc_coding, base = 10)'],
              intercept = coef(pred_growth_rates_PGLS)[[1]] + coef(pred_growth_rates_PGLS)[['multicellularunicellular']] +
              fastest_growth + fast_trna + fast_rrna,
              color = '#D53E4F', size = 0.75) +
  annotate('rect', xmin=4.5e06, xmax=8.25e06, ymin=5200, ymax= 11000, color = 'black', fill = 'white', alpha = 0) +
  annotate('rect', xmin=3.25e05, xmax=5.85e05, ymin=325, ymax= 700, color = 'black', fill = 'white', alpha = 0) +
  scale_y_log10() +
  scale_x_log10() +
  xlab('') +
  ylab('') + 
  theme_linedraw() +
  theme(legend.position = 'NONE') +
  theme(axis.text = element_text(family='Avenir', size=16, color = '#666666'))


#points to annotate
pointz <- data.frame(x=c(input_gc_1, input_gc_1, input_gc_2, input_gc_2), 
                     y=c(fast_motifs_1, slow_motifs_1, fast_motifs_2, slow_motifs_2), 
                     type=c('inset2', 'inset2', 'inset1', 'inset1'),
                     growth=c('fast', 'slow', 'fast', 'slow'))

#inset 1 -- low gc
inset1 <- ggplot() +
  geom_point(data = pretty_data$data, aes(y = total_motifs, x = total_gc_coding), pch = 21, size = 10, alpha = 0.75, color = 'white', fill = '#666666') +
  geom_abline(slope = coef(pred_growth_rates_PGLS)['log(total_gc_coding, base = 10)'],
              intercept = coef(pred_growth_rates_PGLS)[[1]] + coef(pred_growth_rates_PGLS)[['multicellularunicellular']] +
                slowest_growth + slow_trna + slow_rrna,
              color = 'white', size = 3) +
  geom_abline(slope = coef(pred_growth_rates_PGLS)['log(total_gc_coding, base = 10)'],
              intercept = coef(pred_growth_rates_PGLS)[[1]] + coef(pred_growth_rates_PGLS)[['multicellularunicellular']] +
                fastest_growth + fast_trna + fast_rrna,
              color = 'white', size = 3) +
  geom_abline(slope = coef(pred_growth_rates_PGLS)['log(total_gc_coding, base = 10)'],
              intercept = coef(pred_growth_rates_PGLS)[[1]] + coef(pred_growth_rates_PGLS)[['multicellularunicellular']] +
                slowest_growth + slow_trna + slow_rrna,
              color = '#3288BD', size = 2.5) +
  geom_abline(slope = coef(pred_growth_rates_PGLS)['log(total_gc_coding, base = 10)'],
              intercept = coef(pred_growth_rates_PGLS)[[1]] + coef(pred_growth_rates_PGLS)[['multicellularunicellular']] +
                fastest_growth + fast_trna + fast_rrna,
              color = '#D53E4F', size = 2.5) +
  geom_point(data = subset(pointz, type == 'inset1'), aes(y = y, x = x), color = 'white', fill = '#666666', pch = 23, size = 6, alpha = 1) +
  geom_line(data = subset(pointz, type == 'inset1'), aes(y = y, x = x, group = type), color = '#666666', size = 1.75, linetype = 'dashed') +
  geom_text(data = subset(pointz, type == 'inset1' & growth == 'fast'), aes(x=x-5000, y=y+5, label = round(y, 0)), color = '#D53E4F', family='Avenir', size=6, hjust = 1) +
  geom_text(data = subset(pointz, type == 'inset1' & growth == 'slow'), aes(x=x-5000, y=y+5, label = round(y, 0)), color = '#3288BD', family='Avenir', size=6, hjust = 1) +
  scale_x_log10(limits = c(3.25e05, 5.85e05)) +
  scale_y_log10(limits = c(325, 700)) +
  xlab('') +
  ylab('') + 
  theme_linedraw() +
  theme(legend.position = 'NONE') +
  theme(text = element_text(family='Avenir', size=20, color = '#666666')) +
  theme(axis.text = element_text(family='Avenir', size=20, color = '#666666'))

#inset 2 -- high gc
inset2 <- ggplot() +
  geom_point(data = pretty_data$data, aes(y = total_motifs, x = total_gc_coding), pch = 21, size = 10, alpha = 0.75, color = 'white', fill = '#666666') +
  geom_abline(slope = coef(pred_growth_rates_PGLS)['log(total_gc_coding, base = 10)'],
              intercept = coef(pred_growth_rates_PGLS)[[1]] + coef(pred_growth_rates_PGLS)[['multicellularunicellular']] +
                slowest_growth + slow_trna + slow_rrna,
              color = 'white', size = 3) +
  geom_abline(slope = coef(pred_growth_rates_PGLS)['log(total_gc_coding, base = 10)'],
              intercept = coef(pred_growth_rates_PGLS)[[1]] + coef(pred_growth_rates_PGLS)[['multicellularunicellular']] +
                fastest_growth + fast_trna + fast_rrna,
              color = 'white', size = 3) +
  geom_abline(slope = coef(pred_growth_rates_PGLS)['log(total_gc_coding, base = 10)'],
              intercept = coef(pred_growth_rates_PGLS)[[1]] + coef(pred_growth_rates_PGLS)[['multicellularunicellular']] +
                slowest_growth + slow_trna + slow_rrna,
              color = '#3288BD', size = 2.5) +
  geom_abline(slope = coef(pred_growth_rates_PGLS)['log(total_gc_coding, base = 10)'],
              intercept = coef(pred_growth_rates_PGLS)[[1]] + coef(pred_growth_rates_PGLS)[['multicellularunicellular']] +
                fastest_growth + fast_trna + fast_rrna,
              color = '#D53E4F', size = 2.5) +
  geom_point(data = subset(pointz, type == 'inset2'), aes(y = y, x = x), color = 'white', fill = '#666666', pch = 23, size = 6, alpha = 1) +
  geom_line(data = subset(pointz, type == 'inset2'), aes(y = y, x = x, group = type), color = '#666666', size = 1.75, linetype = 'dashed') +
  geom_text(data = subset(pointz, type == 'inset2' & growth == 'fast'), aes(x=x-100000, y=y+200, label = round(y, 0)), color = '#D53E4F', family='Avenir', size=6, hjust = 1) +
  geom_text(data = subset(pointz, type == 'inset2' & growth == 'slow'), aes(x=x-80000, y=y+200, label = round(y, 0)), color = '#3288BD', family='Avenir', size=6, hjust = 1) +
  scale_x_log10(limits = c(4.5e06, 8.25e06)) +
  scale_y_log10(limits = c(5200, 11000)) +
  xlab('') +
  ylab('') + 
  theme_linedraw() +
  theme(legend.position = 'NONE') +
  theme(text = element_text(family='Avenir', size=20, color = '#666666')) +
  theme(axis.text = element_text(family='Avenir', size=20, color = '#666666'))

#assemble full plot
dummy_plot <- ggplot() + theme_void()
overview <- ggarrange(dummy_plot, pgls_plot, dummy_plot, ncol=1, nrow=3, align = 'v',
                   widths = c(1, 1), heights = c(1, 1))
insets <- ggarrange(inset1, inset2, nrow=2, ncol=1, align = 'hv', 
                   widths = c(0.75, 1), heights = c(1, 1))

final_plot_3 <- ggarrange(overview, insets, ncol=2, align = 'v',
                        widths = c(0.8, 1.25), heights = c(1, 1))

# ggsave('figure_3.svg', plot=final_plot_3, width = 13, height = 10)

#clean up to minimize loss of sanity
rm(dummy_plot, final_plot_3, inset1, inset2, insets, overview, pgls_plot, pointz, pretty_data)
rm(b_fast, fast_motifs_1, fast_motifs_2, fast_rrna, fast_trna, fastest_growth, m)
rm(b_slow, slow_motifs_1, slow_motifs_2, slow_rrna, slow_trna, slowest_growth, input_gc_1, input_gc_2)
rm(pred_growth_rates_PGLS)

####Figure S2####
#same as figure 3, just using experimental growth rates of mesophiles only
#use the equation from model F (exp_growth_rates_PGLS_nothermo), plug in values for fast and slow example species
#only doubling_measured and trna count were significant, so can take out rrna count and unicellular
slowest_growth = log(slow_example$doubling_measured, base = 10)*coef(exp_growth_rates_PGLS_nothermo)['log(doubling_measured, base = 10)']
fastest_growth = log(fast_example$doubling_measured, base = 10)*coef(exp_growth_rates_PGLS_nothermo)['log(doubling_measured, base = 10)']

fast_trna = log(fast_example$trna_count, base = 10)*coef(exp_growth_rates_PGLS_nothermo)['log(trna_count, base = 10)']
slow_trna = log(slow_example$trna_count, base = 10)*coef(exp_growth_rates_PGLS_nothermo)['log(trna_count, base = 10)']

#calculations -- both unicellular so add that
m = coef(exp_growth_rates_PGLS_nothermo)['log(total_gc_coding, base = 10)'][[1]]
b_slow = coef(exp_growth_rates_PGLS_nothermo)[[1]] + slowest_growth + slow_trna
b_fast = coef(exp_growth_rates_PGLS_nothermo)[[1]] + fastest_growth + fast_trna

input_gc_1 = 8.00e06
fast_motifs_1 = 10^(m*(log(input_gc_1, base = 10)) + b_fast[[1]])
slow_motifs_1 = 10^(m*(log(input_gc_1, base = 10)) + b_slow[[1]])
input_gc_2 = 3.5e05
fast_motifs_2 = 10^(m*(log(input_gc_2, base = 10)) + b_fast[[1]])
slow_motifs_2 = 10^(m*(log(input_gc_2, base = 10)) + b_slow[[1]])
(slow_motifs_1-fast_motifs_1)/fast_motifs_1
#slow growth 25% increase

#full plot
pgls_plot <- ggplot(pretty_data_minus_thermo$data, aes(y = total_motifs, x = total_gc_coding)) +
  geom_point(pch = 21, size = 2, alpha = 0.75, color = 'white', fill = '#666666') +
  geom_abline(slope = coef(exp_growth_rates_PGLS_nothermo)['log(total_gc_coding, base = 10)'],
              intercept = coef(exp_growth_rates_PGLS_nothermo)[[1]] +
                slowest_growth + slow_trna,
              color = 'white', size = 1.25) +
  geom_abline(slope = coef(exp_growth_rates_PGLS_nothermo)['log(total_gc_coding, base = 10)'],
              intercept = coef(exp_growth_rates_PGLS_nothermo)[[1]] +
                slowest_growth + slow_trna,
              color = '#3288BD', size = 0.75) +
  geom_abline(slope = coef(exp_growth_rates_PGLS_nothermo)['log(total_gc_coding, base = 10)'],
              intercept = coef(exp_growth_rates_PGLS_nothermo)[[1]] +
                fastest_growth + fast_trna,
              color = 'white', size = 1.25) +
  geom_abline(slope = coef(exp_growth_rates_PGLS_nothermo)['log(total_gc_coding, base = 10)'],
              intercept = coef(exp_growth_rates_PGLS_nothermo)[[1]] +
                fastest_growth + fast_trna,
              color = '#D53E4F', size = 0.75) +
  annotate('rect', xmin=4.5e06, xmax=8.25e06, ymin=6500, ymax= 14000, color = 'black', fill = 'white', alpha = 0) +
  annotate('rect', xmin=3.25e05, xmax=5.85e05, ymin=200, ymax= 400, color = 'black', fill = 'white', alpha = 0) +
  scale_y_log10() +
  scale_x_log10() +
  xlab('') +
  ylab('') + 
  theme_linedraw() +
  theme(legend.position = 'NONE') +
  theme(text = element_text(family='Avenir', size=16, color = '#666666')) +
  theme(axis.text = element_text(family='Avenir', size=16, color = '#666666'))


#points to annotate
pointz <- data.frame(x=c(input_gc_1, input_gc_1, input_gc_2, input_gc_2), 
                     y=c(fast_motifs_1, slow_motifs_1, fast_motifs_2, slow_motifs_2), 
                     type=c('inset2', 'inset2', 'inset1', 'inset1'),
                     growth=c('fast', 'slow', 'fast', 'slow'))

#inset 1 -- low gc
inset1 <- ggplot() +
  geom_point(data = pretty_data_minus_thermo$data, aes(y = total_motifs, x = total_gc_coding), pch = 21, size = 10, alpha = 0.75, color = 'white', fill = '#666666') +
  geom_abline(slope = coef(exp_growth_rates_PGLS_nothermo)['log(total_gc_coding, base = 10)'],
              intercept = coef(exp_growth_rates_PGLS_nothermo)[[1]] +
                slowest_growth + slow_trna,
              color = 'white', size = 3) +
  geom_abline(slope = coef(exp_growth_rates_PGLS_nothermo)['log(total_gc_coding, base = 10)'],
              intercept = coef(exp_growth_rates_PGLS_nothermo)[[1]] +
                fastest_growth + fast_trna,
              color = 'white', size = 3) +
  geom_abline(slope = coef(exp_growth_rates_PGLS_nothermo)['log(total_gc_coding, base = 10)'],
              intercept = coef(exp_growth_rates_PGLS_nothermo)[[1]] +
                slowest_growth + slow_trna,
              color = '#3288BD', size = 2.5) +
  geom_abline(slope = coef(exp_growth_rates_PGLS_nothermo)['log(total_gc_coding, base = 10)'],
              intercept = coef(exp_growth_rates_PGLS_nothermo)[[1]] +
                fastest_growth + fast_trna,
              color = '#D53E4F', size = 2.5) +
  geom_point(data = subset(pointz, type == 'inset1'), aes(y = y, x = x), color = 'white', fill = '#666666', pch = 23, size = 6, alpha = 1) +
  geom_line(data = subset(pointz, type == 'inset1'), aes(y = y, x = x, group = type), color = '#666666', size = 1.75, linetype = 'dashed') +
  geom_text(data = subset(pointz, type == 'inset1' & growth == 'fast'), aes(x=x-5000, y=y+5, label = round(y, 0)), color = '#D53E4F', family='Avenir', size=6, hjust = 1) +
  geom_text(data = subset(pointz, type == 'inset1' & growth == 'slow'), aes(x=x-5000, y=y+5, label = round(y, 0)), color = '#3288BD', family='Avenir', size=6, hjust = 1) +
  scale_x_log10(limits = c(3.25e05, 5.85e05)) +
  scale_y_log10(limits = c(200, 400)) +
  xlab('') +
  ylab('') + 
  theme_linedraw() +
  theme(legend.position = 'NONE') +
  theme(text = element_text(family='Avenir', size=20, color = '#666666')) +
  theme(axis.text = element_text(family='Avenir', size=20, color = '#666666'))

#inset 2 -- high gc
inset2 <- ggplot() +
  geom_point(data = pretty_data_minus_thermo$data, aes(y = total_motifs, x = total_gc_coding), pch = 21, size = 10, alpha = 0.75, color = 'white', fill = '#666666') +
  geom_abline(slope = coef(exp_growth_rates_PGLS_nothermo)['log(total_gc_coding, base = 10)'],
              intercept = coef(exp_growth_rates_PGLS_nothermo)[[1]] +
                slowest_growth + slow_trna,
              color = 'white', size = 3) +
  geom_abline(slope = coef(exp_growth_rates_PGLS_nothermo)['log(total_gc_coding, base = 10)'],
              intercept = coef(exp_growth_rates_PGLS_nothermo)[[1]] +
                fastest_growth + fast_trna,
              color = 'white', size = 3) +
  geom_abline(slope = coef(exp_growth_rates_PGLS_nothermo)['log(total_gc_coding, base = 10)'],
              intercept = coef(exp_growth_rates_PGLS_nothermo)[[1]] +
                slowest_growth + slow_trna,
              color = '#3288BD', size = 2.5) +
  geom_abline(slope = coef(exp_growth_rates_PGLS_nothermo)['log(total_gc_coding, base = 10)'],
              intercept = coef(exp_growth_rates_PGLS_nothermo)[[1]] +
                fastest_growth + fast_trna,
              color = '#D53E4F', size = 2.5) +
  geom_point(data = subset(pointz, type == 'inset2'), aes(y = y, x = x), color = 'white', fill = '#666666', pch = 23, size = 6, alpha = 1) +
  geom_line(data = subset(pointz, type == 'inset2'), aes(y = y, x = x, group = type), color = '#666666', size = 1.75, linetype = 'dashed') +
  geom_text(data = subset(pointz, type == 'inset2' & growth == 'fast'), aes(x=x-100000, y=y+200, label = round(y, 0)), color = '#D53E4F', family='Avenir', size=6, hjust = 1) +
  geom_text(data = subset(pointz, type == 'inset2' & growth == 'slow'), aes(x=x-120000, y=y+200, label = round(y, 0)), color = '#3288BD', family='Avenir', size=6, hjust = 1) +
  scale_x_log10(limits = c(4.5e06, 8.25e06)) +
  scale_y_log10(limits = c(6500, 14000)) +
  xlab('') +
  ylab('') + 
  theme_linedraw() +
  theme(legend.position = 'NONE') +
  theme(text = element_text(family='Avenir', size=20, color = '#666666')) +
  theme(axis.text = element_text(family='Avenir', size=20, color = '#666666'))

#assemble full plot
dummy_plot <- ggplot() + theme_void()
overview <- ggarrange(dummy_plot, pgls_plot, dummy_plot, ncol=1, nrow=3, align = 'v',
                      widths = c(1, 1), heights = c(1, 1))
insets <- ggarrange(inset1, inset2, nrow=2, ncol=1, align = 'hv', 
                    widths = c(0.75, 1), heights = c(1, 1))

final_plot_s2 <- ggarrange(overview, insets, ncol=2, align = 'v',
                        widths = c(0.8, 1.25), heights = c(1, 1))

# ggsave('figure_s2.svg', plot=final_plot_s2, width = 13, height = 10)

#clean up to minimize loss of sanity
rm(dummy_plot, final_plot_s2, inset1, inset2, insets, overview, pgls_plot, pointz, pretty_data_minus_thermo)
rm(b_fast, fast_motifs_1, fast_motifs_2, fast_trna, fastest_growth, m)
rm(b_slow, slow_motifs_1, slow_motifs_2, slow_trna, slowest_growth, input_gc_1, input_gc_2)
rm(slow_example, fast_example, exp_growth_rates_PGLS_nothermo)


