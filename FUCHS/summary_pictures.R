#! /usr/bin/env Rscript


# analyse mate pairs

old_cerebellum = read.table('old_cerebellum.mate_status.genes.txt', as.is = T, header = T)
young_cerebellum = read.table('young_cerebellum.mate_status.genes.txt', as.is = T, header = T)
old_hippocampus = read.table('old_hippocampus.mate_status.genes.txt', as.is = T, header = T)
young_hippocampus = read.table('young_hippocampus.mate_status.genes.txt', as.is = T, header = T)

par(mfrow = c(2,2))
boxplot(old_cerebellum$single/old_cerebellum$num_reads, old_cerebellum$double/old_cerebellum$num_reads, old_cerebellum$undefined/old_cerebellum$num_reads, main = 'Old Cerebellum', names = c('single', 'double', 'rolling'))
abline(h = 0.5, lty = 2)
abline(h = 0.75, lty = 2)
abline(h = 0.25, lty = 2)
boxplot(young_cerebellum$single/young_cerebellum$num_reads, young_cerebellum$double/young_cerebellum$num_reads, young_cerebellum$undefined/young_cerebellum$num_reads, main = 'Young Cerebellum', names = c('single', 'double', 'rolling'))
abline(h = 0.5, lty = 2)
abline(h = 0.75, lty = 2)
abline(h = 0.25, lty = 2)
boxplot(old_hippocampus$single/old_hippocampus$num_reads, old_hippocampus$double/old_hippocampus$num_reads, old_hippocampus$undefined/old_hippocampus$num_reads, main = 'Old Hippocampus', names = c('single', 'double', 'rolling'))
abline(h = 0.5, lty = 2)
abline(h = 0.75, lty = 2)
abline(h = 0.25, lty = 2)
boxplot(young_hippocampus$single/young_hippocampus$num_reads, young_hippocampus$double/young_hippocampus$num_reads, young_hippocampus$undefined/young_hippocampus$num_reads, main = 'Young Hippocampus', names = c('single', 'double', 'rolling'))
abline(h = 0.5, lty = 2)
abline(h = 0.75, lty = 2)
abline(h = 0.25, lty = 2)

boxplot(old_cerebellum$min_length, young_cerebellum$min_length, old_hippocampus$min_length, young_hippocampus$min_length, log = 'y', main = 'Circle length', names = c('oCB', 'yCB', 'oHC', 'yHC'))
abline(h = 500, lty = 2)

boxplot(old_cerebellum$max_length, young_cerebellum$max_length, old_hippocampus$max_length, young_hippocampus$max_length, log = 'y', main = 'Circle length', names = c('oCB', 'yCB', 'oHC', 'yHC'))
abline(h = 500, lty = 2)
# significant difference between young and old but not between different brain regions

warnings()