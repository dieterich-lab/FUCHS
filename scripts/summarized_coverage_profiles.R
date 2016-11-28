#! /usr/bin/env Rscript



# script to assign relative position to coverage profiles of circles

# include checks if clusters are empty, it should not just break but just omit the plots, make number of centers flexible

options(echo=FALSE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

folder = args[1]
# all coverage profiles in a folder
files <- list.files(path=folder, pattern="*.txt", full.names=T, recursive=FALSE)

# need to include checks for correct input data

summary_table = data.frame(matrix(rep(NA, length(files) * 103), ncol = 103))

# iterate over all circle coverage profiles found in given folder
for(f in 1:length(files)){
  D <- read.table(files[f], as.is = T, header = T)
  pos = D$relative_pos_in_circle
  pos = pos/pos[length(pos)]
  pos = round(pos, 2)
  coverage = tapply(D$coverage, pos, mean)
  circle_id = files[f]
  circle_id = strsplit(strsplit(circle_id, '/')[[1]][length(strsplit(circle_id, '/')[[1]])], "[.]")[[1]][1]
  summary_table[f,1]<- circle_id
  summary_table[f,2] <- length(D[,1])
  if(length(coverage) == 101) summary_table[f, -c(1,2)] <- as.numeric(coverage)

}
# delete duplicated and empty lines
summary_table = summary_table[!duplicated(summary_table),]
summary_table = subset(summary_table, !is.na(X3))

# group circles by length
summary_table_short_circle = subset(summary_table, summary_table[,2] < 500)
summary_table_medium_circle = subset(summary_table, summary_table[,2] >= 500 & summary_table[,2] < 1000)
summary_table_long_circle = subset(summary_table, summary_table[,2] >= 1000)



# Plot average coverage profile, group circles by circle length
pdf(paste(folder, 'coverage_profiles.all_circles.pdf', sep = '/'))

# if subset is not empty, plot summarized profile
if(length(summary_table[,1]) > 0){
  plot(colMeans(summary_table[,-c(1,2)], na.rm = T), type = 'l',lwd = 2, main = paste('All circles (', length(summary_table[,1]), ')', sep = '' ), ylab = 'avg. coverage', xlab = 'relative position', cex.axis = 1.7, cex.lab = 1.7, cex.main = 2)
} else {
  print('ERROR no circle in the table')
}
if(length(summary_table_short_circle[,1]) > 0){
  plot(colMeans(summary_table_short_circle[,-c(1,2)], na.rm = T), type = 'l',lwd = 2, main =  paste('Short circles (', length(summary_table_short_circle[,1]), ')', sep = '' ), ylab = 'avg. coverage', xlab = 'relative position', cex.axis = 1.7, cex.lab = 1.7, cex.main = 2)
} else {
  print('ERROR no short circle in the table')
}
if(length(summary_table_medium_circle[,1]) > 0){
  plot(colMeans(summary_table_medium_circle[,-c(1,2)], na.rm = T), type = 'l',lwd = 2, main =  paste('Medium circles (', length(summary_table_medium_circle[,1]), ')', sep = '' ), ylab = 'avg. coverage', xlab = 'relative position', cex.axis = 1.7, cex.lab = 1.7, cex.main = 2)
} else {
  print('ERROR no medium circle in the table')
}
if(length(summary_table_long_circle[,1]) > 0){
  plot(colMeans(summary_table_long_circle[,-c(1,2)], na.rm = T), type = 'l',lwd = 2, main =  paste('Long circles (', length(summary_table_long_circle[,1]), ')', sep = '' ), ylab = 'avg. coverage', xlab = 'relative position', cex.axis = 1.7, cex.lab = 1.7, cex.main = 2)
} else {
  print('ERROR no long circle in the table')
}
dev.off()


# write cluster table -> check for strandedness of hostgene maybe the dip-shift correlates with strand.
# clustering
library(amap)
library(Hmisc)
library(gplots)

# define a function to choose number of centers
choose_centers <- function(D){
  number_of_circles = length(D[,1])
  centers = 1
    if (number_of_circles > 2) centers <- 2
    if (number_of_circles >= 10) centers <- 3
    if (number_of_circles >= 20) centers <- 4
    if (number_of_circles >= 100) centers <- round(number_of_circles/20)
    if (centers >= 10) centers <- 10 
  return(centers)
}



# perhaps also define a function for makeing the clustering, plotting the clusters and writing the tables

# cluster all circles 
centers = choose_centers(summary_table)
print(centers)
if(centers > 1){
  kk = Kmeans(summary_table[,-c(1,2)], centers = centers, method = 'correlation', iter.max = 50, nstart = 100)
  
  pdf(paste(folder,'coverage.clusters.all_circles.pdf', sep = '/'))
    for(i in 1:centers) {
      if(kk$size[i] > 0){
	plot(kk[[2]][i,], main = paste('Cluster Size: ', kk$size[i], sep = '') , type = 'l', lwd = 5, col = 'dodgerblue', ylab = 'avg. coverage', xlab = 'relative position', cex.axis = 1.7, cex.lab = 1.7, cex.main = 2 )
      }
    }
  dev.off()
  circle_to_cluster_association <- data.frame(circle_id = summary_table[,1], length = summary_table[,2], cluster_id = kk$cluster)
  write.table(circle_to_cluster_association, paste(folder,'cluster_association.all_circles.tsv', sep = '/'), sep = '\t', row.names = F, quote = F)
  write.table(kk[[2]], paste(folder,'cluster_means.all_circles.tsv', sep = '/'), sep = '\t', quote = F, row.names = T)

} else {
  print('ERROR not enough circles in table to perform a clustering')
}

# cluster short circles
centers = choose_centers(summary_table_short_circle)
print(centers)
if(centers > 1){
  kk.short = Kmeans(summary_table_short_circle[,-c(1,2)], centers = centers, method = 'correlation', iter.max = 50, nstart = 50)
  
  pdf(paste(folder,'coverage.clusters.short_circles.pdf', sep = '/'))
    for(i in 1:centers) {
      if(kk.short$size[i] > 0){
	plot(kk.short[[2]][i,], main = paste('Cluster Size: ', kk.short$size[i], sep = '') , type = 'l', lwd = 5, col = 'dodgerblue', ylab = 'avg. coverage', xlab = 'relative position', cex.axis = 1.7, cex.lab = 1.7, cex.main = 2 )
      }
    }
  dev.off()
  circle_to_cluster_association.short <- data.frame(circle_id = summary_table_short_circle[,1], length = summary_table_short_circle[,2], cluster_id = kk.short$cluster)
  write.table(circle_to_cluster_association.short, paste(folder,'cluster_association.short_circles.tsv', sep = '/'), sep = '\t', row.names = F, quote = F)
  write.table(kk.short[[2]], paste(folder,'cluster_means.short_circles.tsv', sep = '/'), sep = '\t', quote = F, row.names = T)
} else {
  print('ERROR not enough short circles in table to perform a clustering')
}

# cluster medium circles
centers = choose_centers(summary_table_medium_circle)
if(centers > 1){
  kk.medium = Kmeans(summary_table_medium_circle[,-c(1,2)], centers = centers, method = 'correlation', iter.max = 50, nstart = 100)
  
  pdf(paste(folder,'coverage.clusters.medium_circles.pdf', sep = '/'))
    for(i in 1:centers) {
      if(kk.medium$size[i] >= 0){
	plot(kk.medium[[2]][i,], main = paste('Cluster Size: ', kk.medium$size[i], sep = '') , type = 'l', lwd = 5, col = 'dodgerblue', ylab = 'avg. coverage', xlab = 'relative position', cex.axis = 1.7, cex.lab = 1.7, cex.main = 2 )
      }
    }
  dev.off()
circle_to_cluster_association.medium <- data.frame(circle_id = summary_table_medium_circle[,1], length = summary_table_medium_circle[,2], cluster_id = kk.medium$cluster)
write.table(circle_to_cluster_association.medium, paste(folder,'cluster_association.medium_circles.tsv', sep = '/'), sep = '\t', row.names = F, quote = F)
write.table(kk.medium[[2]], paste(folder,'cluster_means.medium_circles.tsv', sep = '/'), sep = '\t', quote = F, row.names = T)
} else {
  print('ERROR not enough medium circles in table to perform a clustering')
}



# cluster long circles
centers = choose_centers(summary_table_long_circle)
print(centers)
if(centers > 1){
  kk.long = Kmeans(summary_table_long_circle[,-c(1,2)], centers = centers, method = 'correlation', iter.max = 50, nstart = 100)
  
  pdf(paste(folder,'coverage.clusters.long_circles.pdf', sep = '/'))
    for(i in 1:centers) {
      	print(kk.long$size[i])
	if(kk.long$size[i] > 0){
	plot(kk.long[[2]][i,], main = paste('Cluster Size: ', kk.long$size[i], sep = '') , type = 'l', lwd = 5, col = 'dodgerblue', ylab = 'avg. coverage', xlab = 'relative position', cex.axis = 1.7, cex.lab = 1.7, cex.main = 2 )
      }
    }
  dev.off()
  circle_to_cluster_association.long <- data.frame(circle_id = summary_table_long_circle[,1], length = summary_table_long_circle[,2], cluster_id = kk.long$cluster)
  write.table(circle_to_cluster_association.long, paste(folder,'cluster_association.long_circles.tsv', sep = '/'), sep = '\t', row.names = F, quote = F)
  write.table(kk.long[[2]], paste(folder,'cluster_means.long_circles.tsv', sep = '/'), sep = '\t', quote = F, row.names = T)
} else {
  print('ERROR not enough long circles in table to perform a clustering')
}

