# script to assign relative position to coverage profiles of circles


options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

folder = args[1]
# all coverage profiles in a folder
files <- list.files(path=folder, pattern="*.txt", full.names=T, recursive=FALSE)

# need to include checks for correct input data

summary_table = data.frame(matrix(rep(NA, length(files) * 103), ncol = 103))

for(f in 1:length(files)){
  D <- read.table(files[f], as.is = T, header = T)
  pos = D$relative_pos_in_circle
  pos = pos/pos[length(pos)]
  pos = round(pos, 2)
  coverage = tapply(D$coverage, pos, mean)
  circle_id = files[f]
  circle_id = strsplit(strsplit(circle_id, '/')[[1]][-1], "[.]")[[1]][1]
  summary_table[f,1]<- circle_id
  summary_table[f,2] <- length(D[,1])
  if(length(coverage) == 101) summary_table[f, -c(1,2)] <- as.numeric(coverage)

}
summary_table = summary_table[!duplicated(summary_table),]
summary_table = subset(summary_table, !is.na(X3))
summary_table_short_circle = subset(summary_table, summary_table[,2] < 500)
summary_table_medium_circle = subset(summary_table, summary_table[,2] >= 500 & summary_table[,2] < 1000)
summary_table_long_circle = subset(summary_table, summary_table[,2] >= 1000)



# make them prettier
pdf(paste(folder, 'coverage_profiles.all_circles.pdf', sep = '/'))
#par(mfrow = c(2,2))
plot(colMeans(summary_table[,-c(1,2)], na.rm = T), type = 'l',lwd = 2, main = 'All circles', ylab = 'avg. coverage', xlab = 'relative position', cex.axis = 1.7, cex.lab = 1.7, cex.main = 2)
plot(colMeans(summary_table_short_circle[,-c(1,2)], na.rm = T), type = 'l',lwd = 2, main = 'Short circles', ylab = 'avg. coverage', xlab = 'relative position', cex.axis = 1.7, cex.lab = 1.7, cex.main = 2)
plot(colMeans(summary_table_medium_circle[,-c(1,2)], na.rm = T), type = 'l',lwd = 2, main = 'Medium circles', ylab = 'avg. coverage', xlab = 'relative position', cex.axis = 1.7, cex.lab = 1.7, cex.main = 2)
plot(colMeans(summary_table_long_circle[,-c(1,2)], na.rm = T), type = 'l',lwd = 2, main = 'Long circles', ylab = 'avg. coverage', xlab = 'relative position', cex.axis = 1.7, cex.lab = 1.7, cex.main = 2)
dev.off()


# write cluster table -> check for strandedness of hostgene maybe the dip-shift correlates with strand.
# clustering
library(amap)
library(Hmisc)
library(gplots)

# good clustering
kk = Kmeans(summary_table[,-c(1,2)], centers = 8, method = 'correlation', iter.max = 50, nstart = 100)

par(mfrow = c(3,3))
pdf(paste(folder,'coverage.clusters.all_circles.pdf', sep = '/'))
for(i in 1:8) plot(kk[[2]][i,], main = paste('Cluster Size: ', kk$size[i], sep = '') , type = 'l', lwd = 5, col = 'dodgerblue', ylab = 'avg. coverage', xlab = 'relative position', cex.axis = 1.7, cex.lab = 1.7, cex.main = 2 )
dev.off()
# good clustering
kk.long = Kmeans(summary_table_long_circle[,-c(1,2)], centers = 5, method = 'correlation', iter.max = 50, nstart = 100)
# par(mfrow = c(3,4))
pdf(paste(folder,'coverage.clusters.long_circles.pdf', sep = '/'))
for(i in 1:5) plot(kk.long[[2]][i,], main = kk.long$size[i] , type = 'l', lwd = 5, col = 'dodgerblue', ylab = 'avg. coverage', xlab = 'relative position', cex.axis = 1.7, cex.lab = 1.7, cex.main = 2 )
dev.off()

# fine
kk.short = Kmeans(summary_table_short_circle[,-c(1,2)], centers = 3, method = 'correlation', iter.max = 50, nstart = 50)
# par(mfrow = c(3,4))
pdf(paste(folder,'coverage.clusters.short_circles.pdf', sep = '/'))
for(i in 1:3) plot(kk.short[[2]][i,], main = kk.short$size[i] , type = 'l', lwd = 5, col = 'dodgerblue', ylab = 'avg. coverage', xlab = 'relative position', cex.axis = 1.7, cex.lab = 1.7, cex.main = 2 )
dev.off()

# ok
kk.medium = Kmeans(summary_table_medium_circle[,-c(1,2)], centers = 5, method = 'correlation', iter.max = 50, nstart = 100)
# par(mfrow = c(3,4))
pdf(paste(folder,'coverage.clusters.medium_circles.pdf', sep = '/'))
for(i in 1:5) plot(kk.medium[[2]][i,], main = kk.medium$size[i] , type = 'l', lwd = 5, col = 'dodgerblue', ylab = 'avg. coverage', xlab = 'relative position', cex.axis = 1.7, cex.lab = 1.7, cex.main = 2 )
dev.off()

# write out cluster tables

circle_to_cluster_association <- data.frame(circle_id = summary_table[,1], length = summary_table[,2], cluster_id = kk$cluster)
circle_to_cluster_association.short <- data.frame(circle_id = summary_table_short_circle[,1], length = summary_table_short_circle[,2], cluster_id = kk.short$cluster)
circle_to_cluster_association.long <- data.frame(circle_id = summary_table_long_circle[,1], length = summary_table_long_circle[,2], cluster_id = kk.long$cluster)
circle_to_cluster_association.medium <- data.frame(circle_id = summary_table_medium_circle[,1], length = summary_table_medium_circle[,2], cluster_id = kk.medium$cluster)

write.table(circle_to_cluster_association, paste(folder,'cluster_association.all_circles.tsv', sep = '/'), sep = '\t', row.names = F, quote = F)
write.table(circle_to_cluster_association.short, paste(folder,'cluster_association.short_circles.tsv', sep = '/'), sep = '\t', row.names = F, quote = F)
write.table(circle_to_cluster_association.medium, paste(folder,'cluster_association.medium_circles.tsv', sep = '/'), sep = '\t', row.names = F, quote = F)
write.table(circle_to_cluster_association.long, paste(folder,'cluster_association.long_circles.tsv', sep = '/'), sep = '\t', row.names = F, quote = F)

write.table(kk[[2]], paste(folder,'cluster_means.all_circles.tsv', sep = '/'), sep = '\t', quote = F, row.names = T)
write.table(kk.short[[2]], paste(folder,'cluster_means.short_circles.tsv', sep = '/'), sep = '\t', quote = F, row.names = T)
write.table(kk.long[[2]], paste(folder,'cluster_means.long_circles.tsv', sep = '/'), sep = '\t', quote = F, row.names = T)
write.table(kk.medium[[2]], paste(folder,'cluster_means.medium_circles.tsv', sep = '/'), sep = '\t', quote = F, row.names = T)


