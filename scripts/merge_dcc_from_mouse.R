# merge and filter DCC calls

FC = read.table('frontal_cortex/CircRNACount', as.is = T, header = T)
OCB = read.table('old_cerebellum/CircRNACount', as.is = T, header = T)
OHE = read.table('old_heart/CircRNACount', as.is = T, header = T)
OHC = read.table('old_hippocampus/CircRNACount', as.is = T, header = T)
OLI = read.table('old_liver/CircRNACount', as.is = T, header = T)
YCB = read.table('young_cerebellum/CircRNACount', as.is = T, header = T)
YHE = read.table('young_heart/CircRNACount', as.is = T, header = T)
YHC = read.table('young_hippocampus/CircRNACount', as.is = T, header = T)
YLI = read.table('young_liver/CircRNACount', as.is = T, header = T)

FC = aggregate(FC$Chimeric.out.junction, by=list(paste(FC$Chr, FC$Start, FC$End)), FUN=sum)
OCB = aggregate(OCB$Chimeric.out.junction, by=list(paste(OCB$Chr, OCB$Start, OCB$End)), FUN=sum)
OHE = aggregate(OHE$Chimeric.out.junction, by=list(paste(OHE$Chr, OHE$Start, OHE$End)), FUN=sum)
OHC = aggregate(OHC$Chimeric.out.junction, by=list(paste(OHC$Chr, OHC$Start, OHC$End)), FUN=sum)
OLI = aggregate(OLI$Chimeric.out.junction, by=list(paste(OLI$Chr, OLI$Start, OLI$End)), FUN=sum)
YCB = aggregate(YCB$Chimeric.out.junction, by=list(paste(YCB$Chr, YCB$Start, YCB$End)), FUN=sum)
YHE = aggregate(YHE$Chimeric.out.junction, by=list(paste(YHE$Chr, YHE$Start, YHE$End)), FUN=sum)
YHC = aggregate(YHC$Chimeric.out.junction, by=list(paste(YHC$Chr, YHC$Start, YHC$End)), FUN=sum)
YLI = aggregate(YLI$Chimeric.out.junction, by=list(paste(YLI$Chr, YLI$Start, YLI$End)), FUN=sum)



CircRNACount = merge(YCB, YHC, by = 'Group.1', all = T)
names(CircRNACount) <- c('Group.1', 'young_cerebellum', 'young_hippocampus')
CircRNACount = merge(CircRNACount, YHE, by = 'Group.1', all = T)
names(CircRNACount) <- c('Group.1', 'young_cerebellum', 'young_hippocampus', 'young_heart')
CircRNACount = merge(CircRNACount, YLI, by = 'Group.1', all = T)
names(CircRNACount) <- c('Group.1', 'young_cerebellum', 'young_hippocampus', 'young_heart', 'young_liver')
CircRNACount = merge(CircRNACount, OCB, by = 'Group.1', all = T)
names(CircRNACount) <- c('Group.1', 'young_cerebellum', 'young_hippocampus', 'young_heart', 'young_liver', 'old_cerebellum')
CircRNACount = merge(CircRNACount, OHC, by = 'Group.1', all = T)
names(CircRNACount) <- c('Group.1', 'young_cerebellum', 'young_hippocampus', 'young_heart', 'young_liver', 'old_cerebellum', 'old_hippocampus')
CircRNACount = merge(CircRNACount, OHE, by = 'Group.1', all = T)
names(CircRNACount) <- c('Group.1', 'young_cerebellum', 'young_hippocampus', 'young_heart', 'young_liver', 'old_cerebellum', 'old_hippocampus', 'old_heart')
CircRNACount = merge(CircRNACount, OLI, by = 'Group.1', all = T)
names(CircRNACount) <- c('Group.1', 'young_cerebellum', 'young_hippocampus', 'young_heart', 'young_liver', 'old_cerebellum', 'old_hippocampus', 'old_heart', 'old_liver')
CircRNACount = merge(CircRNACount, FC, by = 'Group.1', all = T)
names(CircRNACount) <- c('Group.1', 'young_cerebellum', 'young_hippocampus', 'young_heart', 'young_liver', 'old_cerebellum', 'old_hippocampus', 'old_heart', 'old_liver', 'frontal_cortex')


CircRNACount[is.na(CircRNACount)] <- 0
Filter = rowSums(CircRNACount[,2:10]>1)

fCircRNACount = subset(CircRNACount, Filter > 3)

write.table(fCircRNACount, 'CircRNACount', sep = '\t', quote = F, row.names = F)


