# script that loads the coverage profile and plots a pretty picutre, maybe some additional statistics would be nice
# need a wrapper to run in python maybe

# submitted as commandline parameter
# coverage_file = '/home/fmetge/Documents/work/circRNA/exon_usage/test_outputfolder/MiSeq_A_300BP.coverage_profiles/2:120127688|120175004_NM_020909.txt'
# output_folder = '/home/fmetge/Documents/work/circRNA/exon_usage/test_outputfolder/MiSeq_A_300BP.coverage_pictures/'

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

coverage_file = args[1]
output_folder = args[2]

filename = strsplit(coverage_file, '/')
coverage_track = filename[[1]][length(filename[[1]])]

circle_id = strsplit(coverage_track, '_')[[1]][1]
transcript_name = paste(strsplit(coverage_track, '_')[[1]][2:length(strsplit(coverage_track, '_')[[1]])], collapse = '_')
transcript_name = gsub('.txt', '', transcript_name)


D = read.table(coverage_file, header = T, as.is = T)

png(paste(output_folder, circle_id, '_',transcript_name, '.png', sep = ''))
plot(D$coverage, type = 'h', col = D$exon, main = paste(circle_id, transcript_name, sep = '\n'), xlab = paste('Exon:', min(D$exon), '- Exon:', max(D$exon)), ylab = 'number of reads')
dev.off()
