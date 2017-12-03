#! /usr/bin/env Rscript


# script that loads the coverage profile and plots a pretty picutre, maybe some additional statistics would be nice
# need a wrapper to run in python maybe

# submitted as commandline parameter

#Sys.setenv("DISPLAY"=":0")
#options(echo=FALSE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

coverage_file = args[1]
output_folder = args[2]

# define functions

smoothing <- function(x){
  # X being the conservation column smooth conservation track.
  # Take average of region centered around current base +/- xBP depending on gene_length
  w = round(length(x) * 0.01)
  smoothed = rep(NA, length(x))
  for(i in 1:length(x)){
    if(i < w){
      smoothed[i] = mean(x[1: (i+w)])
    }
    else{
      smoothed[i] = mean(x[(i-w): (i+w)], na.rm = TRUE)
    }
  }
  return(smoothed)
}

filename = strsplit(coverage_file, '/')
coverage_track = filename[[1]][length(filename[[1]])]

circle_id = strsplit(coverage_track, '[.]')[[1]][1]
transcript_name = strsplit(coverage_track, '[.]')[[1]][2]

D = read.table(coverage_file, header = T, as.is = T)
smoothed = smoothing(D$coverage)

pdf(paste(output_folder, circle_id, '_',transcript_name, '.pdf', sep = ''), title=paste(circle_id,transcript_name,sep=" "))
  plot(smoothed,
        type = 'h',
        col = D$exon + 1,
        main = paste(circle_id, transcript_name, sep = '\n'),
        xlab = paste('Exon:', min(D$exon), '- Exon:',
        max(D$exon)),
        ylab = 'number of reads'
        )
invisible(dev.off()) # hide NULL device printout