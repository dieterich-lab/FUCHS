# main script to run FUCHS 



import os
import pysam
import argparse

parser = argparse.ArgumentParser(description='')

# input
parser.add_argument('circlefile', metavar = 'sample_circleIDs.txt', help = 'tab separated file chr:start_end(tab)read1,read2,read3.' )
parser.add_argument('bamfile', metavar = 'sample.bam', help = 'bamfile containing chimeric reads, linear reads may be in it but are not required.')
parser.add_argument('bedfile', metavar = 'feature.bed', help = 'bed formatted feature file including exons.' )
# output
parser.add_argument('outputfolder', metavar = 'folder', help = 'outfolder, there will be a subfolder for the sample containing a bam file for each circle.')
parser.add_argument('sample', metavar = 'sample_name', help = 'sample_name to title every thing.')
# options
parser.add_argument('-r', dest = 'reads', default = 5, type = int, help = 'Circle has to have at least <r> reads to be analysed.')
parser.add_argument('-q', dest = 'mapq', default = 3, type = int, help = 'MAPQ cutoff, only reads passing this threshold will be written to circle bamfile.')
parser.add_argument('-e', dest = 'exon_index', default = 3, type = int, help = 'field indicating the exon number after splitting feature name.')
parser.add_argument('-s', dest = 'split_character', default = '_', help = 'feature name separator.')

args = parser.parse_args()

# parse arguments
circles = args.circlefile
bamfile = args.bamfile
bedfile = args.bedfile
outfolder = args.outputfolder 
sample = args.sample
cutoff_reads = args.reads
cutoff_mapq = args.mapq
exon_index = args.exon_index
split_character = args.split_character



# Step1 : extract circle reads from sample bam file
# input: bamfile, circle read file
os.system('python extract_reads.py -r %s -q %s %s %s %s %s' %(cutoff_reads, cutoff_mapq, circles, bamfile, outfolder, sample))
# python extract_reads.py -r 3 -q 5 ../circle_files/MiSeq_test_modules.txt ../alignments/MiSeq_A_300BP_noqual.sorted.bam ../test_outputfolder test_modules

# Step2 : get coverage profile for each circle (one transcript per gene, best if most fitting transcript)
os.system('python get_coverage_profile.py -e %s -s %s %s %s %s' %(exon_index, split_character, bedfile, outfolder, sample))
# python get_coverage_profile.py -e 3 -s _ /home/fmetge/Documents/work/Annotations/hg38/hg38.RefSeq.exons.bed /home/fmetge/Documents/work/circRNA/exon_usage/test_outputfolder/ test_modules


# Step3 : pictures for all circles

files = os.listdir('%s/%s.coverage_profiles' %(outfolder, sample))
os.mkdir('%s/%s.coverage_pictures' %(outfolder, sample))

for f in files:
    os.system('Rscript make_coverage_picture.R %s/%s.coverage_profiles/%s %s/%s.coverage_pictures/' %(outfolder, sample, f, outfolder, sample))


# Step4 : summarize and cluster coverage profiles
os.system('Rscript summarized_coverage_profiles.R %s/%s.coverage_profiles' %(outfolder, sample))
# Step5 : identifiy potential false positives

# Step6 : identify AS if present

## test environment
# python FUCHS.py
