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
parser.add_argument('-p', dest = 'ref_platform', default = 'refseq', help = 'specifies the annotation platform which was used (refseq or ensembl)')

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
platform = args.ref_platform


# test for correct input data
accepted_platforms = ('refseq', 'ensembl')
platform = platform.lower()
if not platform in accepted_platforms:
    print('ERROR please specify an accepted annotation platform. Possible options are: refseq or ensembl')
    quit()

# Step1 : extract circle reads from sample bam file
os.system('python extract_reads.py -r %s -q %s %s %s %s %s' %(cutoff_reads, cutoff_mapq, circles, bamfile, outfolder, sample))

# Step2 : get coverage profile for each circle (one transcript per gene, best if most fitting transcript)
os.system('python get_coverage_profile.py -e %s -s %s -p %s %s %s %s' %(exon_index, split_character, platform, bedfile, outfolder, sample))

# Step3 : circle oriented file containing information about single and double breakpoint reads
os.system('python get_mate_information.py -p %s -s %s -a %s %s/%s %s/%s.mate_status.txt' %(platform, split_character, bedfile, outfolder, sample, outfolder, sample ))

# Step4 : pictures for all circles
files = os.listdir('%s/%s.coverage_profiles' %(outfolder, sample))
os.mkdir('%s/%s.coverage_pictures' %(outfolder, sample))

for f in files:
    os.system('Rscript make_coverage_picture.R %s/%s.coverage_profiles/%s %s/%s.coverage_pictures/' %(outfolder, sample, f, outfolder, sample))

# Step5 : summarize and cluster coverage profiles
os.system('Rscript summarized_coverage_profiles.R %s/%s.coverage_profiles' %(outfolder, sample))

# Step6 : identifiy potential false positives

# Step7 : identify AS if present

## test environment
# python FUCHS.py /home/fmetge/Documents/work/circRNA/FUCHS/testdata/github_testdata.circles.txt /home/fmetge/Documents/work/circRNA/FUCHS/testdata/github_testdata.bam /home/fmetge/Documents/work/Annotations/hg38/hg38.RefSeq.exons.bed /home/fmetge/Documents/work/circRNA/FUCHS/testdata/output/ github_testdata_20160111

