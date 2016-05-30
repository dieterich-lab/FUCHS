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
parser.add_argument('-sS', dest = 'skipped_steps', default = 'none', help = 'comma separated list of steps that should be skipped (e.g. step3,step4,step6)')
parser.add_argument('-j', dest = 'star_junction_file', default = 'none', help = 'if you mapped with star and are using step1 you need to provide the junction file here')
parser.add_argument('-m', dest = 'mates', default = 'none', help = 'if you mapped with star, have paired end data and are using step1 you need to provide the mate2 junction file here')
parser.add_argument('-c', dest = 'circle_ids', default = 'none', help = 'if you mapped with star and are using step1 you need to provide a list of circle ids (CircCoordinates from DCC)')

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
skipped_steps = args.skipped_steps.split(',')
mates = args.mates
junctionfile = args.star_junction_file
circle_ids = args.circle_ids


# test for correct input data
accepted_platforms = ('refseq', 'ensembl')
platform = platform.lower()
if not platform in accepted_platforms:
    print('ERROR please specify an accepted annotation platform. Possible options are: refseq or ensembl')
    quit()

print(skipped_steps)

# Step 1: (optional) if DCC was used, extract circle read names from junction file 
if not 'step1' in skipped_steps:
    circles = '%s.reads.txt' %(junctionfile)
    if not mates == 'none': 
	os.system('python get_readnames_from_DCC.py -m %s %s %s' %(mates, circle_ids, junctionfile))
    else:
	os.system('python get_readnames_from_DCC.py %s %s' %(circle_ids, junctionfile))


# Step2 : extract circle reads from sample bam file
os.system('python extract_reads.py -r %s -q %s %s %s %s %s' %(cutoff_reads, cutoff_mapq, circles, bamfile, outfolder, sample))

# Step3 : (optional) get information about possibly rolling circles 
if not 'step3' in skipped_steps:
    os.system('python get_mate_information.py -p %s -s %s -a %s %s/%s %s/%s.mate_status.txt' %(platform, split_character, bedfile, outfolder, sample, outfolder, sample ))

# Step4 : (optional) find exon skipping events
if not 'step4' in skipped_steps:
    os.system('python detect_skipped_exons.py %s/%s %s %s/%s.skipped_exons.txt' %(outfolder, sample, bedfile, outfolder, sample))

# Step5 : (optional) identify different circles within the same host gene
if not 'step5' in skipped_steps:
    os.system('python detect_splicing_variants.py -s %s -p %s %s %s %s/%s.alternative_splicing.txt' %(split_character, platform ,circles, bedfile, outfolder, sample))

# Step6 : (optional) generate coverage profile for each circle (one transcript per gene, best if most fitting transcript)
if not 'step6' in skipped_steps:
    os.system('python get_coverage_profile.py -e %s -s %s -p %s %s %s %s' %(exon_index, split_character, platform, bedfile, outfolder, sample))

# Step7 : (optionl, reqiures step 5)
if not 'step7' in skipped_steps:
    if not 'step6' in skipped_steps:
	os.system('Rscript summarized_coverage_profiles.R %s/%s.coverage_profiles' %(outfolder, sample))
    else:
	print('You are trying cluster the coverage profiles without generating coverage profiles first, please run step 5')


# Step8 : (optional, requires step6) pictures for all circles
if not 'step8' in skipped_steps:
    if not 'step6' in skipped_steps:
	files = os.listdir('%s/%s.coverage_profiles' %(outfolder, sample))
	folders = os.listdir(outfolder)
	if not '%s.coverage_pictures' %(sample) in folders:
	    os.mkdir('%s/%s.coverage_pictures' %(outfolder, sample))
	for f in files:
	    os.system('Rscript make_coverage_picture.R %s/%s.coverage_profiles/%s %s/%s.coverage_pictures/' %(outfolder, sample, f, outfolder, sample))
    else:
	print('You are trying to generate coverage pictures without generating coverage profiles, please run step 5')


os.system('rm /beegfs/group_dv/home/FMetge/tmp/*')


# Step9: (optional) summary pictures, needs to be implemented ;)

# Step10 : identifiy potential false positives (needs to be implemented and tested in the lab, but potentially very interesting)

## test environment
# python FUCHS.py /home/fmetge/Documents/work/circRNA/FUCHS/testdata/github_testdata.circles.txt /home/fmetge/Documents/work/circRNA/FUCHS/testdata/github_testdata.bam /home/fmetge/Documents/work/Annotations/hg38/hg38.RefSeq.exons.bed /home/fmetge/Documents/work/circRNA/FUCHS/testdata/output/ github_testdata_20160111
