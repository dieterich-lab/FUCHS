# python script, to extract circular reads from a bam file, based on a tab-separated circle_id - reads file and a bam file. writes out circle-specific bam files


# input: tab-separated file with circle id as chr:start-end\tlist with read names
# input: bamfile containing circle reads (as well as others)
# output: circle specific bam files, folder for each sample 



import os
import pysam
import argparse

parser = argparse.ArgumentParser(description='Extracts circular reads based on circle_file from the sample.bam and writes them into circle separated bamfiles.')

# input
parser.add_argument('circlefile', metavar = 'sample_circleIDs.txt', help = 'tab separated file chr:start_end(tab)read1,read2,read3.' )
parser.add_argument('bamfile', metavar = 'sample.bam', help = 'bamfile containing chimeric reads, linear reads may be in it but are not required.')
# output
parser.add_argument('outputfolder', metavar = 'folder', help = 'outfolder, there will be a subfolder for the sample containing a bam file for each circle.')
parser.add_argument('sample', metavar = 'sample_name', help = 'sample_name to title every thing.')
# options
parser.add_argument('-r', dest = 'reads', default = 5, type = int, help = 'Circle has to have at least <r> reads to be analysed.')
parser.add_argument('-q', dest = 'mapq', default = 3, type = int, help = 'MAPQ cutoff, only reads passing this threshold will be written to circle bamfile.')


args = parser.parse_args()

# parse arguments
circles = args.circlefile
bamfile = args.bamfile
outfolder = args.outputfolder 
sample = args.sample
cutoff = args.reads
mapq_cutoff = args.mapq

# include some checks to make sure input was provided correctly


# define functions

def read_circles(infile):
    '''
	for each circle, extracts all reads_ids, but also accumulates all circular reads ids independent of circle id.
    '''
    circle_IDs = {}
    reads = {}
    I = open(infile)
    for circle in I:
	if not circle.startswith('#'):
	    circle_IDs[circle.split('\t')[0]] = circle.replace('\n', '').split('\t')[1].split(',') [:-1]
	    for lola in circle.replace('\n', '').split('\t')[1].split(',') [:-1]:
		reads[lola] = 0
    I.close()
    return(circle_IDs, reads)

def load_alignment(infile, circle_reads, cutoff):
    '''
	loads sample bamfile and extracts all circular reads based on read list obtained from read_circles or any given list of read ids.
    '''
    reads = {}
    if 'bam' in infile.split('.')[-1] :
	samfile = pysam.AlignmentFile(bamfile, 'rb')
    else:
	samfile = pysam.AlignmentFile(bamfile, 'r')
    for i,read in enumerate(samfile.fetch()):
	if i % 1000000 == 0:
	    print('%s reads processed' %(i))
	if read.qname in circle_reads and read.mapq > cutoff:
	    #print(float(i)/17000000)
	    if not read.qname in reads:
		reads[read.qname] = {}
	    reads[read.qname][i] = read
    samfile.close()
    return(reads)

def write_circle_bam(reads, circles, cutoff, template, outfolder):
    '''
	for each circle, writes a bam file containing only reads spanning the circle junction and their mates if mates are present.
    '''
    samfile = pysam.AlignmentFile(template, 'rb')
    for circle in circles:
	if len(circles[circle]) >= cutoff:
	    circle_bam = pysam.AlignmentFile("%s/%s_%sreads.bam" %(outfolder, circle.replace(':', '_').replace('|', '_'), len(circles[circle])), "wb", template=samfile)
	    for read in circles[circle]:
		if read in reads:
		    for part in reads[read]:
			circle_bam.write(reads[read][part])
	    circle_bam.close()
    samfile.close()
    return


# run script

circle_info, circle_reads = read_circles(circles)
print('DONE reading circles, found %s circles' %(len(circle_info)))
reads = load_alignment(bamfile, circle_reads, mapq_cutoff)
print('DONE extracting circular reads')
folders = os.listdir(outfolder)
if not sample in folders:
    os.mkdir('%s/%s' %(outfolder, sample))

write_circle_bam(reads, circle_info, cutoff, bamfile, '%s/%s' %(outfolder, sample))
print('DONE writing circle bam files\n')
files = os.listdir('%s/%s' %(outfolder, sample))
print('%s circles passed your thresholds of at least %s reads with at least a mapq of %s\n\n' %(len(files), cutoff, mapq_cutoff))

for f in files:
    if f.split('.')[-1] == 'bam':
	pysam.sort('%s/%s/%s' %(outfolder, sample, f), '%s/%s/%s' %(outfolder, sample, f.replace('.bam', '.sorted')))
	pysam.index('%s/%s/%s.bam' %(outfolder, sample, f.replace('.bam', '.sorted')))
	os.system('rm %s/%s/%s' %(outfolder, sample, f))





