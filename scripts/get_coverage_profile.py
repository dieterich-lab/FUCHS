# python script to get coverage profile for circle

import os
import argparse
import pybedtools

parser = argparse.ArgumentParser(description='')

# input
parser.add_argument('bedfile', metavar = 'feature.bed', help = 'bed formatted feature file including exons.' )
parser.add_argument('inputfolder', metavar = 'folder', help = 'folder containing all circle bam files. (full path, but without sample name)')
parser.add_argument('sample', metavar = 'sample_name', help = 'sample_name to title every thing.')
# options
parser.add_argument('-e', dest = 'exon_index', default = 3, type = int, help = 'field indicating the exon number after splitting feature name.')
parser.add_argument('-s', dest = 'split_character', default = '_', help = 'feature name separator.')

args = parser.parse_args()

# parse arguments
bedfile = args.bedfile
inputfolder = args.inputfolder 
sample = args.sample
exon_index = args.exon_index
split_character = args.split_character

# include some checks to make sure input was provided correctly


# define functions

def circle_exon_count(bamfile2, bedfile, exon_index, split_character): # does what I think it does, adjust to collapse different transcripts from the same gene, choose transcript describing the circle best
    '''
    '''
    print('before loading bam')
    print(bamfile2)
    x = pybedtools.example_bedtool(bamfile2)
    print('bam loaded')
    b = pybedtools.example_bedtool(bedfile)
    y = x.intersect(b, bed = True, wo = True)
    transcripts = {}
    found_features = []
    for hit in y:
	found_features += [hit[15]]
	transcript = hit[15]
	start = int(hit[6])
	end = int(hit[7])
	length = int((hit[10]).replace(',', ''))
	strand_read = hit[5]
	strand_feature = hit[17]
	transcript_id = split_character.join(transcript.split(split_character)[0:2])
	exon = int(transcript.split(split_character)[exon_index])
	read = hit[3]
	chromosome = hit[0]
	if not transcript_id in transcripts:
	    transcripts[transcript_id] = {}
	if not exon in transcripts[transcript_id]:
	    transcripts[transcript_id][exon] = {'length': length , 'start': start, 'end': end, 'strand_read': [], 'strand_feature': strand_feature, 'reads': [], 'chromosome': chromosome}
	transcripts[transcript_id][exon]['reads'] += [read]
	transcripts[transcript_id][exon]['strand_read'] += [strand_read]
    return(transcripts, found_features)

def write_exon_count(outfile, exon_count, sample, circle_id): #append to existing exon_count file for the sample
    '''
    '''
    out = open(outfile, 'a')
    # sample\tcircle_id\ttranscript_id\texon_id\tchr\tstart\tend\tstrand\texon_length\tunique_reads\tfragments\tnumber+\tnumber-\n
    # sort exon ids per transcript..and then iterate from min to max, if one exon isn't in it, fill with 0's, identify potentially skipped exons
    for transcript in exon_count:
	for exon in range(min(exon_count[transcript]), (max(exon_count[transcript])+1)):
	    if exon in exon_count[transcript]:
		num_plus = exon_count[transcript][exon]['strand_read'].count('+')
		num_minus = exon_count[transcript][exon]['strand_read'].count('-')
		unique_reads = set([w.split('/')[0] for w in exon_count[transcript][exon]['reads']])
		out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(sample, circle_id, transcript, exon, exon_count[transcript][exon]['chromosome'], exon_count[transcript][exon]['start'], exon_count[transcript][exon]['end'], exon_count[transcript][exon]['strand_feature'], exon_count[transcript][exon]['length'], len(unique_reads), len(exon_count[transcript][exon]['reads']), num_plus, num_minus))
	    else:
		out.write('%s\t%s\t%s\t%s\t0\t0\t0\t0\t0\t0\t0\t0\t0\n' %(sample, circle_id, transcript, exon))
    return


def filter_features(bed_features, feature_names):
    '''
    '''
    intervals = []
    for interval in bed_features:
	if interval[3] in feature_names:
	    intervals += [interval]
    return(intervals)

def circle_coverage_profile(bamfile, bedfile, exon_ind, split_character):
    '''
    '''
    x = pybedtools.example_bedtool(bamfile)
    y = x.coverage(bedfile, d = True)
    transcriptwise_coverage = {}
    for position in y:
	transcript = split_character.join(position[3].split(split_character)[0:2])
	exon = int(position[3].split(split_character)[exon_ind])
	if not transcript in transcriptwise_coverage:
	    transcriptwise_coverage[transcript] = {}
	if not exon in transcriptwise_coverage[transcript]:
	    transcriptwise_coverage[transcript][exon] = {'relative_positions': [], 'position_coverage' : [], 'chromosome': position[0], 'start': position[1], 'end': position[2]}
	transcriptwise_coverage[transcript][exon]['position_coverage'] += [position[7]]
	transcriptwise_coverage[transcript][exon]['relative_positions'] += [position[6]]
    return(transcriptwise_coverage)

def write_coverage_profile(inputfolder, coverage_profile, sample, circle_id):
    '''
    '''
    for transcript in coverage_profile:
	out = open('%s/%s.coverage_profiles/%s.%s.txt' %(inputfolder, sample, circle_id, transcript), 'w')
	out.write('exon\trelative_pos_in_circle\trelative_pos_in_exon\tcoverage\n')
	pos_in_circle = 1
	for exon in range(min(coverage_profile[transcript]), (max(coverage_profile[transcript])+1)):
	    if exon in coverage_profile[transcript]:
		for i, position in enumerate(coverage_profile[transcript][exon]['relative_positions']):
		    out.write('%s\t%s\t%s\t%s\n' %(exon, pos_in_circle, position, coverage_profile[transcript][exon]['position_coverage'][i]))
		    pos_in_circle += 1
	    else:
		for i, position in enumerate(range(0,50)):
		    out.write('%s\t%s\t%s\t0\n' %(exon, pos_in_circle, position))
		    pos_in_circle += 1
    out.close()
    return



# circle exon count over all bam files in sample folder, this could easily be parralellised

# initializing the result table file
exon_count_file = '%s/%s.exon_counts.txt' %(inputfolder, sample)
exon_counts_out = open(exon_count_file, 'w')
exon_counts_out.write('sample\tcircle_id\ttranscript_id\texon_id\tchr\tstart\tend\tstrand\texon_length\tunique_reads\tfragments\tnumber+\tnumber-\n')
exon_counts_out.close()

# all circle files in a given folder
files = os.listdir('%s/%s' %(inputfolder, sample))

# create folder for coverage profiles
folders = os.listdir(inputfolder)
print(folders)
if not '%s.coverage_profiles' %(sample) in folders:
    os.mkdir('%s/%s.coverage_profiles' %(inputfolder, sample))

# iterate over all files
for f in files:
    # only consider sorted bam files
    if f.split('.')[-2] == 'sorted':
	# extract circle id from filename, works for files generated by extract_reads.py, consider making this more flexible
	circle_id = '%s_%s_%s' %(f.split('_')[0], f.split('_')[1], f.split('_')[2])
	bamfile2 = '%s/%s/%s' %(inputfolder, sample, f)
	print(bamfile2)
	# open bed feature file
	b = pybedtools.example_bedtool(bedfile)
	# get read counts for each exon in circle
	exon_counts, found_features = circle_exon_count(bamfile2, bedfile, exon_index, split_character)
	# add circle to result table
	write_exon_count(exon_count_file, exon_counts, sample, circle_id)
	filtered_features = filter_features(b, found_features)
	print('.')
	if len(filtered_features) > 0:
	    coverage_track = circle_coverage_profile(bamfile2, filtered_features, exon_index, split_character)
	    write_coverage_profile(inputfolder, coverage_track, sample, circle_id)


## make pictures using rscript

#### test area
#bedfile = '/home/fmetge/Documents/work/Annotations/hg38/hg38.RefSeq.exons.bed'
#bamfile = '/home/fmetge/Documents/work/circRNA/exon_usage/test_outputfolder/MiSeq_A_300BP/2_199368605_199433514_7reads.sorted.bam'
#bamfile = '/home/fmetge/Documents/work/circRNA/exon_usage/test_outputfolder/MiSeq_A_300BP/10_87223718_87355425_9reads.sorted.bam'


#x = pybedtools.example_bedtool(bamfile)
#b = pybedtools.example_bedtool(bedfile)
#exon_counts, found_features = circle_exon_count(bamfile, bedfile, 3)
#filtered_features = filter_features(b, found_features)
##y = x.coverage(filtered_features, d = True)

#coverage_track = circle_coverage_profile(bamfile, filtered_features, 3)


