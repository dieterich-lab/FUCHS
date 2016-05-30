# script to identify circles where both mates map over the junction

import pysam
import os
import argparse
import pybedtools
import tempfile

parser = argparse.ArgumentParser(description='Extracts mate information and identify singe and double breakpoint fragments')

# input
parser.add_argument('bamfolder', metavar = 'PATH', help = 'path to folder containing circle bamfiles' )
parser.add_argument('outfile', metavar = 'outfile', help = 'path and filename to write the output to' )
parser.add_argument('-a', dest = 'bedfile' ,default = 'none', help = 'if specified, the program will try to infer the circle length without internal introns')
parser.add_argument('-p', dest = 'ref_platform', default = 'refseq', help = 'specifies the annotation platform which was used (refseq or ensembl)')
parser.add_argument('-s', dest = 'split_character', default = '_', help = 'specifies the separator within the name column in bedfile')

args = parser.parse_args()

# parse arguments
bamfolder = args.bamfolder
outfile = args.outfile
bedfile = args.bedfile
platform = args.ref_platform
split_character = args.split_character

# define functions
def get_reads_from_bamfile(bamfile, circle_coordinates):
    mates = {}
    non_junction_fragments = []
    circle = pysam.AlignmentFile(bamfile, "rb")
    for lola in circle:
	name = lola.query_name
	reverse = lola.is_reverse
	start = lola.reference_start
	end = lola.reference_end
	if not name in mates:
	    mates[name] = {'forward' : {'start' : [], 'end' : []}, 'reverse' : {'start' : [], 'end' : []} }
	if reverse and end == circle_coordinates[2]:
	    mates[name]['reverse']['end'] += [start]
	elif reverse and start == circle_coordinates[1] - 1:
	    mates[name]['reverse']['start'] += [end]
	elif end == circle_coordinates[2] and not reverse:
	    mates[name]['forward']['end'] += [start]
	elif start == circle_coordinates[1] - 1 and not reverse:
	    mates[name]['forward']['start'] += [end]
	else:
	    non_junction_fragments += [lola]
    circle.close()
    return(mates, non_junction_fragments)

def classify_reads(mates):
    for lola in mates:
	strands = 0
	for strand in mates[lola]:
	    if len(mates[lola][strand]['start']) == 1 and len(mates[lola][strand]['end']) == 1:
		strands += 1
	if strands == 1:
	    mates[lola]['status'] = 'single'
	elif strands == 2:
	    mates[lola]['status'] = 'double'
	else:
	    mates[lola]['status'] = 'undefined'
    return(mates)

def get_statistics(mates):
    stats = {'single': 0, 'double': 0, 'undefined': 0}
    for lola in mates:
	stats[mates[lola]['status']] += 1
    return(stats)

def annotate_circle(circle_coordinates, bedfile, platform, split_character):
    circle = pybedtools.BedTool('%s %s %s' %(circle_coordinates[0], circle_coordinates[1], circle_coordinates[2]), from_string=True)
    exons = pybedtools.example_bedtool(bedfile)
    features = exons.intersect(circle)
    lengths = {}
    for lola in features:
	if platform == 'refseq':
	    transcript_name = split_character.join(lola[3].split(split_character)[0:2])
	elif platform == 'ensembl':
	    transcript_name = lola[3].split(split_character)[0]
	else:
	    transcript_name = 'NA'
	    print('you are using an unkown reference platform. Please choose between refseq or ensembl')
	length = int(lola[2]) - int(lola[1])
	if not transcript_name in lengths:
	    lengths[transcript_name] = 0
	lengths[transcript_name] += length    
    return(lengths)

def iterate_over_folder(inputfolder, bedfile, platform, split_character):
    results = {}
    files = os.listdir(inputfolder)
    for lola in files:
	if lola.split('.')[-1] == 'bam':
	    print(lola)
	    circle_coordinates = ['_'.join(lola.split('_')[0:-3]), int(lola.split('_')[-3]), int(lola.split('_')[-2])]
	    num_reads = int(lola.split('_')[-1].split('.')[0]. replace('reads', ''))
	    MATES, FRAGMENTS = get_reads_from_bamfile('%s/%s' %(inputfolder, lola), circle_coordinates)
	    MATES = classify_reads(MATES)
	    if not bedfile == 'none':
		LENGTH = annotate_circle(circle_coordinates, bedfile, platform, split_character)
	    else:
		LENGTH = {}
	    STATS = get_statistics(MATES)
	    results[lola.split('.')[0]] = STATS  
	    if len(LENGTH) > 0:
		results[lola.split('.')[0]]['min_length'] = min(LENGTH.items(), key=lambda x: x[1])[1]
		results[lola.split('.')[0]]['max_length'] = max(LENGTH.items(), key=lambda x: x[1])[1]
		results[lola.split('.')[0]]['transcript_ids'] = ','.join(LENGTH.keys())
	    else:
		results[lola.split('.')[0]]['min_length'] = circle_coordinates[2] - circle_coordinates[1]
		results[lola.split('.')[0]]['max_length'] = circle_coordinates[2] - circle_coordinates[1]
		results[lola.split('.')[0]]['transcript_ids'] = 'not_annotated'
	    results[lola.split('.')[0]]['circle_id'] = '%s_%s_%s' %(circle_coordinates[0], circle_coordinates[1], circle_coordinates[2])
	    results[lola.split('.')[0]]['num_reads'] = num_reads
    return(results)


def write_results(results, outfile):
    O = open(outfile, 'w')
    O.write('circle_id\ttranscript_ids\tnum_reads\tmin_length\tmax_length\tsingle\tdouble\tundefined\n') # eventually add gene name and length also with exons
    circles = sorted(results.keys())
    for lola in circles:
	O.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(results[lola]['circle_id'],results[lola]['transcript_ids'], results[lola]['num_reads'], results[lola]['min_length'], results[lola]['max_length'], results[lola]['single'], results[lola]['double'], results[lola]['undefined']))
    O.close()
    return

# run script
tempfile.tempdir = '/beegfs/group_dv/home/FMetge/tmp'

RESULTS = iterate_over_folder(bamfolder, bedfile, platform, split_character)
write_results(RESULTS, outfile)


