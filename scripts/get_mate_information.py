# script to identify circles where both mates map over the junction

import pysam
import os
import argparse

parser = argparse.ArgumentParser(description='Extracts mate information and identify singe and double breakpoint fragments')

# input
parser.add_argument('bamfolder', metavar = 'PATH', help = 'path to folder containing circle bamfiles' )
parser.add_argument('outfile', metavar = 'outfile', help = 'path and filename to write the output to' )

args = parser.parse_args()

# parse arguments
bamfolder = args.bamfolder
outfile = args.outfile

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

def iterate_over_folder(inputfolder):
    results = {}
    files = os.listdir(inputfolder)
    for lola in files:
	if lola.split('.')[-1] == 'bam':
	    print(lola)
	    circle_coordinates = [lola.split('_')[0], int(lola.split('_')[1]), int(lola.split('_')[2])]
	    num_reads = lola.split('.')[0].split('_')[3].replace('reads', '')
	    print(circle_coordinates)
	    MATES, FRAGMENTS = get_reads_from_bamfile('%s/%s' %(inputfolder, lola), circle_coordinates)
	    MATES = classify_reads(MATES)
	    STATS = get_statistics(MATES)
	    results[lola.split('.')[0]] = STATS  
	    results[lola.split('.')[0]]['circle_id'] = '%s_%s_%s' %(circle_coordinates[0], circle_coordinates[1], circle_coordinates[2])
	    results[lola.split('.')[0]]['num_reads'] = num_reads
    return(results)


def write_results(results, outfile):
    O = open(outfile, 'w')
    O.write('circle_id\tnum_reads\tsingle\tdouble\tundefined\n') # eventually add gene name and length also with exons
    circles = sorted(results.keys())
    for lola in circles:
	O.write('%s\t%s\t%s\t%s\t%s\n' %(results[lola]['circle_id'], results[lola]['num_reads'], results[lola]['single'], results[lola]['double'], results[lola]['undefined']))
    O.close()
    return

# run script

RESULTS = iterate_over_folder(bamfolder)
write_results(RESULTS, outfile)
