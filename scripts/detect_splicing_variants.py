# script to evaluate circle splicing variants


import argparse

parser = argparse.ArgumentParser(description='Detect genes with different forms of circles')

# input
parser.add_argument('exonfile', metavar = 'exon_count.genes.txt', help = 'Table of exon_counts per circle and transcript information' )
parser.add_argument('outfile', metavar = 'outfile', help = 'File to write genes giving rise to different circles to' )

args = parser.parse_args()

# parse arguments
exon_count_table_file = args.exonfile
outfile = args.outfile


def read_exon_count_table(infile):
    I = open(infile)
    exon_count_table = {}
    header = I.readline()
    for line in I:
	L = line.replace('\n', '').split('\t')
	sample = L[0]
	circle_id = L[1]
	transcript = L[2]
	alternatives = L[3].split(',')
	exon_id = int(L[4])
	exon_chrom = L[5]
	exon_start = int(L[6])
	exon_end = int(L[7])
	strand = L[8]
	exon_length = int(L[9])
	unique_reads = int(L[10])
	fragmented_reads = int(L[11])
	forward_reads = int(L[12])
	reverse_reads = int(L[13])
	if not (circle_id, transcript) in exon_count_table:
	    exon_count_table[(circle_id, transcript)] = {}
	exon_count_table[(circle_id, transcript)][exon_id] = {'alternatives' : alternatives, 'exon_chrom': exon_chrom, 'exon_start' : exon_start, 'strand' : strand, 'unique_reads': unique_reads, 'fragmented_reads': fragmented_reads}
    I.close()
    return(exon_count_table)


def detect_circle_variants_in_transcript(exon_table):
    transcripts = {}
    for lola in exon_table:
	if lola[1] not in transcripts:
	    transcripts[lola[1]] = [lola[0]]
	else:
	    transcripts[lola[1]] += [lola[0]]
    all_transcripts = transcripts.keys()
    for lola in all_transcripts:
	if len(transcripts[lola]) == 1:
	    del transcripts[lola]
    return(transcripts)

def classify_multi_circle_transcripts(transcripts):
    classification = {}
    for lola in transcripts:
	circles = []
	types = {'same_start': {}, 'same_end': {}, 'within': {}, 'overlapping': {}}
	for forrest in transcripts[lola]:
	    circles += [(forrest.split('_'))]
	for count,i in enumerate(circles):
	    for j in circles:
		if i < j:
		    if i[1] == j[1]:
			if not i[1] in types['same_start']:
			    types['same_start'][i[1]] = []
			types['same_start'][i[1]] += [('_'.join(i), '_'.join(j))]
		    elif i[2] == j[2]:
			if not i[2] in types['same_end']:
			    types['same_end'][i[2]] = []
			types['same_end'][i[2]] += [('_'.join(i), '_'.join(j))]
		    elif (i[1] < j[1] and i[2] > j[2]) or (i[1] > j[1] and i[2] < j[2]):
			if not count in types['within']:
			    types['within'][count] = []
			types['within'][count] += [('_'.join(i), '_'.join(j))]
		    elif (i[1] < j[1] and i[2] < j[2] and i[2] > j[1]) or (i[1] > j[1] and i[2] > j[2] and i[1] < j[2]):
			if not count in types['overlapping']:
			    types['overlapping'][count] = []
			types['overlapping'][count] += [('_'.join(i), '_'.join(j))]
	classification[lola] = types
    return(classification)

def write_genes(types, circles, outfile):
    O = open(outfile, 'w')
    O.write('Transcript\tcircles\tsame_start\tsame_end\toverlapping\twithin\n')
    for lola in types:
	O.write('%s\t%s' %(lola, ','.join(circles[lola])))
	if len(types[lola]['same_start']) == 0:
	    O.write('\t.')
	elif len(types[lola]['same_start']) > 0:
	    O.write('\t')
	    for circ in types[lola]['same_start']:
		O.write('%s,' %(':'.join(types[lola]['same_start'][circ][0])))
	if len(types[lola]['same_end']) == 0:
	    O.write('\t.')
	elif len(types[lola]['same_end']) > 0:
	    O.write('\t')
	    for circ in types[lola]['same_end']:
		O.write('%s,' %(':'.join(types[lola]['same_end'][circ][0])))
	if len(types[lola]['overlapping']) == 0:
	    O.write('\t.')
	elif len(types[lola]['overlapping']) > 0:
	    O.write('\t')
	    for circ in types[lola]['overlapping']:
		O.write('%s,' %(':'.join(types[lola]['overlapping'][circ][0])))
	if len(types[lola]['within']) == 0:
	    O.write('\t.')
	elif len(types[lola]['within']) > 0:
	    O.write('\t')
	    for circ in types[lola]['within']:
		O.write('%s,' %(':'.join(types[lola]['within'][circ][0])))
	O.write('\n')
    O.close()
    return

#exon_count_table_file = '/beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/FUCHS/old_cerebellum.exon_counts.genes.txt'

table = read_exon_count_table(exon_count_table_file)
multi_circle_transcripts = detect_circle_variants_in_transcript(table)
classification = classify_multi_circle_transcripts(multi_circle_transcripts)
write_genes(classification, multi_circle_transcripts, outfile)
