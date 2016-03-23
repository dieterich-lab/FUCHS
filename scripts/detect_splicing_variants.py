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
	if not '\n' == line:
	    L = line.replace('\n', '').split('\t')
	    circle_chrom = L[0]
	    circle_start = int(L[1])
	    circle_end = int(L[2])
	    transcript = L[3]
	    num_reads = L[4]
	    strand = L[5]
	    if not transcript in exon_count_table:
		exon_count_table[transcript] = {}
	    exon_count_table[transcript][(circle_chrom, circle_start, circle_end)] = {'score': num_reads, 'strand': strand}
    I.close()
    return(exon_count_table)


def classify_multi_circle_transcripts(transcripts):
    classification = {}
    for lola in transcripts:
	types = {'same_start': {}, 'same_end': {}, 'within': {}, 'overlapping': {}, 'circles':[]}
	circles = sorted(transcripts[lola])
	for i, circle1 in enumerate(circles):
	    types['circles'] += ['%s:%s-%s' %(circle1[0], circle1[1], circle1[2])]
	    for j, circle2 in enumerate(circles):
		if i < j:
		    if circle1[1] == circle2[1]:
			if not circle1[1] in types['same_start']:
			    types['same_start'][circle1[1]] = []
			types['same_start'][circle1[1]] += [('%s:%s-%s' %(circle1[0], circle1[1], circle1[2]), '%s:%s-%s' %(circle2[0], circle2[1], circle2[2]))]
		    elif circle1[2] ==circle2[2]:
			if not circle1[2] in types['same_end']:
			    types['same_end'][circle1[2]] = []
			types['same_end'][circle1[2]] += [('%s:%s-%s' %(circle1[0], circle1[1], circle1[2]), '%s:%s-%s' %(circle2[0], circle2[1], circle2[2]))]
		    elif (circle1[1] < circle2[1] and circle1[2] > circle2[2]) or (circle1[1] > circle2[1] and circle1[2] < circle2[2]):
			if not i in types['within']:
			    types['within'][i] = []
			types['within'][i] += [('%s:%s-%s' %(circle1[0], circle1[1], circle1[2]), '%s:%s-%s' %(circle2[0], circle2[1], circle2[2]))]
		    elif (circle1[1] < circle2[1] and circle1[2] <circle2[2] and circle1[2] > circle2[1]) or (circle1[1] > circle2[1] and circle1[2] > circle2[2] and circle1[1] < circle2[2]):
			if not i in types['overlapping']:
			    types['overlapping'][i] = []
			types['overlapping'][i] += [('%s:%s-%s' %(circle1[0], circle1[1], circle1[2]), '%s:%s-%s' %(circle2[0], circle2[1], circle2[2]))]
	classification[lola] = types
    return(classification)


def write_genes(types, circles, outfile):
    O = open(outfile, 'w')
    O.write('Transcript\tcircles\tsame_start\tsame_end\toverlapping\twithin\n')
    for lola in types:
	O.write('%s\t%s' %(lola, ','.join(types[lola]['circles'])))
	if len(types[lola]['same_start']) == 0:
	    O.write('\t.')
	elif len(types[lola]['same_start']) > 0:
	    O.write('\t')
	    for circ in types[lola]['same_start']:
		O.write('%s,' %('|'.join(types[lola]['same_start'][circ][0])))
	if len(types[lola]['same_end']) == 0:
	    O.write('\t.')
	elif len(types[lola]['same_end']) > 0:
	    O.write('\t')
	    for circ in types[lola]['same_end']:
		O.write('%s,' %('|'.join(types[lola]['same_end'][circ][0])))
	if len(types[lola]['overlapping']) == 0:
	    O.write('\t.')
	elif len(types[lola]['overlapping']) > 0:
	    O.write('\t')
	    for circ in types[lola]['overlapping']:
		O.write('%s,' %('|'.join(types[lola]['overlapping'][circ][0])))
	if len(types[lola]['within']) == 0:
	    O.write('\t.')
	elif len(types[lola]['within']) > 0:
	    O.write('\t')
	    for circ in types[lola]['within']:
		O.write('%s,' %('|'.join(types[lola]['within'][circ][0])))
	O.write('\n')
    O.close()
    return

#exon_count_table_file = '/beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/FUCHS/old_cerebellum.exon_counts.genes.txt'
#exon_count_table_file = 'hek293_A.exon_counts.bed'
#outfile = 'test.txt'
table = read_exon_count_table(exon_count_table_file)
classification = classify_multi_circle_transcripts(table)
write_genes(classification, table, outfile)

