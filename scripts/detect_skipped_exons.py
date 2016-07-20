# script to identify skipped exons in circRNA bam files
import pysam
import pybedtools
import os
import argparse
import tempfile

parser = argparse.ArgumentParser(description='Detect genes with different forms of circles')

# input
parser.add_argument('folder', metavar = 'PATH', help = 'PATH to folder containing circle bam files' )
parser.add_argument('bedfile', metavar = 'exons.bed', help = 'exon annotation file' )
parser.add_argument('outfile', metavar = 'outfile', help = 'File to write genes giving rise to alternatively spliced circles to')
parser.add_argument('-p', dest = 'ref_platform', default = 'refseq', help = 'specifies the annotation platform which was used (refseq or ensembl)')
# add temp folder option

args = parser.parse_args()

# parse arguments
folder = args.folder
outfile = args.outfile
bedfile = args.bedfile
platform = args.ref_platform

# define functions


# define functions
def load_bamfile(bamfile):
    reads = {}
    B = pysam.AlignmentFile(bamfile, 'rb')
    for i, lola in enumerate(B):
	if not lola.get_tag('jI') == [-1]:
	    if not lola.query_name in reads:
		reads[lola.query_name] = {}
	    reads[lola.query_name][i] = {'reference': B.getrname(lola.reference_id), 'breakpoint': lola.get_tag('jI'), 'mapq' : lola.mapping_quality}
    B.close()
    return(reads)


def filter_reads(reads):
    for lola in reads:
	if len(reads[lola]) > 1:
	    occurences = reads[lola]
	    KEYS = occurences.keys()
	    mapq = occurences[KEYS[0]]['mapq']
	    for forrest in KEYS:
		if occurences[forrest]['mapq'] > mapq:
		    mapq = occurences[forrest]['mapq']
	    for forrest in KEYS:
		if occurences[forrest]['mapq'] < mapq:
		    del occurences[forrest]
	    reads[lola] = occurences 
    return(reads)



def get_introns(reads):
    introns = {}
    for lola in reads:
	for forrest in reads[lola]:
	    breakpoints = reads[lola][forrest]['breakpoint']
	    starts = breakpoints[::2]
	    ends = breakpoints[1::2]
	    for i,start in enumerate(starts):
		if not (reads[lola][forrest]['reference'], start, ends[i]) in introns:
		    introns[(reads[lola][forrest]['reference'], start, ends[i])] = {'spanning_reads': [],  'skipped_exons': {}}
		introns[(reads[lola][forrest]['reference'], start, ends[i])]['spanning_reads'] += [lola]
    return(introns)

def intersect_introns_with_bedfile(bedfile, reads, introns, coordinates):
    exons = pybedtools.example_bedtool(bedfile)
    exons = exons.filter(lambda b: b.chrom == coordinates[0] and b.start >= coordinates[1] and b.end <= coordinates[2])
    for lola in introns:
	print(lola)
	intron = pybedtools.BedTool('%s %s %s' %(lola[0], lola[1] + 1, lola[2] -1), from_string=True)
	features = exons.intersect(intron)
	for skipped in features:
	    print(skipped)
	    introns[lola]['skipped_exons'][(skipped[0], int(skipped[1]), int(skipped[2]))] = {'name':skipped[3], 'reads':[]}
    return(introns)


def identify_skipped_exons(bamfile, introns):
    bam = pybedtools.example_bedtool(bamfile)
    for lola in introns:
	if len(introns[lola]['skipped_exons']) > 0:
	    for forrest in introns[lola]['skipped_exons']:
		coordinates = pybedtools.BedTool('%s %s %s' %(forrest[0], forrest[1], forrest[2]), from_string=True)
		reads = bam.intersect(coordinates, split = True)
		for r in reads:
		    introns[lola]['skipped_exons'][forrest]['reads'] += [r[0]]
    return(introns)


def write_exon_skipping(introns, outfile, circle_id, platform):
    O = open(outfile, 'a')
    for intron in introns:
	if len(introns[intron]['skipped_exons']) > 0:
	    for exon in introns[intron]['skipped_exons']:
		if len(introns[intron]['skipped_exons'][exon]['reads']) > 0:
		    if platform == 'refseq':
			name = '_'.join(introns[intron]['skipped_exons'][exon]['name'].split('_')[0:2])
		    else:
			name = introns[intron]['skipped_exons'][exon]['name'].split('_')[0]
		    O.write('%s_%s_%s\t%s\t%s:%s-%s\t%s:%s-%s\t%s\t%s\t%s\n' %(circle_id[0], circle_id[1], circle_id[2], name, exon[0], exon[1], exon[2], intron[0], intron[1], intron[2] ,','.join(set(introns[intron]['spanning_reads'])), len(set(introns[intron]['spanning_reads'])), len(set(introns[intron]['skipped_exons'][exon]['reads']))))
    O.close()    
    return

def write_bed12(introns, outfile, circle_id, platform):
    O = open(outfile, 'a')
    for intron in introns:
	if len(introns[intron]['skipped_exons']) > 0:
	    for exon in introns[intron]['skipped_exons']:
		if len(introns[intron]['skipped_exons'][exon]['reads']) > 0 and (circle_id[2] - circle_id[1] > 100) and (circle_id[2]  >= intron[2]) and (circle_id[1] <= intron[1] ):
		    if platform == 'refseq':
			name = '_'.join(introns[intron]['skipped_exons'][exon]['name'].split('_')[0:2])
		    else:
			name = introns[intron]['skipped_exons'][exon]['name'].split('_')[0]
		    if circle_id[0].startswith('chr'):
			O.write('%s\t%s\t%s\t%s\t%s\t.\t%s\t%s\t255,0,0\t3\t1,%s,1\t0,%s,%s\n' %(circle_id[0], circle_id[1], circle_id[2], name, (float(len(set(introns[intron]['spanning_reads'])))/len(set(introns[intron]['skipped_exons'][exon]['reads'])))*100, intron[1], intron[2], (exon[2]-exon[1]), exon[1]-circle_id[1], circle_id[2]-circle_id[1]-1))
		    else:
			O.write('chr%s\t%s\t%s\t%s\t%s\t.\t%s\t%s\t255,0,0\t3\t1,%s,1\t0,%s,%s\n' %(circle_id[0], circle_id[1], circle_id[2], name, (float(len(set(introns[intron]['spanning_reads'])))/len(set(introns[intron]['skipped_exons'][exon]['reads'])))*100, intron[1], intron[2], (exon[2]-exon[1]), exon[1]-circle_id[1], circle_id[2]-circle_id[1]-1))
    O.close()
    return

# run script

tempfile.tempdir = '/beegfs/group_dv/home/FMetge/tmp'

files = os.listdir(folder)
outfile_bed = outfile.replace('.txt', '.bed')


O = open(outfile, 'w')
O.write('circle_id\ttranscript_id\tskipped_exon\tintron\tread_names\tsplice_reads\texon_reads\n')
O.close()

O = open(outfile_bed, 'w')
O.write('# bed12 format\n')
O.close()


for f in files:
    if f.split('.')[-1] == 'bam':
	print(f)
	circle_id = ('_'.join(f.split('_')[0:-3]), int(f.split('_')[-3]), int(f.split('_')[-2]))
	bamfile = '%s/%s' %(folder, f) 
	READS = load_bamfile(bamfile)
	READS = filter_reads(READS)
	Introns = get_introns(READS)
	Introns = intersect_introns_with_bedfile(bedfile, READS, Introns ,circle_id)
	Introns = identify_skipped_exons(bamfile, Introns)
	write_bed12(Introns, outfile_bed, circle_id, platform)
	write_exon_skipping(Introns, outfile, circle_id, platform)

