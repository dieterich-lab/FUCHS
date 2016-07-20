# file to extract and paste together the exon sequence of a transcripts

import pysam
import argparse
import tempfile

parser = argparse.ArgumentParser(description='')
# input
parser.add_argument('fasta', metavar = 'fasta.fa', help = 'reference fasta file')
parser.add_argument('infile', metavar = 'infile', help = 'bedfile containing exon coordinates')
parser.add_argument('outfile', metavar = 'outfile.fa', help = 'file, sequences should be written to')
parser.add_argument('-A', dest = 'CircCoordinates', default = '', help = 'give a circle Annotation file to add strand and gene name to circles, right now this should be CircCoordinates from DCC')


# options

args = parser.parse_args()

# parse arguments

infile = args.infile
fastafile = args.fasta
outfile = args.outfile
circRNAfile = args.CircCoordinates

def read_bedfile(infile):
    I = open(infile)
    exons = {}
    for line in I:
	if not line.startswith('#'):
	    L = line.replace('\n', '').split('\t')
	    circID = L[3].split('|')[0]
	    transcripts = L[3].split('|')[1].split(',')
	    if not (L[0], int(L[1]), int(L[2])) in exons:
		exons[(L[0], int(L[1]), int(L[2]))] = {'circles' : {}, 'sequence': ''}
	    exons[(L[0], int(L[1]), int(L[2]))]['circles'][circID] = {'transcripts': transcripts,  'score': int(L[4]), 'strand': L[5]}    
    return(exons)

def read_circle_annotation(infile):
    I = open(infile)
    I.readline()
    circle_annotation = {}
    for line in I:
	L = line.replace('\n', '').split('\t')
	circle_annotation['%s:%s-%s' %(L[0], L[1], L[2])] = {'gene' : L[3], 'strand': L[5]}
    I.close()
    return(circle_annotation)

def fetch_fasta(exons, fastafile):
    ref = pysam.Fastafile(fastafile)
    for e in exons:
	exons[e]['sequence'] = ref.fetch(e[0], e[1], e[2]).upper()
    return(exons)

def reconstruct_transcripts(exons):
    transcripts = {}
    for e in exons:
	for c in exons[e]['circles']:
	    if not c in transcripts:
		transcripts[c] = {}
	    for t in exons[e]['circles'][c]['transcripts']:
		if not t in transcripts[c]:
		    transcripts[c][t] = {}
		transcripts[c][t][e] = {'sequence': exons[e]['sequence'] , 'score': exons[e]['circles'][c]['score']}
    return(transcripts)

def paste_sequences(transcripts):
    fasta = {}
    for c in transcripts:
	for t in transcripts[c]:
	    transcripts_fasta = ''
	    average_score = 0
	    sorted_exons = sorted(transcripts[c][t].keys())
	    for e in sorted_exons:
		transcripts_fasta += transcripts[c][t][e]['sequence']
		average_score += transcripts[c][t][e]['score']
	    fasta[(c,t)] = {'sequence': transcripts_fasta, 'score': average_score/len(sorted_exons)}
    return(fasta)

def reverse_complement(sequence):
    letter_dict = {'A':'T', 'C': 'G', 'G':'C', 'T':'A', 'N':'N'}
    revcomp = ''
    for lola in reversed(sequence):
	revcomp += letter_dict[lola]
    return(revcomp)

def write_fasta(outfile, fasta, circAnnotation):
    O = open(outfile, 'w')
    sorted_circles = sorted(fasta.keys())
    for c in sorted_circles:
	if c[0] in circAnnotation:
	    if circAnnotation[c[0]]['strand'] == '+':
		O.write('>%s_%s_%s_score:%s_+\n%s\n' %(c[0], c[1], circAnnotation[c[0]]['gene'], fasta[c]['score'], fasta[c]['sequence']))
	    else:
		O.write('>%s_%s_%s_score:%s_-\n%s\n' %(c[0], c[1], circAnnotation[c[0]]['gene'], fasta[c]['score'], reverse_complement(fasta[c]['sequence'])))
	else:
	    O.write('>%s_%s_NA_score:%s_+\n%s\n' %(c[0], c[1], fasta[c]['score'], fasta[c]['sequence']))
	    O.write('>%s_%s_NA_score:%s_-\n%s\n' %(c[0], c[1], fasta[c]['score'], reverse_complement(fasta[c]['sequence'])))
    O.close()
    return
  


if not circRNAfile == '':
    circAnnotation = read_circle_annotation(circRNAfile)
else:
    circAnnotation = {}

Exons = read_bedfile(infile)
Exons = fetch_fasta(Exons, fastafile)
Transcripts = reconstruct_transcripts(Exons)
Fasta = paste_sequences(Transcripts)
write_fasta(outfile, Fasta, circAnnotation)



