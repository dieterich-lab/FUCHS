# script to get background sequence 

import pysam


def read_bed_files(infile):
    exons = {}
    I = open(infile)
    for line in I:
	L = line.replace('\n', '').split('\t')
	if not (L[0], int(L[1]), int(L[2])) in exons:
	    exons[(L[0], int(L[1]), int(L[2]))] = {'transcripts' : [], 'strand' : L[5], 'sequence': ''}
	exons[(L[0], int(L[1]), int(L[2]))]['transcripts'] += [L[3].split('_')[0]]
    I.close()
    return(exons)


def fetch_fasta(exons, fastafile):
    ref = pysam.Fastafile(fastafile)
    for e in exons:
	exons[e]['sequence'] = ref.fetch(e[0], e[1], e[2]).upper()
    return(exons)

def reconstruct_transcripts(exons):
    transcripts = {}
    for e in exons:
	for t in exons[e]['transcripts']:
	    if not t in transcripts:
		transcripts[t] = {'exons': {}, 'strand' : exons[e]['strand']}
	    transcripts[t]['exons'][e] = exons[e]['sequence']
    return(transcripts)

def paste_sequences(transcripts):
    fasta = {}
    for t in transcripts:
	transcripts_fasta = ''
	sorted_exons = sorted(transcripts[t]['exons'].keys())
	for e in sorted_exons:
	    transcripts_fasta += transcripts[t]['exons'][e]
	fasta[t] = {'sequence': transcripts_fasta, 'strand': transcripts[t]['strand']}
    return(fasta)

def reverse_complement(sequence):
    letter_dict = {'A':'T', 'C': 'G', 'G':'C', 'T':'A', 'N':'N'}
    revcomp = ''
    for lola in reversed(sequence):
	revcomp += letter_dict[lola]
    return(revcomp)

def write_fasta(outfile, fasta):
    O = open(outfile, 'w')
    sorted_transcripts = sorted(fasta.keys())
    for t in sorted_transcripts:
	if fasta[t]['strand'] == '+':
		O.write('>%s_+\n%s\n' %(t, fasta[t]['sequence']))
	else:
		O.write('>%s_-\n%s\n' %(t, reverse_complement(fasta[t]['sequence'])))
    O.close()
    return
  

bedfile = 'mm10.ensembl.exons.bed'
fastafile = 'GRCm38.dna.toplevel.fa'
outfile = 'mm10.exon_sequences.fa'


Exons = read_bed_files(bedfile)
Exons = fetch_fasta(Exons, fastafile)
Transcripts = reconstruct_transcripts(Exons)
Fasta = paste_sequences(Transcripts)
write_fasta(outfile, Fasta)







