# script to evaluate circle splicing variants

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
	
    return(classification)

exon_count_table_file = '/beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/FUCHS/test.exon_counts.txt'

table = read_exon_count_table(exon_count_table_file)
multi_circle_transcripts = detect_circle_variants_in_transcript(table)

