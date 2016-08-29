# script to permutate motif files

# define functions
def read_in_motif_file(infile):
    motifs = {}
    I = open(infile)
    header = ''
    for i in range(0, 9):
	header += I.readline()
    motif = ''
    counter = 0
    for Line in I:
	if 'MOTIF' in Line:
	    motif = Line
	    if motif in motifs:
		motif = motif + str(counter)
		counter += 1
	elif 'letter-probability' in Line:
	    read_motif = True
	    motifs[motif] = [Line]
	elif Line == '\n':
	    read_motif = False
	elif Line.count('\t') == 4 and read_motif:
	    motifs[motif] += [Line]
    I.close()
    return(motifs, header)

def write_permutated_motifs(outfile, motifs, header):
    O = open(outfile, 'w')
    O.write(header)
    for m in motifs:
	tmp = motifs[m][1:]
	random.shuffle(tmp)
	O.write('%s\n\n%s%s\n' %(m.split('\n')[0], motifs[m][0] ,''.join(tmp)))
    O.close()
    return


# run script
if __name__ == '__main__':
    
    # required packages
    import random
    import argparse
    
    # input
    parser = argparse.ArgumentParser(description='Takes a meme formatted motif file and permutates each motif n times.')
    
    parser.add_argument('meme_file', metavar = 'motifs.meme', help = 'Meme formatted files with motifs to permutate for p-value calculations of motif enrichment.' )
    parser.add_argument('output_prefix', metavar = 'prefix_for_path', help = 'prefix for for all permuation files to. Permuations will be written to PREFIX.permuation.i.meme' )
    # options
    parser.add_argument('-n', dest = 'repeats', default = 100, help = 'number of permutations to perform')
    
    
    # parse arguments
    args = parser.parse_args()
    
    meme_file = args.meme_file
    output_prefix = args.output_prefix
    permuations = args.repeats
    
    # run
    MOTIFS, HEADER = read_in_motif_file(meme_file)
    for i in range(0, permuations):
	write_permutated_motifs('%s.permutation.%s.meme' %(output_prefix, i), MOTIFS, HEADER)


