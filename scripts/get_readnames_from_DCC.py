# python script to use CircRNACount file and Chimeric.out.junction of mate1/2 to extract readnames of circle spanning reads

# define function
def read_circRNA_count(infile):
    I = open(infile)
    circIDs = []
    header = I.readline()
    for line in I:
	  L = line.replace('\n', '').split('\t')
	  circIDs += [(L[0], int(L[1]), int(L[2]))]    
    I.close()
    return(circIDs)

def read_junction_file(infile, reads):
    I = open(infile)
    for line in I:
	L = line.replace('\n', '').split('\t')
	chrom = L[0]
	coord = (int(L[1]), int(L[4]))
	start = min(coord)
	end = max(coord)
	if not (chrom, start +1 , end - 1) in reads:
	    reads[(chrom, start +1, end-1)] = [L[9]]
	else:
	    reads[(chrom, start+1, end-1)] += [L[9]]
    I.close()
    return(reads)

def filter_circles(circIDs, reads):
    SK = sorted(reads.keys())
    for lola in SK:
	if not lola in circIDs:
	    del reads[lola]
    return(reads)

def write_circles(reads, outputfile):
    O = open(outputfile, 'w')
    SK = sorted(reads.keys())
    for lola in SK:
	O.write('%s:%s|%s\t%s\n' %(lola[0], lola[1], lola[2], ','.join(set(reads[lola]))))
    O.close()
    return

# run script

if __name__ == '__main__':
    
    # required packages
    import argparse
    
    parser = argparse.ArgumentParser(description='Extracts circular reads based on circle_file from the sample.bam and writes them into circle separated bamfiles.')
    
    # input
    parser.add_argument('circlefile', metavar = 'sample_circleIDs.txt', help = 'DCC CircRNACount file containing the filtered list of circles including the coordinates' )
    parser.add_argument('junctionfile', metavar = 'Chimeric.out.junction', help = 'STAR Chimeric.out.junction file containing the readnames and coordinates of chimerically spliced reads.')
    # option
    parser.add_argument('-m', dest = 'mate2', default = 'none', help = 'If data is paired end and mates were mapped separately, STAR Chimeric.out.junction file containing the readnames and coordinates of chimerically spliced reads of mate2.')
    
    args = parser.parse_args()

    # parse arguments
    circle_file = args.circlefile
    junctionreads_file = args.junctionfile
    mate2 = args.mate2
    
    
    C = read_circRNA_count(circle_file)
    R = {}
    R = read_junction_file(junctionreads_file, R)
    if not mate2 == 'none':
	R = read_junction_file(mate2, R)
    
    R = filter_circles(C,R)
    write_circles(R, '%s.reads.txt' %(junctionreads_file))
