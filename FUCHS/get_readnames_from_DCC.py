#! /usr/bin/env python


# python script to use CircRNACount file and Chimeric.out.junction
# of mate1/2 to extract read names of circle spanning reads

class get_readnames_from_DCC(object):
    def __init__(self, mate2, circlefile, junctionfile, ):

        self.circle_file = circlefile
        self.junctionreads_file = junctionfile
        self.mate2 = mate2

    # define functions
    def read_circrna_count(self, infile):
        input_file = open(infile)
        circIDs = []
        for line in input_file:
	    if not line.startswith('#') and not line.startswith('Chr\t'):
		current_line = line.replace('\n', '').split('\t')
		circIDs += [(current_line[0], int(current_line[1]), int(current_line[2]))]
        input_file.close()
        return circIDs

    def read_junction_file(self, infile, reads):
        input_file = open(infile)
        for line in input_file:
            current_line = line.replace('\n', '').split('\t')
            chrom = current_line[0]
            coord = (int(current_line[1]), int(current_line[4]))
            start = min(coord)
            end = max(coord)
            if not (chrom, start + 1, end - 1) in reads:
                reads[(chrom, start + 1, end - 1)] = [current_line[9]]
            else:
                reads[(chrom, start + 1, end - 1)] += [current_line[9]]
        input_file.close()
        return reads

    def filter_circles(self, circIDs, reads):
        sorted_keys = sorted(reads.keys())
        for key in sorted_keys:
            if not key in circIDs:
                del reads[key]
        return reads

    def write_circles(self, reads, outputfile):
        output_file = open(outputfile, 'w')
        sorted_keys = sorted(reads.keys())
        for key in sorted_keys:
            output_file.write('%s:%s|%s\t%s\n' % (key[0], key[1], key[2], ','.join(set(reads[key]))))
        output_file.close()
        return

    # run script
    def run(self):

        circles = self.read_circrna_count(self.circle_file)
        junctions = {}
        junctions = self.read_junction_file(self.junctionreads_file, junctions)
        if not self.mate2 == 'none':
            junctions = self.read_junction_file(self.mate2, junctions)

        junctions = self.filter_circles(circles, junctions)
        self.write_circles(junctions, '%s.reads.txt' % self.junctionreads_file)
