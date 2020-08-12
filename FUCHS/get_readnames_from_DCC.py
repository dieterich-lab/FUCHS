from functools import reduce

# python script to use CircRNACount file and Chimeric.out.junction
# of mate1/2 to extract read names of circle spanning reads

class get_readnames_from_DCC(object):
    def __init__(self, circlefile, junction_file, mate1, mate2, ):

        self.circle_file = circlefile
        self.junction_file = junction_file
        self.mate1 = mate1
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
                reads[(chrom, start + 1, end - 1)] = {'paired': [], 'mate1': [], 'mate2': []}
            else:
                reads[(chrom, start + 1, end - 1)]['paired'] += [current_line[9]]
        input_file.close()
        return reads

    def read_mate_junction_file(self, infile, reads, mate):
        input_file = open(infile)
        for line in input_file:
            current_line = line.replace('\n', '').split('\t')
            chrom = current_line[0]
            coord = (int(current_line[1]), int(current_line[4]))
            start = min(coord)
            end = max(coord)
            if not (chrom, start + 1, end - 1) in reads:
                reads[(chrom, start + 1, end - 1)] = {'paired': [], 'mate1': [], 'mate2': []}
            else:
                reads[(chrom, start + 1, end - 1)][mate] += [current_line[9]]
        input_file.close()
        return reads

    def filter_circles_by_circID(self, circIDs, reads):
        sorted_keys = sorted(reads.keys())
        for key in sorted_keys:
            if not key in circIDs:
                del reads[key]
        return reads

    def filter_reads_by_mate(self, reads, is_paired):
        unique_reads = {}
        for circ in list(reads.keys()):
            unique_reads[circ] = []
            all_reads = list(set(reduce(lambda x, y: x + y, list(reads[circ].values()), [])))
            if is_paired:
                for read in all_reads:
                    if read in reads[circ]['paired'] and read in reads[circ]['mate1'] and read in reads[circ]['mate2']:
                        print(('false positive read %s in %s' % (read, circ)))
                    else:
                        unique_reads[circ] += [read]
            else:
                unique_reads[circ] += all_reads
        return (unique_reads)

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
        junctions = self.read_junction_file(self.junction_file, junctions)
        if not self.mate1 == 'none':
            junctions = self.read_mate_junction_file(self.mate1, junctions, 'mate1')
        if not self.mate2 == 'none':
            junctions = self.read_mate_junction_file(self.mate2, junctions, 'mate2')
        junctions = self.filter_circles_by_circID(circles, junctions)
        if not self.mate1 == 'none' and not self.mate2 == 'none':
            unique_reads = self.filter_reads_by_mate(junctions, True)
        else:
            unique_reads = self.filter_reads_by_mate(junctions, False)
        self.write_circles(unique_reads, '%s.reads.txt' % self.junction_file)
