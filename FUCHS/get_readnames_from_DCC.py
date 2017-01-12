#! /usr/bin/env python


# python script to use CircRNACount file and Chimeric.out.junction
# of mate1/2 to extract read names of circle spanning reads

# define function
def read_circrna_count(infile):
    input_file = open(infile)
    circIDs = []
    for line in input_file:
        current_line = line.replace('\n', '').split('\t')
        circIDs += [(current_line[0], int(current_line[1]), int(current_line[2]))]
    input_file.close()
    return circIDs


def read_junction_file(infile, reads):
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


def filter_circles(circIDs, reads):
    sorted_keys = sorted(reads.keys())
    for key in sorted_keys:
        if not key in circIDs:
            del reads[key]
    return reads


def write_circles(reads, outputfile):
    output_file = open(outputfile, 'w')
    sorted_keys = sorted(reads.keys())
    for key in sorted_keys:
        output_file.write('%s:%s|%s\t%s\n' % (key[0], key[1], key[2], ','.join(set(reads[key]))))
    output_file.close()
    return


# run script

if __name__ == '__main__':

    # required packages
    import argparse

    parser = argparse.ArgumentParser(
        description='Extracts circular reads based on circle_file from the sample.bam '
                    'and writes them into circle separated bam files.')

    # input
    parser.add_argument('circlefile', metavar='sample_circleIDs.txt',
                        help='DCC CircRNACount file containing the filtered list of circles including the coordinates')
    parser.add_argument('junctionfile', metavar='Chimeric.out.junction',
                        help='STAR Chimeric.out.junction file containing the readnames '
                             'and coordinates of chimerically spliced reads.')
    # option
    parser.add_argument('-m', dest='mate2', default='none',
                        help='If data is paired end and mates were mapped separately, '
                             'STAR Chimeric.out.junction file containing the readnames and coordinates '
                             'of chimerically spliced reads of mate2.')

    args = parser.parse_args()

    # parse arguments
    circle_file = args.circlefile
    junctionreads_file = args.junctionfile
    mate2 = args.mate2

    circles = read_circrna_count(circle_file)
    junctions = {}
    junctions = read_junction_file(junctionreads_file, junctions)
    if not mate2 == 'none':
        junctions = read_junction_file(mate2, junctions)

    junctions = filter_circles(circles, junctions)
    write_circles(junctions, '%s.reads.txt' % junctionreads_file)
