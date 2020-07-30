#! /usr/bin/env python3


# calculate the proportion of circRNA covered by reads. including introns. should probably be independent of mate_status.txt file.


# define functions
def read_mate_status(infile):
    I = open(infile)
    circles = {}
    H = I.readline()
    for line in I:
        L = line.replace('\n', '').split('\t')
        circles[L[0]] = L[1:]
    I.close()
    return (circles, H)


def read_coverage_profile(infile):
    I = open(infile)
    coverage = []
    H = I.readline()
    for line in I:
        coverage += [int(line.replace('\n', '').split('\t')[3])]
    I.close()
    return (coverage)


def iterate_over_circRNAs(mates, coverage_folder):
    coverage_profiles = os.listdir(coverage_folder)
    for lola in mates:
        if not mates[lola][0] == 'not_annotated':
            if mates[lola][0].split(',') > 1:
                T = mates[lola][0].split(',')
                for transcript in T:
                    F = '%s.%s.txt' % (lola, transcript)
                    if F in coverage_profiles:
                        C = read_coverage_profile('%s/%s' % (coverage_folder, F))
            else:
                F = '%s/%s.%s.txt' % (coverage_folder, lola, mates[lola][0])
                C = read_coverage_profile(F)
            mates[lola] += ['%s' % (len(C)), '%s' % (C.count(0)), '%s' % (1 - (float(C.count(0)) / len(C)))]
        else:
            mates[lola] += ['NA', 'NA', 'NA']
    return (mates)


def write_mates(mates, outfile, header):
    O = open(outfile, 'w')
    SK = sorted(mates.keys())
    O.write('%s\tannotated_length\tbases0\tproportion_covered_bases\n' % (header.replace('\n', '')))
    for lola in SK:
        O.write('%s\t%s\n' % (lola, '\t'.join(mates[lola])))
    O.close()
    return


# run script

if __name__ == '__main__':
    # required packages
    import os
    import argparse

    parser = argparse.ArgumentParser(
        description='calculates the proportion of circRNA that was covered by reads, including introns')

    # input
    parser.add_argument('mate_file', metavar='sample.mate_status.txt',
                        help='file listing the number of single and double breakpoint fragments per circle.')
    parser.add_argument('folder', metavar='sample.coverage_profiles',
                        help='Folder containing the coverage profiles of all circles.')

    args = parser.parse_args()

    # parse_arguments
    mate_file = args.mate_file
    folder = args.folder

    # run
    M, H = read_mate_status(mate_file)
    M = iterate_over_circRNAs(M, folder)

    write_mates(M, mate_file.replace('.txt', '.proportion_covered.txt'), H)
