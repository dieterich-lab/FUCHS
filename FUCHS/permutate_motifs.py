#! /usr/bin/env python3

# script to permute motif files


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
                motif = motif.replace(' \n', '_%s \n' % (counter))
                counter += 1
        elif 'letter-probability' in Line:
            read_motif = True
            motifs[motif] = {}
            motifs[motif]['real_motif'] = [Line]
        elif Line == '\n':
            read_motif = False
        elif Line.count('\t') == 4 and read_motif:
            motifs[motif]['real_motif'] += [Line]
    I.close()
    return (motifs, header)


def permutate_motifs(motifs):
    """
    generates a pool of unique permutations
    """
    for m in motifs:
        permutations_temp = itertools.permutations(motifs[m]['real_motif'][1:])
        permutations = []
        for p in permutations_temp:
            permutations += [p]
        motifs[m]['permutations'] = list(set(permutations))
    return (motifs)


def write_permutated_motifs(outfile, motifs, header):
    O = open(outfile, 'w')
    O.write(header)
    for m in motifs:
        if len(motifs[m]['permutations']) > 0:  # if there are permutatations left, write a random permutations
            tmp = motifs[m]['permutations'].pop(random.randint(0, len(motifs[m]['permutations']) - 1))
            if not tmp == motifs[m]['real_motif'][
                          1:]:  # make sure the permutated motif is not the same as the real motif
                O.write('%s\n\n%s%s\n' % (m.split('\n')[0], motifs[m]['real_motif'][0], ''.join(tmp)))
            elif len(motifs[m]['permutations']) > 0:
                tmp = motifs[m]['permutations'].pop(random.randint(0, len(motifs[m]['permutations'])) - 1)
                O.write('%s\n\n%s%s\n' % (m.split('\n')[0], motifs[m]['real_motif'][0], ''.join(tmp)))
            else:  # if there are no more permutations to choose from, write original motif
                O.write('%s\n\n%s\n' % (m.split('\n')[0], ''.join(motifs[m]['real_motif'])))
        else:  # if there are no more permutations to choose from, write original motif
            O.write('%s\n\n%s\n' % (m.split('\n')[0], ''.join(motifs[m]['real_motif'])))
    O.close()
    return


# run script
if __name__ == '__main__':

    # required packages
    import random
    import itertools
    import argparse

    # input
    parser = argparse.ArgumentParser(description='Takes a meme formatted motif file and permutates each motif n times.')

    parser.add_argument('meme_file', metavar='motifs.meme',
                        help='Meme formatted files with motifs to permutate for p-value calculations of motif enrichment.')
    parser.add_argument('output_prefix', metavar='prefix_for_path',
                        help='prefix for for all permuation files to. Permuations will be written to PREFIX.permuation.i.meme')
    # options
    parser.add_argument('-n', dest='repeats', default=100, help='number of permutations to perform')

    # parse arguments
    args = parser.parse_args()

    meme_file = args.meme_file
    output_prefix = args.output_prefix
    num_perm = args.repeats

    # run
    MOTIFS, HEADER = read_in_motif_file(meme_file)
    PERMS = permutate_motifs(MOTIFS)
    O = open('%s.number_of_unique_motifs.txt' % (output_prefix), 'w')
    O.write('motif\tpermutation\n')
    for m in PERMS:
        O.write('%s\t%s\n' % (m.split('\n')[0], len(PERMS[m]['permutations'])))

    O.close()
    for i in range(0, num_perm):
        write_permutated_motifs('%s.permutation.%s.meme' % (output_prefix, i), MOTIFS, HEADER)
