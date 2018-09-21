#! /usr/bin/env python3


# script to replace transcript ids with gene names:

# define functions
def read_names_file(infile):
    id_dict = {}
    I = open(infile)
    for lola in I:
        id_dict[lola.split('\t')[0]] = lola.replace('\n', '').split('\t')[1]
    I.close()
    return (id_dict)


def read_mate_status(infile):
    mates = {}
    I = open(infile)
    header = I.readline()
    for lola in I:
        L = lola.replace('\n', '').split('\t')
        mates[L[0]] = {'line': L, 'transcripts': L[1].split(',')}
    I.close()
    return (mates, header)


def read_exon_counts(infile):
    exons = {}
    I = open(infile)
    header = I.readline()
    for i, lola in enumerate(I):
        L = lola.replace('\n', '').split('\t')
        exons[i] = {'line': L, 'transcripts': [L[2]]}
    I.close()
    return (exons, header)


def replace_names(id_dict, table):
    for lola in table:
        gene_names = []
        for transcript in table[lola]['transcripts']:
            if transcript in id_dict:
                gene_names += [id_dict[transcript]]
            else:
                gene_names += [transcript]
            table[lola]['gene_name'] = ','.join(set(gene_names))
    return (table)


def write_table(table, header, index, outfile):
    O = open(outfile, 'w')
    O.write(header)
    KEYS = sorted(table.keys())
    for lola in KEYS:
        table[lola]['line'][index] = table[lola]['gene_name']
        O.write('\t'.join(table[lola]['line']))
        O.write('\n')
    O.close()
    return


# run script

if __name__ == '__main__':

    # required packages
    import argparse

    parser = argparse.ArgumentParser(description='Replace transcript IDs with gene names.')

    # input
    parser.add_argument('transcript_dictionary', metavar='gene_names.txt',
                        help='tab separated file. Transcript ID in column 1, gene name in column 2.')
    # options
    parser.add_argument('-e', dest='exons', default='',
                        help='FUCHS generated exon_count table to replace transcript IDs with gene names')
    parser.add_argument('-m', dest='mates', default='',
                        help='FUCHS generated mate_status table to replace transcript IDs with gene names.')

    args = parser.parse_args()

    # parse arguments
    id_dict = args.transcript_dictionary
    exons = args.exons
    mates = args.mates

    IDs = read_names_file(id_dict)

    if not exons == '':
        E, H = read_exon_counts(exons)
        E = replace_names(IDs, E)
        write_table(E, H, 2, exons.replace('.txt', '.genes.txt'))

    if not mates == '':
        M, H = read_mate_status(mates)
        M = replace_names(IDs, M)
        write_table(M, H, 1, mates.replace('.txt', '.genes.txt'))
