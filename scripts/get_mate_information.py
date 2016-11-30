#! /usr/bin/env python


# script to identify circles where both mates map over the junction
# length estimate is not correct for circRNAs exceeding the host-gene coordinates. need to fix this!

# define functions
def get_reads_from_bamfile(bamfile, circle_coordinates):
    mates = {}
    non_junction_fragments = []
    circles = pysam.AlignmentFile(bamfile, "rb")
    for circle in circles:
        name = circle.query_name
        reverse = circle.is_reverse
        start = circle.reference_start
        end = circle.reference_end
        if not name in mates:
            mates[name] = {'forward': {'start': [], 'end': []}, 'reverse': {'start': [], 'end': []}}
        if reverse and end == circle_coordinates[2]:
            mates[name]['reverse']['end'] += [start]
        elif reverse and start == circle_coordinates[1] - 1:
            mates[name]['reverse']['start'] += [end]
        elif end == circle_coordinates[2] and not reverse:
            mates[name]['forward']['end'] += [start]
        elif start == circle_coordinates[1] - 1 and not reverse:
            mates[name]['forward']['start'] += [end]
        else:
            non_junction_fragments += [circle]
    circles.close()
    return mates, non_junction_fragments


def classify_reads(mates):
    for mate in mates:
        strands = 0
        for strand in mates[mate]:
            if len(mates[mate][strand]['start']) == 1 and len(mates[mate][strand]['end']) == 1:
                strands += 1
        if strands == 1:
            mates[mate]['status'] = 'single'
        elif strands == 2:
            mates[mate]['status'] = 'double'
        else:
            mates[mate]['status'] = 'undefined'
    return mates


def get_statistics(mates):
    stats = {'single': 0, 'double': 0, 'undefined': 0}
    for mate in mates:
        stats[mates[mate]['status']] += 1
    return stats


def annotate_circle(circle_coordinates, bedfile, platform, split_character):
    circle = pybedtools.BedTool('%s %s %s' % (circle_coordinates[0], circle_coordinates[1], circle_coordinates[2]),
                                from_string=True)
    exons = pybedtools.example_bedtool(bedfile)
    features = exons.intersect(circle)
    lengths = {}
    for feature in features:
        if platform == 'refseq':
            transcript_name = split_character.join(feature[3].split(split_character)[0:2])
        elif platform == 'ensembl':
            transcript_name = feature[3].split(split_character)[0]
        else:
            transcript_name = 'NA'
            print('you are using an unknown reference platform. Please choose between refseq or ensembl')
        length = int(feature[2]) - int(feature[1])
        if not transcript_name in lengths:
            lengths[transcript_name] = 0
        lengths[transcript_name] += length
    return lengths


def iterate_over_folder(inputfolder, bedfile, platform, split_character):
    results = {}
    files = os.listdir(inputfolder)
    for current_file in files:
        if current_file.split('.')[-1] == 'bam':
            print(current_file)
            circle_coordinates = ['_'.join(current_file.split('_')[0:-3]), int(current_file.split('_')[-3]),
                                  int(current_file.split('_')[-2])]
            num_reads = int(current_file.split('_')[-1].split('.')[0].replace('reads', ''))
            mates, fragments = get_reads_from_bamfile('%s/%s' % (inputfolder, current_file), circle_coordinates)
            mates = classify_reads(mates)
            if not bedfile == 'none':
                length = annotate_circle(circle_coordinates, bedfile, platform, split_character)
            else:
                length = {}
            stats = get_statistics(mates)
            results[current_file.split('.')[0]] = stats
            if len(length) > 0:
                results[current_file.split('.')[0]]['min_length'] = min(length.items(), key=lambda x: x[1])[1]
                results[current_file.split('.')[0]]['max_length'] = max(length.items(), key=lambda x: x[1])[1]
                results[current_file.split('.')[0]]['transcript_ids'] = ','.join(length.keys())
            else:
                results[current_file.split('.')[0]]['min_length'] = circle_coordinates[2] - circle_coordinates[1]
                results[current_file.split('.')[0]]['max_length'] = circle_coordinates[2] - circle_coordinates[1]
                results[current_file.split('.')[0]]['transcript_ids'] = 'not_annotated'

            results[current_file.split('.')[0]]['circle_id'] = '%s_%s_%s' % (
                circle_coordinates[0], circle_coordinates[1], circle_coordinates[2])
            results[current_file.split('.')[0]]['num_reads'] = num_reads
    return results


def write_results(results, outfile):
    output_file = open(outfile, 'w')
    # eventually add gene name and length also with exons
    output_file.write('circle_id\ttranscript_ids\tnum_reads\tmin_length\tmax_length\tsingle\tdouble\tundefined\n')
    circles = sorted(results.keys())
    for circle in circles:
        output_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
            results[circle]['circle_id'], results[circle]['transcript_ids'], results[circle]['num_reads'],
            results[circle]['min_length'], results[circle]['max_length'], results[circle]['single'], results[circle]['double'],
            results[circle]['undefined']))
    output_file.close()
    return


# run script
if __name__ == '__main__':
    # required packages
    import pysam
    import os
    import argparse
    import pybedtools
    import tempfile

    parser = argparse.ArgumentParser(
        description='Extracts mate information and identify singe and double breakpoint fragments')

    # input
    parser.add_argument('bamfolder', metavar='PATH', help='path to folder containing circle bamfiles')
    # output
    parser.add_argument('outfile', metavar='outfile', help='path and filename to write the output to')
    # options
    parser.add_argument('-a', dest='bedfile', default='none',
                        help='if specified, the program will try to infer the circle length without internal introns')
    parser.add_argument('-p', dest='ref_platform', default='refseq',
                        help='specifies the annotation platform which was used (refseq or ensembl)')
    parser.add_argument('-s', dest='split_character', default='_',
                        help='specifies the separator within the name column in bedfile')
    parser.add_argument('--tmp', dest='tmp_folder', default='.',
                        help='tempfolder to store tempfiles generated by pybedtools.')

    args = parser.parse_args()

    # parse arguments
    bamfolder = args.bamfolder
    outfile = args.outfile
    bedfile = args.bedfile
    platform = args.ref_platform
    split_character = args.split_character
    tmp_folder = args.tmp_folder

    # set temp folder
    tempfile.tempdir = tmp_folder

    RESULTS = iterate_over_folder(bamfolder, bedfile, platform, split_character)
    write_results(RESULTS, outfile)
