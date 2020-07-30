# required packages
import os
import tempfile

import pybedtools
import pysam


# script to identify circles where both mates map over the junction
# length estimate is not correct for circRNAs exceeding the host-gene coordinates. need to fix this!


class mate_information(object):

    def run_parallel(self, f):

        internal_dict = {}

        if f.split('.')[-1] == 'bam':
            circle_coordinates = ['_'.join(f.split('_')[0:-3]), int(f.split('_')[-3]),
                                  int(f.split('_')[-2])]
            num_reads = int(f.split('_')[-1].split('.')[0].replace('reads', ''))
            mates, fragments = self.get_reads_from_bamfile('%s/%s' % (self.bamfolder, f), circle_coordinates)
            mates = self.classify_reads(mates)
            if not self.bedfile == 'none':
                length = self.annotate_circle(circle_coordinates, self.bedfile, self.platform, self.split_character)
            else:
                length = {}
            stats = self.get_statistics(mates)
            internal_dict[f.split('.')[0]] = stats
            if len(length) > 0:
                internal_dict[f.split('.')[0]]['min_length'] = min(list(length.items()), key=lambda x: x[1])[1]
                internal_dict[f.split('.')[0]]['max_length'] = max(list(length.items()), key=lambda x: x[1])[1]
                internal_dict[f.split('.')[0]]['transcript_ids'] = ','.join(list(length.keys()))
            else:
                internal_dict[f.split('.')[0]]['min_length'] = circle_coordinates[2] - circle_coordinates[1]
                internal_dict[f.split('.')[0]]['max_length'] = circle_coordinates[2] - circle_coordinates[1]
                internal_dict[f.split('.')[0]]['transcript_ids'] = 'not_annotated'

            internal_dict[f.split('.')[0]]['circle_id'] = '%s_%s_%s' % (
                circle_coordinates[0], circle_coordinates[1], circle_coordinates[2])
            internal_dict[f.split('.')[0]]['num_reads'] = num_reads

        return internal_dict


    def __init__(self, platform, split_character, bedfile, outfolder, sample, tmp_folder, cpus):

        # parse arguments
        self.bamfolder = outfolder + sample
        self.outfile = outfolder + sample + ".mate_status.txt"
        self.bedfile = bedfile
        self.platform = platform
        self.split_character = split_character
        self.tmp_folder = tmp_folder
        self.internal_dict = {}
        self.cpus = cpus

        # set temp folder
        tempfile.tempdir = tmp_folder
        pybedtools.set_tempdir(tmp_folder)

    # define functions
    def get_reads_from_bamfile(self, bamfile, circle_coordinates):
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

    def classify_reads(self, mates):
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

    def get_statistics(self, mates):
        stats = {'single': 0, 'double': 0, 'undefined': 0}
        for mate in mates:
            stats[mates[mate]['status']] += 1
        return stats

    def annotate_circle(self, circle_coordinates, bedfile, platform, split_character):
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

    def iterate_over_folder(self, inputfolder, bedfile, platform, split_character):
        files = os.listdir(inputfolder)

        from pathos.multiprocessing import ProcessingPool as Pool

        p = Pool(self.cpus)
        tmp = p.map(self.run_parallel, files)

        new_dict = {}
        for item in tmp:
            for entry in item:
                new_dict[entry] = item[entry]

        return new_dict

    def write_results(self, results, outfile):
        output_file = open(outfile, 'w')
        # eventually add gene name and length also with exons
        output_file.write('circle_id\ttranscript_ids\tnum_reads\tmin_length\tmax_length\tsingle\tdouble\tundefined\n')
        circles = sorted(results.keys())
        for circle in circles:
            output_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                results[circle]['circle_id'], results[circle]['transcript_ids'], results[circle]['num_reads'],
                results[circle]['min_length'], results[circle]['max_length'], results[circle]['single'],
                results[circle]['double'],
                results[circle]['undefined']))
        output_file.close()
        return

    # run script
    def run(self):

        RESULTS = self.iterate_over_folder(self.bamfolder, self.bedfile, self.platform, self.split_character)
        self.write_results(RESULTS, self.outfile)
