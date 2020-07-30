# script to evaluate circle splicing variants

import pybedtools
import tempfile


class detect_splicing_variants(object):
    def __init__(self, split_character, platform, circles, bedfile, outfolder, sample, tmp_folder, cpus):

        # parse arguments
        self.circlefile = circles
        self.bedfile = bedfile
        self.outfile = outfolder + sample + ".alternative_splicing.txt"
        self.split_character = split_character
        self.platform = platform
        self.tmp_folder = tmp_folder
        self.cpus = cpus
        self.bed = pybedtools.example_bedtool(self.bedfile)

        # set temp folder
        tempfile.tempdir = tmp_folder
        pybedtools.set_tempdir(tmp_folder)

    def run_parallel(self, c):
        annotated_circles = {c: []}
        coordinates = pybedtools.BedTool('%s %s %s' % (c[0], c[1], c[2]), from_string=True)
        transcripts = self.bed.intersect(coordinates)
        for t in transcripts:
            if self.platform == 'refseq':
                tname = self.split_character.join(t[3].split(self.split_character)[0:2])
            else:
                tname = t[3].split(self.split_character)[0]
            annotated_circles[c] += [tname]
        annotated_circles[c] = set(annotated_circles[c])
        return annotated_circles

    def read_circle_file(self, infile):
        I = open(infile)
        circ_table = []
        for line in I:
            circle_id = line.split('\t')[0]
            chromosome = circle_id.split(':')[0]
            start = int(circle_id.split(':')[1].split('|')[0])
            end = int(circle_id.split(':')[1].split('|')[1])
            circ_table += [(chromosome, start, end)]
        I.close()
        return (circ_table)

    def annotate_circles(self, circles):

        from pathos.multiprocessing import ProcessingPool as Pool

        pool = Pool(self.cpus)
        tmp = pool.map(self.run_parallel, circles)

        new_dict = {}
        for item in tmp:
            for entry in item:
                new_dict[entry] = item[entry]

        return new_dict

    def accumulate_over_transcripts(self, circles):
        transcripts = {}
        for c in circles:
            for t in circles[c]:
                if not t in transcripts:
                    transcripts[t] = []
                transcripts[t] += [c]
        return (transcripts)

    def classify_multi_circle_transcripts(self, transcripts):
        classification = {}
        for lola in transcripts:
            types = {'same_start': {}, 'same_end': {}, 'within': {}, 'overlapping': {}, 'circles': []}
            circles = sorted(transcripts[lola])
            for i, circle1 in enumerate(circles):
                types['circles'] += ['%s:%s-%s' % (circle1[0], circle1[1], circle1[2])]
                for j, circle2 in enumerate(circles):
                    if i < j:
                        if circle1[1] == circle2[1]:
                            if not circle1[1] in types['same_start']:
                                types['same_start'][circle1[1]] = []
                            types['same_start'][circle1[1]] += ['%s:%s-%s' % (circle1[0], circle1[1], circle1[2]),
                                                                '%s:%s-%s' % (circle2[0], circle2[1], circle2[2])]
                        elif circle1[2] == circle2[2]:
                            if not circle1[2] in types['same_end']:
                                types['same_end'][circle1[2]] = []
                            types['same_end'][circle1[2]] += ['%s:%s-%s' % (circle1[0], circle1[1], circle1[2]),
                                                              '%s:%s-%s' % (circle2[0], circle2[1], circle2[2])]
                        elif (circle1[1] < circle2[1] and circle1[2] > circle2[2]) or (
                                        circle1[1] > circle2[1] and circle1[2] < circle2[2]):
                            if not i in types['within']:
                                types['within'][i] = []
                            types['within'][i] += ['%s:%s-%s' % (circle1[0], circle1[1], circle1[2]),
                                                   '%s:%s-%s' % (circle2[0], circle2[1], circle2[2])]
                        elif (circle1[1] < circle2[1] and circle1[2] < circle2[2] and circle1[2] > circle2[1]) or (
                                            circle1[1] > circle2[1] and circle1[2] > circle2[2] and circle1[1] <
                                    circle2[2]):
                            if not i in types['overlapping']:
                                types['overlapping'][i] = []
                            types['overlapping'][i] += ['%s:%s-%s' % (circle1[0], circle1[1], circle1[2]),
                                                        '%s:%s-%s' % (circle2[0], circle2[1], circle2[2])]
            classification[lola] = types
        return (classification)

    def write_genes(self, types, outfile):
        O = open(outfile, 'w')
        O.write('Transcript\tcircles\tsame_start\tsame_end\toverlapping\twithin\n')
        for lola in types:
            O.write('%s\t%s' % (lola, ','.join(types[lola]['circles'])))
            if len(types[lola]['same_start']) == 0:
                O.write('\t.')
            elif len(types[lola]['same_start']) > 0:
                O.write('\t')
                for circ in types[lola]['same_start']:
                    O.write('%s,' % ('|'.join(set(types[lola]['same_start'][circ]))))
            if len(types[lola]['same_end']) == 0:
                O.write('\t.')
            elif len(types[lola]['same_end']) > 0:
                O.write('\t')
                for circ in types[lola]['same_end']:
                    O.write('%s,' % ('|'.join(set(types[lola]['same_end'][circ]))))
            if len(types[lola]['overlapping']) == 0:
                O.write('\t.')
            elif len(types[lola]['overlapping']) > 0:
                O.write('\t')
                for circ in types[lola]['overlapping']:
                    O.write('%s,' % ('|'.join(set(types[lola]['overlapping'][circ]))))
            if len(types[lola]['within']) == 0:
                O.write('\t.')
            elif len(types[lola]['within']) > 0:
                O.write('\t')
                for circ in types[lola]['within']:
                    O.write('%s,' % ('|'.join(set(types[lola]['within'][circ]))))
            O.write('\n')
        O.close()
        return

    def run(self):

        # run
        circles = self.read_circle_file(self.circlefile)
        annotated_circles = self.annotate_circles(circles)
        transcripts = self.accumulate_over_transcripts(annotated_circles)
        circle_types = self.classify_multi_circle_transcripts(transcripts)
        self.write_genes(circle_types, self.outfile)
