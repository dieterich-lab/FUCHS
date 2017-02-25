#! /usr/bin/env python

# script to identify skipped exons in circRNA bam files


import pysam
import pybedtools
import tempfile
import os

class detect_skipped_exons(object):
    def __init__(self, outfolder, sample, bedfile, tmp_folder, platform):

        self.folder = outfolder + sample
        self.sample = sample
        self.outfile = outfolder + sample + ".skipped_exons.txt"
        self.bedfile = bedfile
        self.tmp_folder = tmp_folder
        self.platform = platform

        tempfile.tempdir = tmp_folder
        pybedtools.set_tempdir(tmp_folder)

        # define functions
    def load_bamfile(self, bamfile):
        reads = {}
        B = pysam.AlignmentFile(bamfile, 'rb')
        for i, lola in enumerate(B):
            if not -1 in lola.get_tag('jI'):
                if not lola.query_name in reads:
                    reads[lola.query_name] = {}
                reads[lola.query_name][i] = {'reference': B.getrname(lola.reference_id), 'breakpoint': lola.get_tag('jI'),
                                             'mapq': lola.mapping_quality}
        B.close()
        return (reads)


    def filter_reads(self, reads):
        for lola in reads:
            if len(reads[lola]) > 1:
                occurences = reads[lola]
                KEYS = occurences.keys()
                mapq = occurences[KEYS[0]]['mapq']
                for forrest in KEYS:
                    if occurences[forrest]['mapq'] > mapq:
                        mapq = occurences[forrest]['mapq']
                for forrest in KEYS:
                    if occurences[forrest]['mapq'] < mapq:
                        del occurences[forrest]
                reads[lola] = occurences
        return (reads)


    def intersect_introns_with_bedfile(self, bedfile, reads, coordinates):
        exons = pybedtools.example_bedtool(bedfile)
        exons = exons.filter(lambda b: b.chrom == coordinates[0] and b.start >= coordinates[1] and b.end <= coordinates[2])
        skipped_exons = {}
        for lola in reads:
            for forrest in reads[lola]:
                reads[lola][forrest]['intron'] = {}
                breakpoints = reads[lola][forrest]['breakpoint']
                starts = breakpoints[::2]
                ends = breakpoints[1::2]
                for i, start in enumerate(starts):
                    intron = pybedtools.BedTool('%s %s %s' % (reads[lola][forrest]['reference'], start, ends[i]),
                                                from_string=True)
                    exons = pybedtools.example_bedtool(bedfile)
                    features = exons.intersect(intron, stream = True)
                    for skipped in features:
                        if not (skipped[0], int(skipped[1]), int(skipped[2])) in skipped_exons:
                            skipped_exons[(skipped[0], int(skipped[1]), int(skipped[2]))] = {'reads': [], 'intron': [],
                                                                                             'name': skipped[3]}
                        skipped_exons[(skipped[0], int(skipped[1]), int(skipped[2]))]['reads'] += [lola]
                        skipped_exons[(skipped[0], int(skipped[1]), int(skipped[2]))]['intron'] += [
                            (reads[lola][forrest]['reference'], start, ends[i])]
        return (skipped_exons)


    def identify_skipped_exons(self, bamfile, skipped_exons):
        bam = pybedtools.example_bedtool(bamfile)
        for lola in skipped_exons:
            exon_readcount = []
            coordinates = pybedtools.BedTool('%s %s %s' % (lola[0], lola[1], lola[2]), from_string=True)
            reads = bam.intersect(coordinates, stream = True)
            for r in reads:
                exon_readcount += [r[0]]
            exon_readcount = set(exon_readcount)
            skipped_exons[lola]['exon_readcount'] = len(exon_readcount)
        return (skipped_exons)


    def write_exon_skipping(self,skipped, outfile, circle_id, platform):
        O = open(outfile, 'a')
        for exon in skipped:
            if not len(set(skipped[exon]['reads'])) == skipped[exon]['exon_readcount']:
                if platform == 'refseq':
                    name = '_'.join(skipped[exon]['name'].split('_')[0:2])
                else:
                    name = skipped[exon]['name'].split('_')[0]
                O.write('%s_%s_%s\t%s\t%s:%s-%s\t%s\t%s\t%s\t%s\n' % (
                circle_id[0], circle_id[1], circle_id[2], name, exon[0], exon[1], exon[2], set(skipped[exon]['intron']),
                ','.join(set(skipped[exon]['reads'])), len(set(skipped[exon]['reads'])), skipped[exon]['exon_readcount']))
        O.close()
        return


    def write_bed12(self, skipped, outfile, circle_id, platform):
        O = open(outfile, 'a')
        for exon in skipped:
            if not len(set(skipped[exon]['reads'])) == skipped[exon]['exon_readcount'] and (
                    circle_id[2] - circle_id[1] > 100) and (circle_id[2] >= skipped[exon]['intron'][0][2]) and (
                circle_id[1] <= skipped[exon]['intron'][0][1]):
                if platform == 'refseq':
                    name = '_'.join(skipped[exon]['name'].split('_')[0:2])
                else:
                    name = skipped[exon]['name'].split('_')[0]
                if circle_id[0].startswith('chr'):
                    O.write('%s\t%s\t%s\t%s\t%s\t.\t%s\t%s\t255,0,0\t3\t1,%s,1\t0,%s,%s\n' % (
                    circle_id[0], circle_id[1], circle_id[2], name,
                    (float(len(set(skipped[exon]['reads']))) / skipped[exon]['exon_readcount']) * 100,
                    skipped[exon]['intron'][0][1], skipped[exon]['intron'][0][2], (exon[2] - exon[1]),
                    exon[1] - circle_id[1], circle_id[2] - circle_id[1] - 1))
                else:
                    O.write('chr%s\t%s\t%s\t%s\t%s\t.\t%s\t%s\t255,0,0\t3\t1,%s,1\t0,%s,%s\n' % (
                    circle_id[0], circle_id[1], circle_id[2], name,
                    (float(len(set(skipped[exon]['reads']))) / skipped[exon]['exon_readcount']) * 100,
                    skipped[exon]['intron'][0][1], skipped[exon]['intron'][0][2], (exon[2] - exon[1]),
                    exon[1] - circle_id[1], circle_id[2] - circle_id[1] - 1))
        O.close()
        return


    def run(self):

        files = os.listdir(self.folder)
        outfile_bed = self.outfile.replace('.txt', '.bed')

        O = open(self.outfile, 'w')
        O.write('circle_id\ttranscript_id\tskipped_exon\tintron\tread_names\tsplice_reads\texon_reads\n')
        O.close()

        O = open(outfile_bed, 'w')
        O.write('# bed12 format\n')
        O.close()

        for f in files:
            if f.split('.')[-1] == 'bam':
                print(f)
                circle_id = ('_'.join(f.split('_')[0:-3]), int(f.split('_')[-3]), int(f.split('_')[-2]))
                bamfile = '%s/%s' % (self.folder, f)
                READS = self.load_bamfile(bamfile)
                READS = self.filter_reads(READS)
                SKIPPED = self.intersect_introns_with_bedfile(self.bedfile, READS, circle_id)
                if len(SKIPPED) > 0:
                    SKIPPED = self.identify_skipped_exons(bamfile, SKIPPED)
                    self.write_bed12(SKIPPED, outfile_bed, circle_id, self.platform)
                    self.write_exon_skipping(SKIPPED, self.outfile, circle_id, self.platform)
