# python script to get coverage profile for circle
# include some checks to make sure input was provided correctly

# required packages
import os
import argparse
import pybedtools
import tempfile


class get_coverage_profile(object):
    def __init__(self, exon_index, split_character, platform, bedfile, outfolder, sample, tmp_folder, cpus):

        # parse arguments
        self.bedfile = bedfile
        self.inputfolder = outfolder
        self.sample = sample
        self.exon_index = exon_index
        self.split_character = split_character
        self.platform = platform
        self.tmp_folder = tmp_folder
        self.exon_count_file = ""
        self.cpus = cpus

        # set temp folder
        tempfile.tempdir = tmp_folder
        pybedtools.set_tempdir(tmp_folder)

    # define functions
    def circle_exon_count(self, bamfile2, bedfile, exon_index, split_character, platform, coordinates):

        # does what I think it does, adjust to collapse different transcripts from the same gene,
        # choose transcript describing the circle best

        """
        """
        x = pybedtools.BedTool(bamfile2)
        b = pybedtools.BedTool(bedfile)

        b = b.filter(
            lambda b: b.chrom == coordinates[0] and b.start >= coordinates[1] - 1000 and b.end <= coordinates[2] + 1000)
        y = x.intersect(b, bed=True, wo=True, split=True)
        transcripts = {}
        found_features = []
        y = y.remove_invalid()

        for hit in y:

            if len(str(hit).split("\t")) < 19:
                print(("Malformed BED line: " + str(hit)))
                continue

            found_features += [hit[15]]
            transcript = hit[15]
            start = int(hit[13])
            end = int(hit[14])
            length = end - start
            strand_read = hit[5]
            strand_feature = hit[17]

            if platform == 'refseq':
                transcript_id = split_character.join(transcript.split(split_character)[0:2])
            elif platform == 'ensembl':
                transcript_id = transcript.split(split_character)[0]
            else:
                transcript_id = 'NA'
                print('you are using an unkown annotation platform, please use refseq or ensembl')

            if transcript.split(split_character)[exon_index].isdigit():
                exon = int(transcript.split(split_character)[exon_index])
            else:
                exon = 0

            read = hit[3]
            chromosome = hit[0]

            if not transcript_id in transcripts:
                transcripts[transcript_id] = {}
            if not exon in transcripts[transcript_id]:
                transcripts[transcript_id][exon] = {'length': length, 'start': start, 'end': end, 'strand_read': [],
                                                    'strand_feature': strand_feature, 'reads': [],
                                                    'chromosome': chromosome}

            transcripts[transcript_id][exon]['reads'] += [read]
            transcripts[transcript_id][exon]['strand_read'] += [strand_read]

        return transcripts, found_features

    def write_exon_count(self, outfile, exon_count, sample, circle_id,
                         transcript):  # append to existing exon_count file for the sample
        """
        """
        out = open(outfile, 'a')
        # sample\tcircle_id\ttranscript_id\texon_id\tchr\tstart\tend\tstrand\texon_length\tunique_reads\tfragments\tnumber+\tnumber-\n
        # sort exon ids per transcript..and then iterate from min to max, if one exon isn't in it, fill with 0's, identify potentially skipped exons
        for t in exon_count:
            if t == transcript:
                for exon in range(min(exon_count[transcript]), (max(exon_count[transcript]) + 1)):
                    if exon in exon_count[transcript]:
                        num_plus = exon_count[transcript][exon]['strand_read'].count('+')
                        num_minus = exon_count[transcript][exon]['strand_read'].count('-')
                        unique_reads = set([w.split('/')[0] for w in exon_count[transcript][exon]['reads']])
                        out.write('%s\t%s:%s-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                            sample, circle_id[0], circle_id[1], circle_id[2], transcript, ','.join(list(exon_count.keys())),
                            exon,
                            exon_count[transcript][exon]['chromosome'], exon_count[transcript][exon]['start'],
                            exon_count[transcript][exon]['end'], exon_count[transcript][exon]['strand_feature'],
                            exon_count[transcript][exon]['length'], len(unique_reads),
                            len(exon_count[transcript][exon]['reads']), num_plus, num_minus))
                    else:
                        out.write('%s\t%s:%s-%s\t%s\t%s\t%s\t0\t0\t0\t0\t0\t0\t0\t0\t0\n' % (
                            sample, circle_id[0], circle_id[1], circle_id[2], transcript, ','.join(list(exon_count.keys())),
                            exon))
        return

    def filter_features(self, bed_features, feature_names):
        """
        """
        intervals = ""
        for interval in bed_features:
            if interval[3] in feature_names:
                intervals += str(interval)
        return intervals

    def choose_transcript(self, exon_counts):
        """
        """
        if len(exon_counts) > 0:
            transcript = list(exon_counts.keys())[0]
            missing_exons_transcript = 100
            max_length_transcript = 0
            for t in exon_counts:
                missing_exons = 0
                max_length = 0
                for e in range(min(exon_counts[t]), max(exon_counts[t]) + 1):
                    if not e in exon_counts[t]:
                        missing_exons += 1
                    else:
                        max_length += len(exon_counts[t][e]['strand_read'])
                if len(exon_counts[t]) > len(exon_counts[transcript]):
                    transcript = t
                    missing_exons_transcript = missing_exons
                    max_length_transcript = max_length
                elif missing_exons < missing_exons_transcript:
                    transcript = t
                    missing_exons_transcript = missing_exons
                    max_length_transcript = max_length
                elif max_length_transcript < max_length:
                    transcript = t
                    max_length_transcript = max_length
                elif 'NR' in transcript and 'NM' in t:
                    transcript = t
                    missing_exons_transcript = missing_exons
                    max_length_transcript = max_length
        else:
            transcript = ''
        return transcript

    def circle_coverage_profile(self, bamfile, bedfile, exon_ind, split_character, platform):
        """
        """
        virtual_bed = pybedtools.BedTool(bedfile, from_string=True)
        bam = pybedtools.BedTool(bamfile)
        coverage = virtual_bed.coverage(bam, d=True, split=True)
        transcriptwise_coverage = {}
        for position in coverage:
            if platform == 'refseq':
                transcript = split_character.join(position[3].split(split_character)[0:2])
            elif platform == 'ensembl':
                transcript = position[3].split(split_character)[0]
            else:
                transcript = 'NA'
                print('you are using an unknown annotation platform, please use refseq or ensembl like formats')
            exon = int(position[3].split(split_character)[exon_ind])
            if not transcript in transcriptwise_coverage:
                transcriptwise_coverage[transcript] = {}
            if not exon in transcriptwise_coverage[transcript]:
                transcriptwise_coverage[transcript][exon] = {'relative_positions': [], 'position_coverage': [],
                                                             'chromosome': position[0], 'start': position[1],
                                                             'end': position[2]}
            transcriptwise_coverage[transcript][exon]['position_coverage'] += [position[7]]
            transcriptwise_coverage[transcript][exon]['relative_positions'] += [position[6]]
        return transcriptwise_coverage

    def write_coverage_profile(self, inputfolder, coverage_profile, sample, circle_id, transcript):
        """
        """
        for t in coverage_profile:
            if t == transcript:
                out = open('%s/%s.coverage_profiles/%s_%s_%s.%s.txt' % (
                    inputfolder, sample, circle_id[0], circle_id[1], circle_id[2], transcript), 'w')
                out.write('exon\trelative_pos_in_circle\trelative_pos_in_exon\tcoverage\n')
                pos_in_circle = 1
                for exon in range(min(coverage_profile[transcript]), (max(coverage_profile[transcript]) + 1)):
                    if exon in coverage_profile[transcript]:
                        for i, position in enumerate(coverage_profile[transcript][exon]['relative_positions']):
                            out.write('%s\t%s\t%s\t%s\n' % (
                                exon, pos_in_circle, position,
                                coverage_profile[transcript][exon]['position_coverage'][i]))
                            pos_in_circle += 1
                    else:
                        for i, position in enumerate(range(0, 50)):
                            out.write('%s\t%s\t%s\t0\n' % (exon, pos_in_circle, position))
                            pos_in_circle += 1
        out.close()
        return

    def remove_exons_outside_circle(self, exon_count, transcript, circle_id):
        exons = exon_count[transcript]
        keys_to_delete = []
        for exon in exons:
            start_exon = int(exon_count[transcript][exon]['start'])
            start_circle = circle_id[1]
            end_circle = circle_id[2]
            end_exon = int(exon_count[transcript][exon]['end'])
            if not (start_circle <= start_exon + 10 and end_exon - 10 <= end_circle):
                keys_to_delete += [exon]
        for exon in keys_to_delete:
            del exon_count[transcript][exon]
        return exon_count

    def format_to_bed12(self, exon_count, transcript, circle_id, number_of_reads,
                        outfile):  # correctly formatted now :)
        bed12 = {}
        for t in exon_count:
            if t == transcript and len(exon_count[t]) > 0:
                if circle_id[0].startswith('chr'):
                    bed12['01_chrom'] = circle_id[0]
                else:
                    bed12['01_chrom'] = 'chr%s' % (circle_id[0])
                bed12['02_start'] = '%s' % (circle_id[1])
                bed12['03_end'] = '%s' % (circle_id[2])
                bed12['04_name'] = t
                bed12['05_score'] = '%s' % (number_of_reads)
                bed12['06_strand'] = '.'
                bed12['07_thick_start'] = '%s' % (circle_id[1])
                bed12['08_thick_end'] = '%s' % (circle_id[2])
                bed12['09_itemRGB'] = '0,255,0'
                bed12['10_blockCount'] = '%s' % (len(exon_count[t]))
                bed12['11_block_sizes'] = []
                bed12['12_block_starts'] = []
                for e in sorted(exon_count[t]):
                    bed12['06_strand'] = exon_count[t][e]['strand_feature']
                    bed12['11_block_sizes'] += ['%s' % (exon_count[t][e]['length'] - 1)]
                    bed12['12_block_starts'] += ['%s' % (int(exon_count[t][e]['start'] + 1 - circle_id[1]))]
                # print(bed12)
                # if the exon doesn't start where the circle starts (circle starts in intron)
                if bed12['12_block_starts'][0] > '0':
                    bed12['12_block_starts'] = ['0'] + bed12['12_block_starts']
                    bed12['11_block_sizes'] = ['1'] + bed12['11_block_sizes']
                    bed12['10_blockCount'] = '%s' % (int(bed12['10_blockCount']) + 1)
                # if the exon starts before the circle starts (circle starts in exon)
                if bed12['12_block_starts'][0] < '0':
                    bed12['11_block_sizes'][0] = '%s' % (
                    int(bed12['11_block_sizes'][0]) + int(bed12['12_block_starts'][0]))
                    bed12['12_block_starts'][0] = '0'
                # if the last exon extends over the circle boundaries:
                if int(bed12['12_block_starts'][-1]) + int(bed12['02_start']) + int(bed12['11_block_sizes'][-1]) > int(
                        bed12['03_end']):
                    actual_end = int(bed12['12_block_starts'][-1]) + int(bed12['02_start']) + int(
                        bed12['11_block_sizes'][-1])
                    difference = actual_end - int(bed12['03_end'])
                    new_blocksize = int(bed12['11_block_sizes'][-1]) - difference
                    bed12['11_block_sizes'][-1] = '%s' % (new_blocksize)
                if int(bed12['12_block_starts'][-1]) + int(bed12['02_start']) + int(bed12['11_block_sizes'][-1]) < int(
                        bed12['03_end']):
                    bed12['11_block_sizes'] += '1'
                    bed12['12_block_starts'] += ['%s' % (int(bed12['03_end']) - int(bed12['02_start']) - 1)]
                    bed12['10_blockCount'] = '%s' % (int(bed12['10_blockCount']) + 1)
                bed12['12_block_starts'] = ','.join(bed12['12_block_starts'])
                bed12['11_block_sizes'] = ','.join(bed12['11_block_sizes'])
        Bed12 = []
        for i in sorted(bed12):
            Bed12 += [bed12[i]]
        o = open(outfile, 'a')
        o.write('%s\n' % ('\t'.join(Bed12)))
        o.close()
        return bed12

    def run_parallel(self, f):
        if f.split('.')[-2] == 'sorted':
            # extract circle id from filename, works for files
            # generated by extract_reads.py, consider making this more flexible
            circle_id = ('_'.join(f.split('_')[0:-3]), int(f.split('_')[-3]), int(f.split('_')[-2]))
            number_of_reads = int(f.split('_')[-1].split('.')[0].replace('reads', ''))
            bamfile2 = '%s/%s/%s' % (self.inputfolder, self.sample, f)
            # open bed feature file
            b = pybedtools.BedTool(self.bedfile)
            # get read counts for each exon in circle

            exon_counts, found_features = self.circle_exon_count(bamfile2, self.bedfile, self.exon_index,
                                                                 self.split_character, self.platform, circle_id)
            if len(exon_counts) > 0:
                # choose best fitting transcript
                transcript_id = self.choose_transcript(exon_counts)
                # add circle to result table
                self.write_exon_count(self.exon_count_file, exon_counts, self.sample, circle_id, transcript_id)
                exon_counts = self.remove_exons_outside_circle(exon_counts, transcript_id, circle_id)
                self.format_to_bed12(exon_counts, transcript_id, circle_id, number_of_reads,
                                     '%s/%s.exon_counts.bed' % (self.inputfolder, self.sample))
                filtered_features = self.filter_features(b, found_features)
                if len(filtered_features) > 0:
                    coverage_track = self.circle_coverage_profile(bamfile2, filtered_features, self.exon_index,
                                                                  self.split_character,
                                                                  self.platform)
                    self.write_coverage_profile(self.inputfolder, coverage_track, self.sample, circle_id, transcript_id)

    # circle exon count over all bam files in sample folder, this could easily be parralellised
    # also think about using our denovo recontruction as bases for coverage profiles

    def run(self):

        # initializing the result table file
        self.exon_count_file = '%s/%s.exon_counts.txt' % (self.inputfolder, self.sample)
        exon_counts_out = open(self.exon_count_file, 'w')
        exon_counts_out.write('sample\tcircle_id\ttranscript_id\tother_ids\texon_id\tchr\tstart'
                              '\tend\tstrand\texon_length\tunique_reads\tfragments\tnumber+\tnumber-\n')
        exon_counts_out.close()

        output_file = open('%s/%s.exon_counts.bed' % (self.inputfolder, self.sample), 'w')
        output_file.write('# BED12\n')
        output_file.close()

        # all circle files in a given folder
        files = os.listdir('%s/%s' % (self.inputfolder, self.sample))

        # create folder for coverage profiles
        folders = os.listdir(self.inputfolder)
        if not '%s.coverage_profiles' % (self.sample) in folders:
            os.mkdir('%s/%s.coverage_profiles' % (self.inputfolder, self.sample))

        from pathos.multiprocessing import ProcessingPool as Pool

        pool = Pool(self.cpus)
        pool.map(self.run_parallel, files)
