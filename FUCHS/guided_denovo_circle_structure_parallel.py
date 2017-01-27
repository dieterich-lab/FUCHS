#! /usr/bin/env python


# define functions

def load_bamfile(bamfile, coordinates):
    reads = {}
    B = pysam.AlignmentFile(bamfile, 'rb')
    for i, lola in enumerate(B):
        if not lola.get_tag('jI')[0] == -1 and B.getrname(lola.reference_id) == coordinates[0]:
            if not lola.query_name in reads:
                reads[lola.query_name] = {}
            breakpoints = []
            for b in lola.get_tag('jI'):
                breakpoints += [b]
            reads[lola.query_name][i] = {'reference': B.getrname(lola.reference_id), 'breakpoint': breakpoints,
                                         'mapq': lola.mapping_quality}
    B.close()
    return (reads)


def filter_reads(reads, coordinates):
    reads_to_delte = []
    for lola in reads:
        if len(reads[lola]) > 1:
            occurences = reads[lola]
            KEYS = occurences.keys()
            mapq = occurences[KEYS[0]]['mapq']
            for forrest in KEYS:
                if not occurences[forrest]['reference'] == coordinates[0]:
                    del occurences[forrest]
            KEYS = occurences.keys()
            for forrest in KEYS:
                if occurences[forrest]['mapq'] > mapq:
                    mapq = occurences[forrest]['mapq']
            for forrest in KEYS:
                if occurences[forrest]['mapq'] < mapq:
                    del occurences[forrest]
            reads[lola] = occurences
    return (reads)


def get_introns(reads):
    introns = {}
    for lola in reads:
        for forrest in reads[lola]:
            breakpoints = reads[lola][forrest]['breakpoint']
            starts = breakpoints[::2]
            ends = breakpoints[1::2]
            for i, start in enumerate(starts):
                if not (reads[lola][forrest]['reference'], start, ends[i]) in introns:
                    introns[(reads[lola][forrest]['reference'], start, ends[i])] = {'spanning_reads': [],
                                                                                    'skipped_exons': {}}
                introns[(reads[lola][forrest]['reference'], start, ends[i])]['spanning_reads'] += [lola]
    return (introns)


def connect_introns(introns, circ_coordinates):
    transcripts = {0: [(circ_coordinates[0], circ_coordinates[1] - 1, circ_coordinates[1])]}
    sorted_introns = sorted(introns.keys())
    for i in sorted_introns:
        # print(i)
        connected = False
        for t in transcripts:
            if i[1] > transcripts[t][-1][2]:
                transcripts[t] += [i]
                connected = True
        if not connected:
            st = sorted(transcripts)
            for t in st:
                ti = max(transcripts) + 1
                transcripts[ti] = transcripts[t][:-1] + [i]
    for t in transcripts:
        transcripts[t] += [(circ_coordinates[0], circ_coordinates[2] - 1, circ_coordinates[2])]
    tmp = {}
    for key, value in transcripts.items():
        if value not in tmp.values():
            tmp[key] = value
    transcripts = tmp
    return (transcripts)


def get_coverage_profile(bamfile, circ_coordinates, transcripts):
    bam = pybedtools.example_bedtool(bamfile)
    coordinates = pybedtools.BedTool('%s %s %s' % (circ_coordinates[0], circ_coordinates[1], circ_coordinates[2]),
                                     from_string=True)
    circle_coverage = coordinates.coverage(bam, d=True)
    split_circle_coverage = coordinates.coverage(bam, d=True, split=True)
    coverage = []
    split_coverage = []
    for lola in circle_coverage:
        coverage += [int(lola[4])]
    for lola in split_circle_coverage:
        split_coverage += [int(lola[4])]
    transcript_coverage = {}
    for t in transcripts:
        transcript_coverage[t] = {'exons': {}, 'introns': {}}
        if 0 in coverage:
            transcript_coverage[t]['coverage_breaks'] = [(circ_coordinates[0], circ_coordinates[1] + coverage.index(0),
                                                          circ_coordinates[1] + len(coverage) - (
                                                          list(reversed(coverage)).index(0)))]
        else:
            transcript_coverage[t]['coverage_breaks'] = []
        for i, intron in enumerate(transcripts[t]):
            if i > 0 and (i + 1) < len(transcripts[t]):
                transcript_coverage[t]['introns'][intron] = sum(
                    coverage[intron[1] - circ_coordinates[1]:intron[2] - circ_coordinates[1]]) / float(
                    intron[2] - intron[1])
            if not (i + 1) == len(transcripts[t]):
                transcript_coverage[t]['exons'][(intron[0], intron[2], transcripts[t][i + 1][1])] = sum(
                    coverage[intron[2] - circ_coordinates[1]:transcripts[t][i + 1][1] - circ_coordinates[1]]) / float(
                    transcripts[t][i + 1][1] - intron[2])
    return (transcript_coverage, coverage, split_coverage)


def filter_out_exons(transcript_coverage, split_coverage, circ_coordinates):  # also keep right side of exon
    for t in transcript_coverage:
        exons_to_remove = []
        tmp_exons = {}
        for exon in transcript_coverage[t]['exons']:
            relative_start = exon[1] - circ_coordinates[1]
            relative_end = exon[2] - circ_coordinates[1] - 1
            coverage = split_coverage[relative_start:relative_end]
            if sum(coverage) > 0:
                if 0 in coverage and not (coverage[0] == 0 or coverage[-1] == 0):
                    breakpoint = coverage.index(0)
                    new_end = exon[1] + coverage.index(0) - 1
                    new_start = exon[2] - (list(reversed(coverage)).index(0))
                    avg_coverage = sum(coverage[:coverage.index(0)]) / float(coverage.index(0))
                    avg_coverage_right = sum(list(reversed(coverage))[:list(reversed(coverage)).index(0)]) / float(
                        list(reversed(coverage)).index(0))
                    exons_to_remove += [exon]
                    tmp_exons[(exon[0], exon[1], new_end)] = avg_coverage
                    tmp_exons[(exon[0], new_start, exon[2])] = avg_coverage_right
                elif 0 in coverage and coverage[0] == 0:
                    breakpoint = coverage.index(0)
                    new_end = exon[1] + coverage.index(0) - 1
                    new_start = exon[2] - (list(reversed(coverage)).index(0))
                    avg_coverage_right = sum(list(reversed(coverage))[:list(reversed(coverage)).index(0)]) / float(
                        list(reversed(coverage)).index(0))
                    exons_to_remove += [exon]
                    tmp_exons[(exon[0], new_start, exon[2])] = avg_coverage_right
                elif 0 in coverage and coverage[-1] == 0:
                    breakpoint = coverage.index(0)
                    new_end = exon[1] + coverage.index(0) - 1
                    new_start = exon[2] - (list(reversed(coverage)).index(0))
                    avg_coverage = sum(coverage[:coverage.index(0)]) / float(coverage.index(0))
                    exons_to_remove += [exon]
                    tmp_exons[(exon[0], exon[1], new_end)] = avg_coverage
                else:
                    avg_coverage = sum(coverage) / float(len(coverage))
                    transcript_coverage[t]['exons'][exon] = avg_coverage
            else:
                exons_to_remove += [exon]
        for exon in exons_to_remove:
            del transcript_coverage[t]['exons'][exon]
        for exon in tmp_exons:
            transcript_coverage[t]['exons'][exon] = tmp_exons[exon]
    return (transcript_coverage)


def collapse_transcripts(transcript_coverage):
    duplicated_transcripts = []
    for t, transcript in enumerate(transcript_coverage):
        for i, transcript2 in enumerate(transcript_coverage):
            if i > t:
                # print(sorted(transcript_coverage[transcript]['exons'].keys()))
                # print(sorted(transcript_coverage[transcript2]['exons'].keys()))
                # print('\n\n')
                if sorted(transcript_coverage[transcript]['exons'].keys()) == sorted(
                        transcript_coverage[transcript2]['exons'].keys()):
                    duplicated_transcripts += [transcript]
    for t in duplicated_transcripts:
        del transcript_coverage[t]
    return (transcript_coverage)


def infer_missing_structure(transcripts, coordinates, bedfile):
    B = pybedtools.example_bedtool(bedfile)
    # B = B.filter(lambda b: b.chrom == coordinates[0] and b.start >= coordinates[1] and b.end <= coordinates[2])
    for t in transcripts:
        if len(transcripts[t]['coverage_breaks']) > 0:
            for unsupported in transcripts[t]['coverage_breaks']:
                missing_region = pybedtools.BedTool('%s %s %s' % (unsupported[0], unsupported[1], unsupported[2]),
                                                    from_string=True)
                features = B.intersect(missing_region)
                unsuported_exons = {}
                for f in features:
                    transcripts[t]['exons'][(f[0], int(f[1]), int(f[2]))] = 0
        merged = [1]
        counter = 0
        while len(merged) > 0 or counter < 100:
            transcripts[t]['exons'], merged = merge_exons(transcripts[t]['exons'])
            counter += 1
    return (transcripts)


def merge_exons(exons):
    sorted_exons = sorted(exons)
    new_exons = {}
    merged = []
    if len(sorted_exons) > 1:
        for i, e in enumerate(sorted_exons):
            if not e == sorted_exons[-1] and not e in merged:
                next_exon = sorted_exons[i + 1]
                if e[2] + 1 >= next_exon[1]:
                    new_exons[(e[0], e[1], next_exon[2])] = ((exons[e] * (e[2] - e[1])) + (
                    exons[next_exon] * (next_exon[2] - next_exon[1]))) / float(next_exon[2] - e[1])
                    merged = [e, next_exon]
                else:
                    new_exons[e] = exons[e]
            elif e == sorted_exons[-1] and e not in merged:
                prev_exon = sorted_exons[i - 1]
                if e[1] - 1 <= prev_exon[2]:
                    new_exons[(e[0], prev_exon[1], e[2])] = ((exons[e] * (e[2] - e[1])) + (
                    exons[prev_exon] * (prev_exon[2] - prev_exon[1]))) / float(e[2] - prev_exon[1])
                    merged = [e, prev_exon]
                else:
                    new_exons[e] = exons[e]
    else:
        new_exons = exons
    return (new_exons, merged)


def write_bed12(outfile, transcript_coverage, circ_coordinates, coverage, introns):
    O = open(outfile, 'a')
    for t in transcript_coverage:
        O.write('%s\t%s\t%s\t%s:%s-%s|%s|%s\t' % (
        circ_coordinates[0], circ_coordinates[1], circ_coordinates[2], circ_coordinates[0], circ_coordinates[1],
        circ_coordinates[2], t, 1 - (coverage.count(0) / float(len(coverage)))))
        transcript_confidence = []
        if len(transcript_coverage[t]['introns']) > 0:
            for i in transcript_coverage[t]['introns']:
                transcript_confidence += [len(introns[i]['spanning_reads'])]
            O.write('%s\t.\t%s\t%s\t255,0,0\t%s\t' % (
            int(sum(transcript_confidence) / float(len(transcript_confidence))), circ_coordinates[1],
            circ_coordinates[2], len(transcript_coverage[t]['exons'])))
        else:
            if len(transcript_coverage[t]['exons']) > 0:
                O.write('%s\t.\t%s\t%s\t255,0,0\t%s\t' % (
                int(sum(transcript_coverage[t]['exons'].values()) / float(len(transcript_coverage[t]['exons']))),
                circ_coordinates[1], circ_coordinates[2], len(transcript_coverage[t]['exons'])))
            else:
                O.write('0\t.\t%s\t%s\t255,0,0\t%s\t' % (
                circ_coordinates[1], circ_coordinates[2], len(transcript_coverage[t]['exons'])))
        exon_length = []
        exon_location = []
        sorted_exons = sorted(transcript_coverage[t]['exons'])
        for e in sorted_exons:
            exon_length += ['%s' % (e[2] - e[1] + 1)]
            exon_location += ['%s' % (e[1] - circ_coordinates[1])]
        O.write('%s\t%s\n' % (','.join(exon_length), ','.join(exon_location)))
    O.close()
    return


def write_bed6(transcript_coverage, outfile, circ_coordinates, coverage):
    exons = {}
    O = open(outfile, 'a')
    for t in transcript_coverage:
        for e in transcript_coverage[t]['exons']:
            if not e in exons:
                exons[e] = ['%s' % (t)]
            else:
                exons[e] += ['%s' % (t)]
    sorted_exons = sorted(exons.keys())
    for e in sorted_exons:
        O.write('%s\t%s\t%s\t%s:%s-%s|%s\t%s\t.\n' % (
        e[0], e[1], e[2], circ_coordinates[0], circ_coordinates[1], circ_coordinates[2], ','.join(exons[e]),
        int(sum(coverage[e[1] - circ_coordinates[1]:e[2] - circ_coordinates[1]]) / float(e[2] - e[1]))))
    O.close()
    return


def get_coverage(circ_coordinates, bamfile):
    bam = pybedtools.example_bedtool(bamfile)
    coordinates = pybedtools.BedTool('%s %s %s' % (circ_coordinates[0], circ_coordinates[1], circ_coordinates[2]),
                                     from_string=True)

    circle_coverage = coordinates.coverage(bam, d=True)
    coverage = []
    for lola in circle_coverage:
        coverage += [int(lola[4])]
    return (coverage)


def write_single_exon(outfile, coverage, circ_coordinates):
    O12 = open('%sinferred_12.bed' % (outfile), 'a')
    O6 = open('%sinferred_6.bed' % (outfile), 'a')
    if 0 in coverage and sum(coverage) > 0:
        breakpoints = (
        circ_coordinates[1] + coverage.index(0), circ_coordinates[2] - (list(reversed(coverage)).index(0)))
        exon1 = (circ_coordinates[0], circ_coordinates[1], breakpoints[0])
        exon2 = (circ_coordinates[0], breakpoints[1], circ_coordinates[2])
        O12.write('%s\t%s\t%s\t%s:%s-%s|0|%s\t%s\t.\t%s\t%s\t255,0,0\t2\t%s,%s\t0,%s\n' % (
        circ_coordinates[0], circ_coordinates[1], circ_coordinates[2], circ_coordinates[0], circ_coordinates[1],
        circ_coordinates[2], 1 - (coverage.count(0) / float(len(coverage))), int((sum(
            coverage[:coverage.index(0)]) / float(coverage.index(0)) + sum(
            coverage[breakpoints[1] - circ_coordinates[1]:]) / float(circ_coordinates[2] - breakpoints[1])) / 2),
        circ_coordinates[1], circ_coordinates[2], coverage.index(0), circ_coordinates[2] - breakpoints[1],
        breakpoints[1] - circ_coordinates[1]))
        O6.write('%s\t%s\t%s\t%s:%s-%s|0\t%s\t.\n' % (
        exon1[0], exon1[1], exon1[2], circ_coordinates[0], circ_coordinates[1], circ_coordinates[2],
        int(sum(coverage[:coverage.index(0)]) / float(coverage.index(0)))))
        O6.write('%s\t%s\t%s\t%s:%s-%s|0\t%s\t.\n' % (
        exon2[0], exon2[1], exon2[2], circ_coordinates[0], circ_coordinates[1], circ_coordinates[2],
        int(sum(coverage[breakpoints[1] - circ_coordinates[1]:]) / float(circ_coordinates[2] - breakpoints[1]))))
    elif not 0 in coverage:
        O12.write('%s\t%s\t%s\t%s:%s-%s|0|%s\t%s\t.\t%s\t%s\t255,0,0\t1\t%s\t0\n' % (
        circ_coordinates[0], circ_coordinates[1], circ_coordinates[2], circ_coordinates[0], circ_coordinates[1],
        circ_coordinates[2], 1 - (coverage.count(0) / float(len(coverage))), int(sum(coverage) / float(len(coverage))),
        circ_coordinates[1], circ_coordinates[2], circ_coordinates[2] - circ_coordinates[1]))
        O6.write('%s\t%s\t%s\t%s:%s-%s|0\t%s\t.\n' % (
        circ_coordinates[0], circ_coordinates[1], circ_coordinates[2], circ_coordinates[0], circ_coordinates[1],
        circ_coordinates[2], int(sum(coverage) / float(len(coverage)))))
    else:
        print('no reads for circle %s:%s-%s' % (circ_coordinates[0], circ_coordinates[1], circ_coordinates[2]))
    O6.close()
    O12.close()
    return


def run_denovo_exon_chain_reconstruction(f, folder, annotation, outfile):
    if f.split('.')[-1] == 'bam':
        print(f)
        # extract generic information from filename
        circ_coordinates = ('_'.join(f.split('_')[0:-3]), int(f.split('_')[-3]), int(f.split('_')[-2]))
        num_reads = int(f.split('_')[-1].split('.')[0].replace('reads', ''))
        bamfile = '%s/%s' % (folder, f)
        # load bamfile
        READS = load_bamfile(bamfile, circ_coordinates)
        READS = filter_reads(READS, circ_coordinates)
        Introns = get_introns(READS)
        if len(Introns) > 0:
            # reconstruct transcripts
            T = connect_introns(Introns, circ_coordinates)
            # get coverage for Transcripts
            TC, Cov, splitCov = get_coverage_profile(bamfile, circ_coordinates, T)
            TC = filter_out_exons(TC, splitCov, circ_coordinates)
            if not annotation == '.':
                TC = infer_missing_structure(TC, circ_coordinates, annotation)
            # write out results to 3 different files, this will probably have to be adjuste
            write_bed12('%sinferred_12.bed' % (outfile), TC, circ_coordinates, Cov, Introns)
            write_bed6(TC, '%sinferred_6.bed' % (outfile), circ_coordinates, splitCov)
        else:
            Cov = get_coverage(circ_coordinates, bamfile)
            if not 0 in Cov or annotation == '.':
                write_single_exon(outfile, Cov, circ_coordinates)
            elif sum(Cov) > 0:
                breakpoints = (circ_coordinates[1] + Cov.index(0), circ_coordinates[2] - (list(reversed(Cov)).index(0)))
                TC = {
                    0: {'introns': Introns, 'coverage_breaks': [(circ_coordinates[0], breakpoints[0], breakpoints[1])],
                        'exons': {
                            (circ_coordinates[0], circ_coordinates[1], breakpoints[0]): sum(Cov[:Cov.index(0)]) / float(
                                Cov.index(0)), (circ_coordinates[0], breakpoints[1], circ_coordinates[2]): sum(
                                Cov[breakpoints[1] - circ_coordinates[1]:]) / float(
                                circ_coordinates[2] - breakpoints[1])}}}
                TC = infer_missing_structure(TC, circ_coordinates, annotation)
                write_bed12('%sinferred_12.bed' % (outfile), TC, circ_coordinates, Cov, Introns)
                write_bed6(TC, '%sinferred_6.bed' % (outfile), circ_coordinates, Cov)
            else:
                print('no reads mapped for %s' % (f))
    return (f, len(Introns))


# Run script

if __name__ == '__main__':
    # Handle command line options

    # required packages
    import pysam
    import pybedtools
    import argparse
    import tempfile

    import os
    import multiprocessing

    parser = argparse.ArgumentParser(
        description='file to run a denonvo exon chain reconstruction. circRNAs are independent from each other and script can be run on multiple cores')
    # input
    parser.add_argument('inputfolder', metavar='folder',
                        help='folder containing all circle bam files. (full path, but without sample name)')
    parser.add_argument('sample', metavar='sample_name', help='sample_name to title every thing.')
    # options
    parser.add_argument('-A', dest='annotation', default='.',
                        help='if specified, FUCHS will use the annotation to reconstruct the internal structure of unsupported circRNA regions')
    parser.add_argument('-c', dest='num_cpus', required=False, type=int, default=multiprocessing.cpu_count(),
                        help='Number of processors to use')
    parser.add_argument('--tmp', dest='tmp_folder', default='.',
                        help='tempfolder to store tempfiles generated by pybedtools.')

    # parse arguments
    args = parser.parse_args()

    infolder = args.inputfolder
    sample = args.sample
    annotation_file = args.annotation
    num_cpus = args.num_cpus
    tmp_folder = args.tmp_folder

    # set temp folder. foldr needs to exist!
    tempfile.tempdir = tmp_folder
    pybedtools.set_tempdir(tmp_folder)

    folder = '%s/%s/' % (infolder, sample)
    outfile = '%s/%s_exon_chain_' % (infolder, sample)

    output = open('%sinferred_12.bed' % (outfile), 'w')
    output.write('#bed12\n')
    output.close()

    output = open('%sinferred_6.bed' % (outfile), 'w')
    output.write('#bed6\n')
    output.close()

    # Start my pool
    pool = multiprocessing.Pool(num_cpus)

    files = os.listdir(folder)
    files = sorted(files)
    # Build task list
    tasks = []
    plotNum = 0
    for f in files:
        if f.split('.')[-1] == 'bam':
            tasks.append((f, folder, annotation_file, outfile))

    print(len(files))
    print("Processing %d circRNAs using %d processors..." % (len(tasks), num_cpus))

    # Run tasks
    results = [pool.apply_async(run_denovo_exon_chain_reconstruction, t) for t in tasks]

    # Process results
    for result in results:
        (filename, introns) = result.get()
        print("Result: circRNA %s has %s introns" % (filename, introns))

    pool.close()
    pool.join()
    pybedtools.helpers.cleanup()
