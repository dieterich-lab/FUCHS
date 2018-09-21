#! /usr/bin/env python3

# main script to run FUCHS


def main():
    # required packages
    import os
    import argparse
    import datetime
    import time

    parser = argparse.ArgumentParser(description='Main script of the FUCHS pipeline.'
                                                 'For a detailed help see https://github.com/dieterich-lab/FUCHS '
                                                 'or the included README.rst file.')

    # input
    parser.add_argument('-C', '--circIDs', dest='circlefile', default='none',
                        help='Tab-separated file chr:start_end(tab)read1,read2,read3.')
    parser.add_argument('-D', '--DCC', dest='CircRNACount', default='none',
                        help='If you mapped with STAR and are using step1 you need to provide a list'
                             ' of circle ids (CircRNACount or CircCoordinates from DCC)'
                             'You must supply either -C or -DCC')
    parser.add_argument('-J', '--chimericJunctions', dest='chimeric_junction', default='none',
                        help='If you mapped with STAR and are using step1 you need to provide the paired end Chimeric.junction.out file here')
    parser.add_argument('-F', '--mate1', dest='mate1', default='none',
                        help='If you mapped with STAR and are using step1 you need to provide the mate1.Chimeric.junction.out file here (optional if ends were mapped separately)')
    parser.add_argument('-R', '--mate2', dest='mate2', default='none',
                        help='If you mapped with STAR and are using step1 you need to provide the mate2.Chimeric.junction.out file here (optional if ends were mapped separately)')
    parser.add_argument('-B', '--bamfile', dest='bamfile', required=True,
                        help='BAM file containing chimeric reads, linear reads may be in it but are not required.')
    parser.add_argument('-A', '--annotation', dest='bedfile', required=True,
                        help='bed formatted feature file including exons.')
    # output
    parser.add_argument('-O', '--outFolder', dest='out_folder', default='.',
                        help='Output folder. There will be a sub folder for the sample containing a BAM file '
                             'for each circle.')
    parser.add_argument('-N', '--sampleName', dest='sample', required=True,
                        help='sample name to title every thing.')

    # options
    parser.add_argument('-r', '--thresholdReads', dest='reads', default=5, type=int,
                        help='Circle has to have at least <r> reads to be analysed.')

    # TODO: default: no multi map
    parser.add_argument('-q', '--thresholdMapq', dest='mapq', default=3, type=int,
                        help='MAPQ cutoff, only reads passing this threshold will be written to circle BAM file.')
    # TODO: add 0 based info
    parser.add_argument('-c', '--splitCharacter', dest='split_character', default='_', help='feature name separator.')
    parser.add_argument('-e', '--exonIndex', dest='exon_index', default=3, type=int,
                        help='Field indicating the exon number after splitting feature name by split_character (for the annotation file).')
    parser.add_argument('-p', '--annotationFormat', dest='ref_platform', default='refseq',
                        help='Specifies the annotation platform which was used (refseq or ensembl)')
    parser.add_argument('-s', '--skipSteps', dest='skipped_steps', default='none',
                        help='Comma separated list of steps that should be skipped (e.g. step3,step4,step6)')
    parser.add_argument('-T', '--tmp', dest='tmp_folder', default='/tmp/',
                        help='Folder to store temporary files generated by pybedtools.')

    parser.add_argument('-P', '--cpus', dest='num_cpus', default=4, type=int,
                        help='Number of CPUs used.')

    args = parser.parse_args()

    # parse arguments
    circles = os.path.expanduser(args.circlefile)
    circle_ids = os.path.expanduser(args.CircRNACount)
    paired = os.path.expanduser(args.chimeric_junction)  # not the greatest naming scheme
    mate1 = os.path.expanduser(args.mate1)
    mate2 = os.path.expanduser(args.mate2)
    bamfile = os.path.expanduser(args.bamfile)
    bedfile = os.path.expanduser(args.bedfile)
    outfolder = os.path.expanduser(args.out_folder) + '/'
    sample = args.sample
    num_cpus = args.num_cpus

    cutoff_reads = args.reads
    cutoff_mapq = args.mapq
    exon_index = args.exon_index
    split_character = args.split_character
    platform = args.ref_platform
    skipped_steps = args.skipped_steps.split(',')
    tmp_folder = os.path.expanduser(args.tmp_folder) + '/'

    # start writing down FUCHS time for retracing
    print(('Started FUCHS at %s' % (datetime.datetime.now())))
    dt = str(datetime.datetime.now())
    start_time = time.time()
    # make log file
    # TODO
    # test if command line was correct
    if circles == 'none' and circle_ids == 'none':
        print(
            'ERROR, you need to specify either a -C or -DCC.\nIf you mapped and detected your circRNAs with STAR/DCC you may indicate \n-DCC CircRNACount, -CJ Chimeric.junction.out, -m1 mate1.Chimeric.junction.out and -m2 mate2.Chimeric.junction.out\nif you used a different program, please supply a circID list using -C.\n')
        quit()

    if not circles == 'none' and not circle_ids == 'none':
        print(
            'You have indicated both -C and -DCC, this is not necessary, we will skip step1 (read extraction from the STAR output) and proceed with the circID file\n')
        circle_ids == 'none'
        skipped_steps += ['step1']

    if not circle_ids == 'none' and paired == 'none':
        print(
            'You have indicated that you detected your circRNAs using STAR/DCC with the -DCC flag, \nhowever you did not specify a Chimeric.junction.out file, this is necessary, \nplease specify at least -CJ, if you have paired end data also specify -m1/-m2\n')
        quit()

    # convert relative paths names to absolute path names
    working_dir = os.getcwd()
    if not circles == 'none' and not os.path.isabs(circles):
        circles = os.path.abspath(os.path.join(os.getcwd(), circles))
        print(('changed circID file to %s\n' % (circles)))
    if not circles == 'none' and not os.path.exists(circles):
        print(('ERROR, no such file or directory: %s' % (circles)))
        quit()

    if not circle_ids == 'none' and not os.path.isabs(circle_ids):
        circle_ids = os.path.abspath(os.path.join(os.getcwd(), circle_ids))
        print(('changed CircRNACount file to %s\n' % (circle_ids)))
    if not circle_ids == 'none' and not os.path.exists(circle_ids):
        print(('ERROR, no such file or directory: %s' % (circle_ids)))
        quit()

    if not paired == 'none' and not os.path.isabs(paired):
        paired = os.path.abspath(os.path.join(os.getcwd(), paired))
        print(('changed Chimeric.junction.out file to %s\n' % (paired)))
    if not paired == 'none' and not os.path.exists(paired):
        print(('ERROR, no such file or directory: %s' % (paired)))
        quit()

    if not mate2 == 'none' and not os.path.isabs(mate2):
        mate2 = os.path.abspath(os.path.join(os.getcwd(), mate2))
        print(('changed mate2.Chimeric.junction.out file to %s\n' % (mate2)))
    if not mate2 == 'none' and not os.path.exists(mate2):
        print(('ERROR, no such file or directory: %s' % (mate2)))
        quit()

    if not mate1 == 'none' and not os.path.isabs(mate1):
        mate1 = os.path.abspath(os.path.join(os.getcwd(), mate1))
        print(('changed mate1.Chimeric.junction.out file to %s\n' % (mate1)))
    if not mate1 == 'none' and not os.path.exists(mate1):
        print(('ERROR, no such file or directory: %s' % (mate1)))
        quit()

    if not os.path.isabs(bamfile):
        bamfile = os.path.abspath(os.path.join(os.getcwd(), bamfile))
        print(('changed bamfile file to %s\n' % (bamfile)))
    if not os.path.exists(bamfile):
        print(('ERROR, no such file or directory: %s' % (bamfile)))
        quit()

    if not os.path.isabs(outfolder):
        outfolder = os.path.abspath(os.path.join(os.getcwd(), outfolder))
        print(('changed output folder to %s\n' % (outfolder)))
    if not os.path.isdir(outfolder):
        os.mkdir(outfolder)

    if not os.path.isabs(tmp_folder):
        tmp_folder = os.path.abspath(os.path.join(os.getcwd(), tmp_folder))
        print(('changed tmp folder to %s\n' % (tmp_folder)))
    if not os.path.isdir(tmp_folder):
        os.mkdir(tmp_folder)

    if not os.path.isabs(bedfile):
        bedfile = os.path.abspath(os.path.join(os.getcwd(), bedfile))
        print(('changed bedfile file to %s\n' % (bedfile)))
    if not os.path.exists(bedfile):
        print(('ERROR, no such file or directory: %s' % (bedfile)))
        quit()

    accepted_platforms = ('refseq', 'ensembl')
    platform = platform.lower()
    if not platform in accepted_platforms:
        print('ERROR please specify an accepted annotation platform. Possible options are: refseq or ensembl')
        quit()

    print("The following analysis steps will be skipped: " + '%s' % ', '.join(map(str, skipped_steps)))

    # Step 1: (optional) if DCC was used, extract circle read names from junction file    
    output_file = open('%s/%s.logfile.%s' % (outfolder, sample, dt.replace(' ', '_')), 'w')
    output_file.write('FUCHS is starting at %s\n\n' % (dt))
    output_file.write('%s: starting to get readnames from Chimeric.junction.out\n' % (datetime.datetime.now()))
    output_file.close()
    if not 'step1' in skipped_steps:

        circles = "%s.reads.txt" % paired
        if not os.path.isfile(circles):
            import get_readnames_from_DCC as get_readnames
            names = get_readnames.get_readnames_from_DCC(circle_ids, paired, mate1, mate2)
            names.run()
        else:
            output_file = open('%s/%s.logfile.%s' % (outfolder, sample, dt.replace(' ', '_')), 'a')
            output_file.write('\tskipping get_readnames_from_DCC because %s exists already\n' % (circles))
            output_file.close()

    # Step2 : extract circle reads from sample bam file
    output_file = open('%s/%s.logfile.%s' % (outfolder, sample, dt.replace(' ', '_')), 'a')
    output_file.write('\tfinished\n\n%s: starting to extract chimeric reads from bamfile\n' % (datetime.datetime.now()))
    output_file.close()
    if not 'step2' in skipped_steps:
        import extract_reads as extract_reads
        er = extract_reads.extract_reads(cutoff_reads, cutoff_mapq, circles, bamfile, outfolder, sample, tmp_folder,
                                         num_cpus)
        er.run()

    # Step3 : (optional) get information about possibly rolling circles 
    output_file = open('%s/%s.logfile.%s' % (outfolder, sample, dt.replace(' ', '_')), 'a')
    output_file.write('\tfinished\n\n%s: starting to get mate pair information\n' % (datetime.datetime.now()))
    output_file.close()
    if not 'step3' in skipped_steps:

        if not os.path.isfile('%s/%s.mate_status.txt' % (outfolder, sample)):
            import get_mate_information as mateinformation
            mi = mateinformation.mate_information(platform, split_character, bedfile, outfolder, sample, tmp_folder,
                                                  num_cpus)
            mi.run()
        else:
            output_file = open('%s/%s.logfile.%s' % (outfolder, sample, dt.replace(' ', '_')), 'a')
            output_file.write(
                '\tskipping get_mate_information because %s/%s.mate_status.txt exists already\n' % (outfolder, sample))
            output_file.close()

    # Step4 : (optional) find exon skipping events
    output_file = open('%s/%s.logfile.%s' % (outfolder, sample, dt.replace(' ', '_')), 'a')
    output_file.write('\tfinished\n\n%s: starting to detect skipped exons\n' % (datetime.datetime.now()))
    output_file.close()
    if not 'step4' in skipped_steps:

        if not os.path.isfile('%s/%s.skipped_exons.bed' % (outfolder, sample)):
            import detect_skipped_exons as skipped_exons
            se = skipped_exons.detect_skipped_exons(outfolder, sample, bedfile, tmp_folder, platform, num_cpus)
            se.run()
        else:
            output_file = open('%s/%s.logfile.%s' % (outfolder, sample, dt.replace(' ', '_')), 'a')
            output_file.write('\tskipping detect_skipped_exons because %s/%s.skipped_exons.bed exists already\n' % (
                outfolder, sample))
            output_file.close()

    # Step5 : (optional) identify different circles within the same host gene
    output_file = open('%s/%s.logfile.%s' % (outfolder, sample, dt.replace(' ', '_')), 'a')
    output_file.write('\tfinished\n\n%s: starting to detect alternative splicing\n' % (datetime.datetime.now()))
    output_file.close()
    if not 'step5' in skipped_steps:

        if not os.path.isfile('%s/%s.alternative_splicing.txt' % (outfolder, sample)):
            import detect_splicing_variants as splicing_variants
            sv = splicing_variants.detect_splicing_variants(split_character, platform, circles, bedfile,
                                                            outfolder, sample, tmp_folder, num_cpus)
            sv.run()
        else:
            output_file = open('%s/%s.logfile.%s' % (outfolder, sample, dt.replace(' ', '_')), 'a')
            output_file.write(
                '\tskipping detect_splicing_variants because %s/%s.alternative_splicing.txt exists already\n' % (
                    outfolder, sample))
            output_file.close()

    # Step6 : (optional) generate coverage profile for each circle
    # (one transcript per gene, best if most fitting transcript)
    output_file = open('%s/%s.logfile.%s' % (outfolder, sample, dt.replace(' ', '_')), 'a')
    output_file.write('\tfinished\n\n%s: starting to generate coverage profiles\n' % (datetime.datetime.now()))
    output_file.close()
    if not 'step6' in skipped_steps:
        if not os.path.isfile('%s/%s.exon_counts.bed' % (outfolder, sample)) and not os.path.isdir(
                        '%s/%s.coverage_profiles/' % (outfolder, sample)):
            import get_coverage_profile as coverage_profile
            sv = coverage_profile.get_coverage_profile(exon_index, split_character, platform, bedfile,
                                                       outfolder, sample, tmp_folder, num_cpus)
            sv.run()
        else:
            output_file = open('%s/%s.logfile.%s' % (outfolder, sample, dt.replace(' ', '_')), 'a')
            output_file.write(
                '\tskipping get_coverage_profile because %s/%s.exon_counts.bed exists already\n' % (outfolder, sample))
            output_file.close()

    # Step7 : (optional, requires step 5)
    output_file = open('%s/%s.logfile.%s' % (outfolder, sample, dt.replace(' ', '_')), 'a')
    output_file.write('\tfinished\n\n%s: starting to summarize the coverage profiles\n' % (datetime.datetime.now()))
    output_file.close()
    if not 'step7' in skipped_steps:
        if not os.path.isfile('%s/%s.coverage_profiles/coverage_profiles.all_circles.pdf' % (outfolder, sample)):
            if os.path.isdir('%s/%s.coverage_profiles/' % (outfolder, sample)):
                os.system('summarized_coverage_profiles.R %s/%s.coverage_profiles' % (outfolder, sample))
            else:
                output_file = open('%s/%s.logfile.%s' % (outfolder, sample, dt.replace(' ', '_')), 'a')
                output_file.write('\tYou are trying cluster the coverage profiles without '
                                  'generating coverage profiles first, please run step 5 (get_coverage_profile)\n')
                output_file.close()
        else:
            output_file = open('%s/%s.logfile.%s' % (outfolder, sample, dt.replace(' ', '_')), 'a')
            output_file.write(
                '\tskipping summarized_coverage_profiles.R because %s/%s.coverage_profiles/coverage_profiles.all_circles.pdf already exists\n' % (
                    outfolder, sample))
            output_file.close()

    # Step8 : (optional, requires step6) pictures for all circles
    output_file = open('%s/%s.logfile.%s' % (outfolder, sample, dt.replace(' ', '_')), 'a')
    output_file.write('\tfinished\n\n%s: starting to visualize the coverage profiles\n' % (datetime.datetime.now()))
    output_file.close()
    if not 'step8' in skipped_steps:
        if os.path.isdir('%s/%s.coverage_profiles/' % (outfolder, sample)):
            files = os.listdir('%s/%s.coverage_profiles' % (outfolder, sample))
            folders = os.listdir(outfolder)
            if not '%s.coverage_pictures' % (sample) in folders:
                os.mkdir('%s/%s.coverage_pictures' % (outfolder, sample))

            def run_r_parallel(f):
                    if f.endswith('.txt'):
                        os.system('make_coverage_picture.R %s/%s.coverage_profiles/%s %s/%s.coverage_pictures/' %
                                  (outfolder, sample, f, outfolder, sample))

            from pathos.multiprocessing import ProcessingPool as Pool

            pool = Pool(num_cpus)
            pool.map(run_r_parallel, files)

        else:
            output_file = open('%s/%s.logfile.%s' % (outfolder, sample, dt.replace(' ', '_')), 'a')
            output_file.write('\tYou are trying to generate coverage pictures '
                              'without generating coverage profiles, please run step 5 (get_coverage_profile)\n')
            output_file.close()

    output_file = open('%s/%s.logfile.%s' % (outfolder, sample, dt.replace(' ', '_')), 'a')
    output_file.write('\tfinished\n\n\nFUCHS finished at %s\n\n' % (datetime.datetime.now()))
    output_file.write("FUCHS took --- %s minutes ---\n\n" % (round((time.time() - start_time) / 60.0)))
    output_file.close()


if __name__ == '__main__':
    main()
