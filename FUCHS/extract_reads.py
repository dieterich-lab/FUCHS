# python script, to extract circular reads from a bam file, based on a
# tab-separated circle_id - reads file and a bam file. writes out circle-specific bam files

# required packages
import os
import pysam
import tempfile


def run_parallel(f):
    # no do re-sort sorted files and create duplicates
    if f.split('.')[-1] == 'bam' and "sorted" not in f:
        pysam.sort("-o",
                   '%s' % (f.replace('.bam', '.sorted.bam')),
                   '%s' % (f)
                   )

        pysam.index('%s' % (f.replace('.bam', '.sorted.bam')))
        os.system('rm %s' % (f))


class extract_reads(object):
    def __init__(self, reads, mapq, circlefile, bamfile, outputfolder, sample, tmp_folder, cpus):

        self.circles = circlefile
        self.bamfile = bamfile
        self.outfolder = outputfolder
        self.sample = sample
        self.cutoff = reads
        self.mapq_cutoff = mapq
        self.tmp_folder = tmp_folder
        self.cpus = cpus

    def read_circles(self, infile):
        """
        for each circle, extracts all reads_ids, but also accumulates all circular reads ids independent of circle id.
        """
        circle_IDs = {}
        reads = {}
        input_file = open(infile)
        for line in input_file:
            if not line.startswith('#'):
                circle_IDs[line.split('\t')[0]] = line.replace('\n', '').split('\t')[1].split(',')[:-1]
                for current_circle in line.replace('\n', '').split('\t')[1].split(',')[:-1]:
                    reads[current_circle] = 0
        input_file.close()
        return circle_IDs, reads

    def load_alignment(self, infile, circle_reads, cutoff):
        """
        loads sample bam file and extracts all circular reads based
        on read list obtained from read_circles or any given list of read ids.
        """
        reads = {}
        if 'bam' in infile.split('.')[-1]:
            samfile = pysam.AlignmentFile(self.bamfile, 'rb')
        else:
            samfile = pysam.AlignmentFile(self.bamfile, 'r')
        for i, read in enumerate(samfile.fetch()):
            if read.query_name in circle_reads and read.mapq > cutoff and read.get_tag("HI") == 1:
                if not read.query_name in reads:
                    reads[read.query_name] = {}
                reads[read.query_name][(read.reference_start, read.cigarstring, read.is_reverse)] = read
        samfile.close()
        return reads

    def write_circle_bam(self, reads, circles, cutoff, template, outfolder):
        """
        for each circle, writes a bam file containing only reads spanning the
        circle junction and their mates if mates are present.
        """
        samfile = pysam.AlignmentFile(template, 'rb')
        for circle in circles:
            if len(circles[circle]) >= cutoff:
                circle_bam = pysam.AlignmentFile(
                    "%s/%s_%sreads.bam" % (outfolder, circle.replace(':', '_').replace('|', '_'), len(circles[circle])),
                    "wb", template=samfile)
                has_content = 0
                for read in circles[circle]:
                    if read in reads:
                        for part in reads[read]:
                            circle_bam.write(reads[read][part])
                            has_content += 1
                circle_bam.close()
                # nothing was written to the bam file, delete it so that following steps don't operate
                # on an empty file
                if has_content == 0:
                    os.remove("%s/%s_%sreads.bam" % (
                        outfolder, circle.replace(':', '_').replace('|', '_'), len(circles[circle])))

        samfile.close()
        return

    def run(self):

        tempfile.tempdir = self.tmp_folder  # set global tmp dir

        circle_info, circle_reads = self.read_circles(self.circles)
        print(('DONE reading circles, found %s circles' % (len(circle_info))))
        reads = self.load_alignment(self.bamfile, circle_reads, self.mapq_cutoff)
        print('DONE extracting circular reads')
        folders = os.listdir(self.outfolder)
        if not self.sample in folders:
            os.mkdir('%s/%s' % (self.outfolder, self.sample))
        self.write_circle_bam(reads, circle_info, self.cutoff, self.bamfile, '%s/%s' % (self.outfolder, self.sample))
        print('DONE writing circle bam files\n')
        # files = os.listdir('%s/%s' % (self.outfolder, self.sample))
        import glob

        files = glob.glob('%s/%s/*.bam' % (self.outfolder, self.sample))

        # possible sorted files from previous run
        sorted_bams = glob.glob('%s/%s/*.sorted.bam' % (self.outfolder, self.sample))

        # fix the file / circle count
        actual_bams = len(files) - len(sorted_bams)

        print(('%s circles passed your thresholds of at least %s reads with at least a mapq of %s\n\n' % (
            actual_bams, self.cutoff, self.mapq_cutoff)))

        from pathos.multiprocessing import ProcessingPool as Pool

        pool = Pool(self.cpus)
        pool.map(run_parallel, files)
