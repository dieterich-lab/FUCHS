*****************************************
FUCHS - FUll circle CHaracterization from rna-Seq
*****************************************
FUCHS is a python pipeline desigend to fully characterize circular RNAs. It uses a list of circular RNAs and reads spanning the back-splice junction as well 
as a BAM file containing the mapping of all reads (alternatively of all chimeric reads).

The reads from one circle are extracted by FUCHS and saved in an individual BAM file. Based on these BAM files, FUCHS will detect alternative splicing within the same 
circle boundaries, summarize different circular isoforms from the same host-gene and generates coverage plots for each circle. It will also cluster circles based on their
coverage profile. These results can be used to identify potential false positive circles.

=============
Installation
=============

FUCHS dependes on pysam, pybedtools, numpy and R (amap, Hmisc, gplots).

Right now FUCHS is an assembly of scripts which is called through one main scripts executing all steps. No installation is neccessary right now!

You can clone the scripts from git using git clone:

.. code-block:: bash

  $ git clone git@github.com:dieterich-lab/FUCHS.git
  
  $ cd FUCHS/scripts
  
  $ python FUCHS.py [options] <sample_circleIDs> <sample.bam> <feature.bed> <PATH/to/output/folder> <sample_name>

========
Usage
========
To characterize circRNAs from RNA-seq data you have to:

1. Map RNAseq data from quality checked fastq files with either STAR, BWA, TopHat-Fusion.

2. Detect circRNAs using DCC, CIRI, CIRCfinder or CIRCexplorer depending on the program you used for mapping.

3. Run FUCHS (right now only the combination STAR-DCC has been tested, the rest is under development.)

========================
Step by step tutorial
========================
In this tutorial I will be using HEK293 data (published with the paper?) and use STAR with DCC to detect circular RNAs

1. Map RNA-seq data with STAR (Dobin et al., 2013). Note that --alignSJoverhangMin and --chimJunctionOverhangMin should use the same value, to make the circRNA expression and linear gene expression level comparable. 
Note that STARlong is not mapping chimeric reads correctly. 

* Do pairs joined mapping first. If your data are paired end, do additional mates separate mapping (not mandatory, but will increase the sensitivity of DCC detection, because it collect small circRNAs which appear with one chimeric junction point at each read mate). If the data is single end, only one mapping step is needed. In this case, we have PE sequencing data.

.. code-block:: bash

  $ mkdir Sample1
  $ cd Sample1
  $ STAR --runThreadN 10   --genomeDir [genome]  --outSAMtype BAM Unsorted --readFilesIn Sample1_1.fastq.gz  Sample1_2.fastq.gz   --readFilesCommand zcat  --outFileNamePrefix [sample prefix] --outReadsUnmapped Fastx  --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --outFilterMultimapNmax 20   --outFilterScoreMin 1   --outFilterMatchNmin 1   --outFilterMismatchNmax 2  --chimSegmentMin 15    --chimScoreMin 15   --chimScoreSeparation 10  --chimJunctionOverhangMin 15


* (Skip when you have single end data). Mates separate mapping. Be careful that, what you define as first mate (mate1) should also appears the first in the joined mapping. In this case, SamplePairedRead_1.fastq.gz is the first mate which came first above.

.. code-block:: bash

  # Make a directory for mate1
  $ mkdir mate1
  $ STAR --runThreadN 10   --genomeDir [genome]  --outSAMtype None --readFilesIn Sample1_1.fastq.gz  --readFilesCommand zcat   --outFileNamePrefix [sample prefix] --outReadsUnmapped Fastx  --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30  --outFilterMultimapNmax 20   --outFilterScoreMin 1   --outFilterMatchNmin 1   --outFilterMismatchNmax 2  --chimSegmentMin 15    --chimScoreMin 15   --chimScoreSeparation 10  --chimJunctionOverhangMin 15

  $ cd ..
  $ mkdir mate2
  # Do the same mapping as mate1 for mate2

2. Detect circRNAs from chimeric.out.junction files with DCC

- It is strongly recommended to specify a repetitive region file in GTF format for filtering. You can obtain this file through UCSC table browser: http://genome.ucsc.edu/cgi-bin/hgTables. Select your genome, select group as "Repeats" or "Variation and Repeats". For the track, I recommend chose all possible repeats and combine the results. **NOTE**: the output file needs to comply with GTF format specification.

- Prepare path files to specify where is your chimeric.junction.out files are. 

  First, "samplesheet" file, in which you specify your chimeric.out.junction file's absolute paths (mates joined mapping chimeric.out.junction files, for paired end data), one line per sample. 

  Second (only if you have paired end sequencing data), "mate1" and "mate2" files. As with the "samplesheet" file, you specify where your mate1 and mate2 separately mapped chimeric.junction.out files are.

  You can find a example of these files for HEK293 data at:
  
.. code-block:: bash

  $ <FUCHS directory>/testdata/dcc/samplesheet # Mates jointly mapped chimeric.junction.out files
  $ <FUCHS directory>/testdata/dcc/mate1 # Mate1 independently mapped chimeric.junction.out files
  $ <FUCHS directory>/testdata/dcc/mate1 # Mate2 independently mapped chimeric.junction.out files

- After all the preparation steps, you can now run DCC for circRNA detection. 


.. code-block:: bash

  # Call DCC to detect circRNAs, using HEK293 data as example.
  $ DCC @samplesheet -mt1 @mate1 -mt2 @mate2 -D -R [Repeats].gtf -an [Annotation].gtf -Pi -F -M -Nr 5 6 -fg -G -A [Reference].fa

  # Details of parameters please refer to the help page of DCC:
  $ DCC -h

By default, DCC assume the data are stranded, for non-stranded data, use -N flag.
NOTE: -F flag is mandatory, if you want to filter on the results. All filtering steps are not mandatory, but strongly recommended.

--------------------

The output of DCC include: CircRNACount, CircCoordinates, LinearCount and CircSkipJunctions.

**CircRNACount:** a table containing read counts for circRNAs detected. First three columns are chr, circRNA start, circRNA end. From fourth column on are the circRNA read counts, one sample per column, shown in the order given in your samplesheet.

**CircCoordinates:** CircRNA annotation in BED format. The columns are chr, start, end, genename, junctiontype (come from STAR, 1 for GT/AG, 2 for CT/AC), strand, circRNA region (startregion-endregion), overall regions (the genomic features circRNA coordinates interval covers).

**LinearCount:** host gene expression count table, same setup with CircRNACount file.

**CircSkipJunctions:** CircSkip junctions. First three columns are the same with LinearCount/CircRNACount, the rest columns are circSkip junctions found for each sample. circSkip junction shows in the format: chr:start-end:count (chr1:1787-6949:10 for example. It's possible that for one circRNA multiple circSkip junctions are found, because circRNA possible come from multiple RNA isoforms. In this case, multiple circSkip junctions are delimited with semicolon). 0 implies no circSkip junction found for this circRNA.

-----------------------------------

3. Merge mate1.chimeric.sam and mate2.chimeric.sam files for FUCHS (This is not neccessary if circles were detected using BWA/CIRI)

.. code-block:: bash

  $ samtools view -Sb -o hek293.1 hek293.1/Chimeric.out.sam
  $ samtools view -Sb -o hek293.2 hek293.2/Chimeric.out.sam
  
  $ samtools sort hek293.1 hek293.1.sorted
  $ samtools sort hek293.2 hek293.2.sorted
   
  $ samtools index hek293.1.sorted.bam
  $ samtools index hek293.2.sorted.bam
   
  $samtools merge hek293.sorted.bam hek293.1.sorted.bam hek293.2.sorted.bam
   
  $samtools index hek293.sorted.bam

4. Run FUCHS.py to start the pipeline which will extract reads, check mate status, detect alternative splicing events, classify different isoforms, generate coverage profiles and cluster circRNAs based on coverage profiles

.. code-block:: bash
  $ python FUCHS.py -r 2 -q 2 -p refseq -e 3 -c CircRNACount -m hek293.mate1.Chimeric.out.junction.fixed -j hek293.mate2.Chimeric.out.junction.fixed mock hek293.sorted.bam hg38.refseq.bed FUCHS/ hek293
  
  # if you used BWA/CIRI you can skip -c, -m, and -j, specify to skip the first step -sS step1 and specify the circIDs file

**Finished!!!**

========================================================================
INPUT
========================================================================
**circIDs:** 
.. code-block:: bash
  $ 1:3740233|3746181	MISEQ:136:000000000-ACBC6:1:2107:10994:20458,MISEQ:136:000000000-ACBC6:1:1116:13529:8356
  $ 1:8495063|8557523	MISEQ:136:000000000-ACBC6:1:2117:11302:22227,MISEQ:136:000000000-ACBC6:1:1111:4979:10994,MISEQ:136:000000000-ACBC6:1:2117:14163:16664,MISEQ:136:000000000-ACBC6:1:1103:13343:14303
  $ 1:8495063|8614686	MISEQ:136:000000000-ACBC6:1:2118:9328:9926

The first column contains the circle id formated as folllowed **chr:start|end**. The second column is a comma separated list of read names spanning the back-splice junction.

**bamfile:** Alignment file produced by any mapper. This file must contain all chimerically mapped reads and may contain also linearly mapped reads.

**bedfile:** 
.. code-block:: bash
  $ 1	67092175	67093604	NR_075077_exon_0_0_chr1_67092176_r	0	-
  $ 1	67096251	67096321	NR_075077_exon_1_0_chr1_67096252_r	0	-
  $ 1	67103237	67103382	NR_075077_exon_2_0_chr1_67103238_r	0	-

Normal BED file in BED6 format. The name should contain a gene name or gene ID and the exon_number. You can specify how the name should be processed using -p (platform), -s (character used to separate name and exon number) and -e (exon_index). 

========================================================================
OUTPUT
========================================================================
**hek293.alternative_splicing.txt:** 
**hek293.exon_counts.bed:** 
**hek293.exon_counts.txt:** 
**hek293.mate_status.txt:** 
**hek293.skipped_exons.bed:** 
**hek293.skipped_exons.txt:** 
--------------------

hek293:
20_41533050_41551360_5reads.sorted.bam
--------------------

hek293.coverage_pictures:
20_41533050_41551360_NM_032221.png
cluster_means_all_circles.png
--------------------

hek293.coverage_profiles:
20_41533050_41551360.NM_032221.txt
coverage.clusters.all_circles.pdf
coverage_profiles.all_circles.pdf

--------------------


========================
Error and solutions
========================

 - ERROR: File <file> has inconsistent naming convention for record:
   CHR_MG132_PATCH 124291803 124294101 ENSMUSG00000098810 . - protein_coding exon CAAA01180111.1
Please update your bedtools at least to 2.24.0, and make sure the new version is included in your path.
