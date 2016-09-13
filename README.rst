*****************************************
FUCHS - FUll circular RNA CHaracterization from rna-Seq
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

==================== ==========================================================================================
 circID               read1,read2,read3                                                                        
==================== ==========================================================================================
 1:3740233\|3746181  MISEQ:136:000000000-ACBC6:1:2107:10994:20458,MISEQ:136:000000000-ACBC6:1:1116:13529:8356 
 1:8495063\|8614686  MISEQ:136:000000000-ACBC6:1:2118:9328:9926                                               
==================== ==========================================================================================


The first column contains the circle id formated as folllowed **chr:start|end**. The second column is a comma separated list of read names spanning the back-splice junction.

**bamfile:** Alignment file produced by any mapper. This file must contain all chimerically mapped reads and may contain also linearly mapped reads.

**bedfile:** 

====   ===========    =============     ===================================   =======  ======
Chr      Start            End               Name                               Score   Strand
====   ===========    =============     ===================================   =======  ======
 1      67092175        67093604         NR_075077_exon_0_0_chr1_67092176_r     0       \-
 1      67096251        67096321         NR_075077_exon_1_0_chr1_67096252_r     0       \-
 1      67103237        67103382         NR_075077_exon_2_0_chr1_67103238_r     0       \-
====   ===========    =============     ===================================   =======  ======

Normal BED file in BED6 format. The name should contain a gene name or gene ID and the exon_number. You can specify how the name should be processed using -p (platform), -s (character used to separate name and exon number) and -e (exon_index). 

========================================================================
OUTPUT
========================================================================

**hek293.alternative_splicing.txt:** 

This file summarizes the relationship of different circRNAs derived from the same host-gene. 

=============  ============================================================    =========================================  =========   ===========  =============================================
Transcript      circles                                                        same_start                                 same_end    overlapping  within
=============  ============================================================    =========================================  =========   ===========  =============================================
NM_016287	1:20749723-20773610                                            .                                           .          .            .
NM_005095	1:35358925-35361789,1:35381259-35389082,1:35381259-35390098    1:35381259-35389082|1:35381259-35390098,    .          .            .
NM_001291940    1:236803428-236838599,1:236806144-236816543                    .                                           .          .            1:236803428-236838599|1:236806144-236816543,
=============  ============================================================    =========================================  =========   ===========  =============================================

| *Transcript*: Transcript name as defined by the bed-annotation file
| *circles*: Comma-separated list of circRNA ids derived from this transcript
| *same_start*: Comma-seprated list of circRNA pairs separated by |. Pairs in this column share the same start coordinates. A "." indicates that there are no circle pairs that share the same start coordinates.
| *same_end*: Same as *same_start*, only now, circle pairs share the same end coordinates.
| *overlapping*: Comma-seprated list of circRNA pairs separated by |. Pairs in this column share neither start nor end coordinates, but their relation is such that: start.x < start.y && end.x < end.y && start.y < end.x
| *within*: Same as *overlapping*, only now, circle pairs have the follwoing relation: start.x < start.y && end.x > end.y
| 

**hek293.exon_counts.bed:** 
This file is a bed-formatted file that describes the exon-structure and can be loaded into any genome browser. Each line corresponds to a circRNA.

=====  ============  =============    ============    =============    =======   ======== =========  ======= ===========  ==============  =====================
Chr    Circle Start   Circle  End      Transcript     Num of Reads     Strand      Start   End        Color  Num of Exon  Exon Lengths     Relative Exon Starts   
=====  ============  =============    ============    =============    =======   ======== =========  ======= ===========  ==============  =====================
chr1    35358925        35361789        NM_005095       9               \+       35358925 35361789   0,255,0  3           521,61,170      0,2269,2694
chr1    20749723        20773610        NM_016287       4               \-       20749723 20773610   0,255,0  4           159,90,143,159  0,7443,21207,23728
=====  ============  =============    ============    =============    =======   ======== =========  ======= ===========  ==============  =====================

| *Chr*: Chromosome of circRNA
| *Circle Start*: The 5' site of the chimeric junction. This is relative to the reference strand, i.e. start < end! The location is 1-index based
| *Cirlce End*: The 3' site of the chimeric junction. This is relative to the reference strand, i.e. start < end! The location is 0-index based
| *Transcript*: Transcript name as defined by the bed-annotation file
| *Num of Reads* : Number of reads supporting this chimeric junction, in other words, reads that are chimerically mapped to this junction
| *Strand*: Strand of the host-gene
| *Start*: Copied *Circle Start* to stay conform with BED12 format
| *End*: Copied *Circle End* to stay conform with BED12 format
| *Color*: pre defined color the exons will show up in the genome viewer (0,255,0 -> green)
| *Num of Exon*: Number of exons in this circRNA consists of
| *Exon Lengths*: Comma-seprated list of the length of each exon
| *Relative Exon Starts*: Comma-separated list of the relative starting positions of the exons within the circle boundaries.
| 
**hek293.exon_counts.txt:** 
This file contains similar information as the previous file, just more detailed inforamtion on the exons. Each line corresponds to one exon.

======= =====================  ================ ============  ========== =====  ============   ============= ======= =============   ==============  ===========     ========= ========
sample   circle_id               transcript_id   other_ids       exon_id chr     start           end          strand  exon_length     unique_reads    fragments       number\+ number\-
======= =====================  ================ ============  ========== =====  ============   ============= ======= =============   ==============  ===========     ========= ========
hek293   1:35358925-35361789     NM_005095       NM_005095       2       1       35358924        35359446        \+       522          9               9               4        5
hek293   1:35358925-35361789     NM_005095       NM_005095       3       1       35361193        35361255        \+       62           3               3               1        2
hek293   1:35358925-35361789     NM_005095       NM_005095       4       1       35361618        35361789        \+       171          9               9               4        5
hek293   1:20749723-20773610     NM_016287       NM_016287       3       1       20749722        20749882        \-       160          4               4               4        0
hek293   1:20749723-20773610     NM_016287       NM_016287       4       1       20757165        20757256        \-       91           1               1               1        0
hek293   1:20749723-20773610     NM_016287       NM_016287       5       0       0               0               \0       0            0               0               0        0
hek293   1:20749723-20773610     NM_016287       NM_016287       6       0       0               0               \0       0            0               0               0        0
hek293   1:20749723-20773610     NM_016287       NM_016287       7       1       20770929        20771073        \-       144          1               1               1        0
hek293   1:20749723-20773610     NM_016287       NM_016287       8       1       20773450        20773610        \-       160          4               4               4        0
======= =====================  ================ ============  ========== =====  ============   ============= ======= =============   ==============  ===========     ========= ========

| *sample*: Sample name as specified by the user. This is useful if the user wants to merge files from different samples
| *circle_id*: circRNA-ID. The circleID is formatted to be copy and pasted to a genome browser for easy access
| *transcript_id*: Transcript name as defined by the bed-annotation file. This is the best fitting transcript. i.e. the splicing variants that contains the most exons that are actually covered
| *other_ids*: Alternative Transcript names that are either just as fitting, or contain more or less exons as supported by reads
| *exon_id*: Exon number relative to the host-gene of the circularized exon. One circle may have more than one exon. These will be listed as consecutive lines
| *chr*: Chromosome the circRNA is located on
| *start*: 5' start of the exon, relative to the reference strand, 0-based
| *end*: 3' end of the exon, relative to the reference start, 0-based
| *strand*: Strand of the host-gene
| *exon_length*: Length of the current exon
| *unique_reads*: Number of unique reads associated with the chimeric junction. When the data is paired end, then both ends are considered as separate reads.
| *fragments*: Number of broken fragments aligning to the circle
| *number\+*: Number of reads spanning the chimeric junction on the forward strand
| *number\-*: Number of reads spanning the chimeric junction on the reverse strand (if reads are only from one strand, it could indicate, that there is a sequencing bias.)
| 

**hek293.mate_status.txt:** 
This output file contains the results of analysing the amount of how often each fragment spans a chimeric junction. A fragment can either span the chimeric junction once (single), only one end spans the junction, 
twice (double) both ends span the chimeric junction, or more than twice (undefined).

=====================  ================ =============   ============   ============    ======= ======== ==========
circle_id               transcript_ids  num_reads       min_length      max_length      single  double  undefined
=====================  ================ =============   ============   ============    ======= ======== ==========
1_20749723_20773610     NM_016287       4               790              790             4       0       0
1_35358925_35361789     NM_005095       9               754              754             9       0       0
=====================  ================ =============   ============   ============    ======= ======== ==========

| *circle_id*: 
| *transcript_ids*: 
| *num_reads*: 
| *min_length*: 
| *max_length*: 
| *single*: 
| *double*: 
| *undefined*: 
| 

**hek293.skipped_exons.bed:** 

=====  ==============  ============    ==============  ======= ======= =============== ============   ========= ========== ============ =============
Chr     Circle-Start    Circle-End      Transcript      Ratio  Strand   Intron-Start    Intron-End     Color    NumExon\-2 IntronLength RelativeStart
=====  ==============  ============    ==============  ======= ======= =============== ============   ========= ========== ============ =============
chr5    178885614       178931326       NM_030613       60.0    .       178913072       178931236      255,0,0  3          1,146,1      0,30950,45711
chr6    161034259       161049979       NM_001291958    40.0    .       161049332       161049852      255,0,0  3          1,520,1      0,15073,15719
=====  ==============  ============    ==============  ======= ======= =============== ============   ========= ========== ============ =============


**hek293.skipped_exons.txt:** 

=====================   ==============  ======================  =============================================   ======================================================================================================================================   =============   ===========    
circle_id               transcript_id   skipped_exon            intron                                          read_names                                                                                                                               splice_reads    exon_reads
=====================   ==============  ======================  =============================================   ======================================================================================================================================   =============   ===========    
5_178885614_178931326   NM_030613       5:178916564-178916710   set\(\[\(\'5\', 178913072, 178931236\)\]\)      MISEQ:136:000000000-ACBC6:1:2103:10044:24618,MISEQ:136:000000000-ACBC6:1:2115:19571:6931,MISEQ:136:000000000-ACBC6:1:1119:25537:8644     3               5
6_161034259_161049979   NM_001291958    6:161049332-161049852   set\(\[\(\'6\', 161049332, 161049852\)\]\)      MISEQ:136:000000000-ACBC6:1:1113:25288:9067,MISEQ:136:000000000-ACBC6:1:2116:11815:3530                                                  2               5
=====================   ==============  ======================  =============================================   ======================================================================================================================================   =============   ===========    


--------------------

**hek293:**

1_35358925_35361789_9reads.sorted.bam
1_35358925_35361789_9reads.sorted.bam.bai
1_20749723_20773610_4reads.sorted.bam
1_20749723_20773610_4reads.sorted.bam.bai

--------------------

**hek293.coverage_pictures:**

1_35358925_35361789_NM_005095.png
1_20749723_20773610_NM_016287.png
cluster_means_all_circles.png

--------------------

**hek293.coverage_profiles:**

1_35358925_35361789.NM_005095.txt
1_20749723_20773610.NM_016287.txt
coverage.clusters.all_circles.pdf
coverage_profiles.all_circles.pdf

--------------------


========================
Error and solutions
========================
