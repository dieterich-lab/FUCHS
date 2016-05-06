#!/usr/bin/bash
#SBATCH --time=2-10:00:00
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH -p blade
#SBATCH --cpus-per-task=1


# 
# samtools view -Sb -o /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/FUCHS/frontal_cortex.1.bam /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/alignments/star/frontal_cortex/mate1.Chimeric.out.sam
# samtools view -Sb -o /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/FUCHS/frontal_cortex.2.bam /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/alignments/star/frontal_cortex/mate2.Chimeric.out.sam
# 
# samtools sort /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/FUCHS/frontal_cortex.1.bam /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/FUCHS/frontal_cortex.1.sorted
# samtools sort /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/FUCHS/frontal_cortex.2.bam /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/FUCHS/frontal_cortex.2.sorted
# 
# samtools index /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/FUCHS/frontal_cortex.1.sorted.bam
# samtools index /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/FUCHS/frontal_cortex.2.sorted.bam
# 
# samtools merge /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/FUCHS/frontal_cortex.bam /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/alignments/star/frontal_cortex/mate1.Aligned.sortedByCoord.out.bam /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/alignments/star/frontal_cortex/mate2.Aligned.sortedByCoord.out.bam /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/FUCHS/frontal_cortex.1.sorted.bam /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/FUCHS/frontal_cortex.2.sorted.bam
# 
# samtools index /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/FUCHS/frontal_cortex.bam

python FUCHS.py -r 4 -q 4 -p ensembl -e 2 -c /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/FUCHS/CircRNACount -m /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/FUCHS/frontal_cortex.mate1.Chimeric.out.junction.fixed -j /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/FUCHS/frontal_cortex.mate2.Chimeric.out.junction.fixed mock /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/FUCHS/frontal_cortex.bam /beegfs/group_dv/home/FMetge/genomes/mus_musculus/GRCm38_79/mm10.ensembl.exons.bed /beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/FUCHS/ frontal_cortex

