#!/bin/bash -l
#SBATCH -A b2015054
#SBATCH --mail-type=all
#SBATCH --mail-user=alexander.vergara@slu.se
#SBATCH -N 1
#SBATCH --time=12:00:00

module load bioinfo-tools
module load htseq

htseq-count -f bam -r pos -m union -s no -t exon -i Parent /proj/b2015054/nobackup/Merged_technical_replicates/merge_P1406_102_control.bam Eugene.gff3 > htseq_count_Eugene_vs_P1406_102_control_sortmerna_trimmomatic_STAR.txt
