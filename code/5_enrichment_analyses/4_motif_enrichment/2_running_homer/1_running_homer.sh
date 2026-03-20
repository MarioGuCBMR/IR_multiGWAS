#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=1200

module load perl/5.38.0
module load samtools/1.20
module load gcc/13.2.0
module load openjdk/20.0.0
module load R
module load homer/1.0

findMotifsGenome=/opt/software/homer/1.0/bin/findMotifsGenome.pl

##############################################
#Let's do this first with the 282 IR variants#
##############################################

#Adipose tissue peaks (done)

input=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/4_homer/input/282_adipose_consensus
output=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/4_homer/output/282_adipose_consensus

$findMotifsGenome $input/IR_peaks.txt hg19 $output -bg $input/non_IR_peaks.txt -size given

#SGBS day 0 (done)

input=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/4_homer/input/282_sgbs_day_0
output=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/4_homer/output/282_sgbs_day_0

$findMotifsGenome $input/IR_peaks.txt hg19 $output -bg $input/non_IR_peaks.txt -size given

#SGBS day 4 (done)

input=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/4_homer/input/282_sgbs_day_4
output=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/4_homer/output/282_sgbs_day_4

$findMotifsGenome $input/IR_peaks.txt hg19 $output -bg $input/non_IR_peaks.txt -size given

#SGBS day 14 (done)

input=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/4_homer/input/282_sgbs_day_14
output=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/4_homer/output/282_sgbs_day_14

$findMotifsGenome $input/IR_peaks.txt hg19 $output -bg $input/non_IR_peaks.txt -size given



####################################################################
#Let's run this for SGBS at day 14 and our BMI-independent variants# (done)
####################################################################

input=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/4_homer/input/141_sgbs_day_14
output=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/4_homer/output/141_sgbs_day_14

$findMotifsGenome $input/IR_peaks.txt hg19 $output -bg $input/non_IR_peaks.txt -size given

#########################################################
#Let's run this for SGBS at day 4 and our BMI-decreasing# (in process)
#########################################################

input=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/4_homer/input/63_sgbs_day_4
output=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/4_homer/output/63_sgbs_day_4

$findMotifsGenome $input/IR_peaks.txt hg19 $output -bg $input/non_IR_peaks.txt -size given
