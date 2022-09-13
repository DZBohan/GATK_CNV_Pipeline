#!/bin/bash

module load gcc/11.2.0
module load jdk/17.0.1

source config_gatk_cnv_prepare.txt

# total two output files of these two pre-steps:
# targets.preprocessed.interval_list
# annotated_intervals.tsv
# they are required for the gatk cnv analysis steps.

# pre-step1 preprocess intervals

gatk PreprocessIntervals \
        -L ${bed} \
        -R ${refer} \
        --bin-length 0 \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O ${output_path}/targets.preprocessed.interval_list

# pre-step2 annotate intervals

gatk AnnotateIntervals \
        -R ${refer} \
        -L ${output_path}/targets.preprocessed.interval_list \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O ${output_path}/annotated_intervals.tsv
