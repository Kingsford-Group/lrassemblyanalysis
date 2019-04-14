#!/bin/bash
#
# Laura Tung
#
# Create bai index files for subreads bam files in the directory for a SRA Run ID.
# Run this script under the directory for a SRA Run ID, after run_bax2bam.sh is successfully run.
#
# Usage: create_bai.sh 


# create bai
for f1 in *.subreads.bam
do
    samtools index $f1
done


