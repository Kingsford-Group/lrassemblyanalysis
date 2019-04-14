#!/bin/bash
#
# Laura Tung
#
# Convert bax.h5 files to subreads bam file.
#
# Copy this script to the directory for a SRA Run ID.
# Edit this file according to the current <Run_ID> and movie (bax.h5 files).
# For example, SRR6380196 is the <Run_ID>. The following three bax.h5 files belong to the current SRA Run (one movie).
# This will create subreads and scraps bam files and the corresponding index pbi files for this movie.
#
# Usage: run_bax2bam.sh

bax2bam -o SRR6380196_m151225_010700 \
 m151225_010700_42177R_c100908062550000001823201904301665_s1_p0.1.bax.h5 \
 m151225_010700_42177R_c100908062550000001823201904301665_s1_p0.2.bax.h5 \
 m151225_010700_42177R_c100908062550000001823201904301665_s1_p0.3.bax.h5
