#!/bin/bash
#
# Laura Tung
#
# Process rnaQUAST results.
#
# Usage: process_rnaQUAST_results.sh <full_path_run_dir> <merge_dir> <install_dir>
#
# <full_path_run_dir>: The directory of a dataset (full-path all the way to BioSample ID)
# <merge_dir>: ccs_flnc_and_nfl or ccs_merge
# <install_dir>: Installation directory where this repo resides

if [ "$#" != "3" ]; then
        echo "Usage: process_rnaQUAST_results.sh <full_path_run_dir> <merge_dir> <install_dir>"
        echo "<full_path_run_dir>: The directory of a dataset (full-path all the way to BioSample ID)"
        echo "<merge_dir>: ccs_flnc_and_nfl or ccs_merge"
        echo "<install_dir>: Installation directory where this repo resides"
        echo "Run this script under the directory for a dataset."
        exit
fi

dir=$1
merge_dir=$2
curr_dir="$PWD"

bin_dir=$3/lrassemblyanalysis/src

minimap2_dir=$dir/$merge_dir/minimap2

# take care of bins for isoform and transcript
cd $minimap2_dir/rnaQUAST_output
sed -n '34,35 p' short_report.txt > isoform_data
sed -n '42,43 p' short_report.txt > matched_data

cd $minimap2_dir/rnaQUAST_output_1
sed -n '34,35 p' short_report.txt > isoform_data
sed -n '42,43 p' short_report.txt > matched_data

cd $minimap2_dir
python $bin_dir/compute_rnaQuast.py $minimap2_dir > rnaQUAST_bins

# take care of mean fractions for isoform and transcript
sed -n '39 p' rnaQUAST_output/short_report.txt > rnaQUAST_output/mean_isoform
sed -n '45 p' rnaQUAST_output/short_report.txt > rnaQUAST_output/mean_transcript

cd $curr_dir

