#!/bin/bash
#
# Laura Tung
#
# Top-level script for Minimap2 + Scallop-LR assembly pipeline.
#
# Usage: minimap2_scallop_isoseq_allreads_pipeline.sh <Full_path_run_dir> <Organism> <Install_dir>
#
# <Full_path_run_dir>: The directory of a dataset (full-path all the way to BioSample ID)
# <Organism>: human or mouse
# <Install_dir>: Installation directory where this repo resides

if [ "$#" != "3" ]; then
        echo "Usage: minimap2_scallop_isoseq_allreads_pipeline.sh <Full_path_run_dir> <Organism> <Install_dir>"
        echo "<Full_path_run_dir>: The directory of a dataset (full-path all the way to BioSample ID)"
        echo "<Organism>: human or mouse"
        echo "<Install_dir>: Installation directory where this repo resides"
        echo "Run this script under the directory for a dataset, after the isoseq full analysis for this dataset is successfully completed."
        exit
fi

run_dir=$1
organism=$2

bin_dir=$3/lrassemblyanalysis/src

# Run scallop-isoseq on combined flnc and nfl reads
header_file1=$run_dir/ccs_flnc_and_nfl/ccsread_info
if [ ! -f $header_file1 ]
then
    cat $run_dir/ccs_flnc_and_nfl/flnc_and_nfl.fasta | grep ">" > $header_file1
fi
$bin_dir/minimap2_scallop_isoseq_pipeline.sh $run_dir/ccs_flnc_and_nfl flnc_and_nfl.fasta $header_file1 $organism $bin_dir


