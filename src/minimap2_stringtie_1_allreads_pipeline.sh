#!/bin/bash
#
# Laura Tung
#
# Top-level script for Minimap2 + StringTie assembly pipeline.
#
# Usage: minimap2_stringtie_1_allreads_pipeline.sh <Full_path_run_dir> <Organism> <Install_dir>
#
# <Full_path_run_dir>: The directory of a dataset (full-path all the way to BioSample ID)
# <Organism>: human or mouse
# <Install_dir>: Installation directory where this repo resides

if [ "$#" != "3" ]; then
        echo "Usage: minimap2_stringtie_1_allreads_pipeline.sh <Full_path_run_dir> <Organism> <Install_dir>"
        echo "<Full_path_run_dir>: The directory of a dataset (full-path all the way to BioSample ID)"
        echo "<Organism>: human or mouse"
        echo "<Install_dir>: Installation directory where this repo resides"
        echo "Run this script under the directory for a dataset, after the isoseq full analysis for this dataset is successfully completed."
        exit
fi

run_dir=$1
organism=$2

bin_dir=$3/lrassemblyanalysis/src


# Run stringtie_1 on combined flnc and nfl
$bin_dir/minimap2_stringtie_1_pipeline.sh $run_dir/ccs_flnc_and_nfl flnc_and_nfl.fasta $organism $bin_dir

