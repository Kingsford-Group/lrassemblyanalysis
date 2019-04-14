#!/bin/bash
#
# Laura Tung
#
# Perform isoseq full analysis for a dataset corresponding to a BioSample. 
# First, create the merged subreads dataset from all the SRA Runs of current BioSample before performing the full analysis.
#
# Usage: biosample_isoseq.sh <BioSample_ID> <Top_dir> <Organism>
#
# <BioSample_ID>: BioSample ID. Under the BioSamples/ directory of the SRA Study directory, each BioSample has a directory named by the BioSample ID. 
#                 The directory for this BioSample ID should have a 'SRA_Runs' file that contains all the SRA Run ID's of this BioSample, and each line is one SRA Run ID.
# <Top_dir>: Top-level directory of SRA Study (all the SRA Runs of this SRA Study are under this directory).
# <Organism>: human or mouse
#
# Run this script under the directory for this BioSample ID.

if [ "$#" != "3" ]; then
	echo "Usage: biosample_isoseq.sh <BioSample_ID> <Top_dir> <Organism>"
        echo "<BioSample_ID>: BioSample ID. Under the BioSamples/ directory of the SRA Study directory, each BioSample has a directory named by the BioSample ID."
        echo "                The directory for this BioSample ID should have a 'SRA_Runs' file that contains all the SRA Run ID's of this BioSample, and each line is one SRA Run ID."
        echo "<Top_dir>: Top-level directory of SRA Study (all the SRA Runs of this SRA Study are under this directory)"
        echo "<Organism>: human or mouse"
        echo "Run this script under the directory for this BioSample ID"
	exit
fi

run_id=$1
# note: run_id here is the BioSample_ID.
top_dir=$2
organism=$3
curr_dir="$PWD"

# GMAP reference sets (pre-generated) locations for human and mouse
human_gmap_refset="/mnt/disk27/user/ltung/longreadscallop/data/datasets/PacBio/human/ERX1468898/ERR1397639/GMAP_Ref_GRCh38/GmapReferenceSet_GRCh38"
mouse_gmap_refset="/mnt/disk27/user/ltung/longreadscallop/data/datasets/PacBio/mouse/GMAP_Ref_GRCm38/GmapReferenceSet_GRCm38"

# make soft links to the bam files in all the SRA Runs of current BioSample
filename="SRA_Runs"
while read -r line
do
    sra_run_id=$line
    ln -s $top_dir/$sra_run_id/*.bam* .
done < $filename

# create dataset
dataset create --type SubreadSet ${run_id}.subreadset.xml *.subreads.bam

# perform full analysis
analysis_dir=${run_id}_full_analysis

mkdir $analysis_dir

if [ $organism == 'mouse' ]
then
    echo "mouse"
    pbsmrtpipe pipeline-id pbsmrtpipe.pipelines.sa3_ds_isoseq2_with_genome -e eid_subread:$curr_dir/${run_id}.subreadset.xml eid_gmapref_dataset:$mouse_gmap_refset/gmapreferenceset.xml -o $analysis_dir --local-only --force-chunk-mode > full_analysis.log
else
    echo "human"
    pbsmrtpipe pipeline-id pbsmrtpipe.pipelines.sa3_ds_isoseq2_with_genome -e eid_subread:$curr_dir/${run_id}.subreadset.xml eid_gmapref_dataset:$human_gmap_refset/gmapreferenceset.xml -o $analysis_dir --local-only --force-chunk-mode > full_analysis.log 
fi

