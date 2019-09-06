#!/bin/bash
#
# Laura Tung
#
# Pipeline of Minimap2 + Scallop-LR assembly. Lower-level script not to be directly run by users.
#
# Usage: minimap2_scallop_isoseq_pipeline.sh <longreads_dir> <longreads_data_filename> <full_path_header_file> <organism> <bin_dir>
#
# <longreads_dir>: The full-path directory where the long reads data file resides
# <longreads_data_filename>: The long reads data filename (fasta or fastq)
# <full_path_header_file>: The full-path ccs info header filename.
# <organism>: human or mouse
# <bin_dir>: The directory where the scripts and python files reside.

dir=$1
longreads_data=$2
header_file=$3
organism=$4
curr_dir="$PWD"

bin_dir=$5

if [ $organism == 'mouse' ]
then
    echo "mouse"
    # REPLACE WITH YOUR ACTUAL PATH TO THE REFERENCE DATA (for the following 3 lines)
    ref_genome=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/mouse/GRCm38/GRCm38.fa
    ref_annotation=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/mouse/gtf/Mus_musculus.GRCm38.92.gtf
    ref_transcriptome=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/mouse/cDNA/Mus_musculus.GRCm38.cdna.all.fa
    local_ref_genome=$ref_genome
else
    echo "human"
    # REPLACE WITH YOUR ACTUAL PATH TO THE REFERENCE DATA (for the following 4 lines)
    ref_genome=/home/mingfus/data/transcriptomics/ensembl/human/GRCh38/GRCh38.fa
    ref_annotation=/home/mingfus/data/transcriptomics/ensembl/human/gtf/Homo_sapiens.GRCh38.90.gtf
    ref_transcriptome=/mnt/disk27/user/ltung/longreadscallop/data/longreads/mashmap/Homo_sapiens.GRCh38.cdna.all.fa
    local_ref_genome=/mnt/disk27/user/ltung/longreadscallop/data/longreads/longreads_using_enhanced_Scallop/GRCh38.fa
fi

minmap2_dir=$dir/minimap2
if [ ! -d $minmap2_dir ]
then
    mkdir $minmap2_dir
fi

# Run Minimap2
if [ ! -f $minmap2_dir/longreads.sorted.bam ]
then
    echo "Running Minimap2..."
    { time minimap2 -ax splice $ref_genome $dir/$longreads_data > $minmap2_dir/longreads.sam; } 2> $minmap2_dir/minimap2.time
    echo "Done Minimap2"

    # convert sam to bam and sort bam file
    samtools view -b -S $minmap2_dir/longreads.sam > $minmap2_dir/longreads.bam
    samtools sort $minmap2_dir/longreads.bam -o $minmap2_dir/longreads.sorted.bam
fi

scallop_dir=$minmap2_dir/scallop_isoseq
if [ ! -d $scallop_dir ]
then
    mkdir $scallop_dir
fi

# Run Scallop-LR
echo "Running Scallop-LR..."
{ time scallop -i $minmap2_dir/longreads.sorted.bam -o $scallop_dir/longreads.gtf -c $header_file --min_num_hits_in_bundle 1 | grep Bundle > $scallop_dir/test.output; } 2> $scallop_dir/scallop.time
echo "Done Scallop-LR"

# Run gffcompare
cd $scallop_dir
gffcompare -M -N -r $ref_annotation $scallop_dir/longreads.gtf

# Get ROC
if [ $organism == 'mouse' ]
then
    echo "mouse"
    gtfcuff roc $scallop_dir/gffcmp.longreads.gtf.tmap 108855 > $scallop_dir/scallop.roc
else
    echo "human"
    gtfcuff roc $scallop_dir/gffcmp.longreads.gtf.tmap 174075 > $scallop_dir/scallop.roc
fi

# Get PR-AUC
python $bin_dir/PR_auc.py $scallop_dir/scallop.roc > $scallop_dir/pr_auc

# Get potential novel isoforms
grep "class_code \"j\"" $scallop_dir/gffcmp.annotated.gtf > $scallop_dir/class_code_j_transcripts

# Fetch transcripts sequences
gffread $scallop_dir/longreads.gtf -g $local_ref_genome -w $scallop_dir/longreads.fa

# Run Transrate
if [ ! -d $scallop_dir/transrate_scallop ]
then
    mkdir $scallop_dir/transrate_scallop
fi

echo "Running Transrate..."
transrate --assembly $scallop_dir/longreads.fa --reference $ref_transcriptome --threads 8 --output $scallop_dir/transrate_scallop > $scallop_dir/transrate_scallop/transrate.log
echo "Done Transrate"
cd $curr_dir
