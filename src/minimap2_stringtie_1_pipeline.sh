#!/bin/bash
#
# Laura Tung
#
# Pipeline of Minimap2 + StringTie assembly. Lower-level script not to be directly run by users.
#
# Usage: minimap2_stringtie_1_pipeline.sh <longreads_dir> <longreads_data_filename> <organism> <bin_dir>
#
# <longreads_dir>: The full-path directory where the long reads data file resides
# <longreads_data_filename>: The long reads data filename (fasta or fastq)
# <organism>: human or mouse
# <bin_dir>: The directory where the scripts and python files reside.

dir=$1
longreads_data=$2
organism=$3
curr_dir="$PWD"

bin_dir=$4

if [ $organism == 'mouse' ]
then
    echo "mouse"
    ref_genome=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/mouse/GRCm38/GRCm38.fa
    ref_annotation=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/mouse/gtf/Mus_musculus.GRCm38.92.gtf
    ref_transcriptome=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/mouse/cDNA/Mus_musculus.GRCm38.cdna.all.fa
    local_ref_genome=$ref_genome
else
    echo "human"
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
if [ ! -f $minmap2_dir/longreads.sorted.XS.bam ]
then
    /home/mingfus/data/repositories/bamkit/build/src/bamkit ts2XS $minmap2_dir/longreads.sorted.bam $minmap2_dir/longreads.sorted.XS.bam
fi

stringtie_dir=$minmap2_dir/stringtie_1
if [ ! -d $stringtie_dir ]
then
    mkdir $stringtie_dir
fi

# Run StringTie
echo "Running StringTie..."
{ time /home/mingfus/data/tools/bin/stringtie $minmap2_dir/longreads.sorted.XS.bam -o $stringtie_dir/longreads.gtf -c 1.0; } 2> $stringtie_dir/stringtie.time
echo "Done StringTie"

# Run gffcompare
cd $stringtie_dir
/home/mingfus/data/tools/bin/gffcompare -M -N -r $ref_annotation $stringtie_dir/longreads.gtf

# Get ROC
if [ $organism == 'mouse' ]
then
    echo "mouse"
    gtfcuff roc $stringtie_dir/gffcmp.longreads.gtf.tmap 108855 > $stringtie_dir/stringtie.roc
else
    echo "human"
    gtfcuff roc $stringtie_dir/gffcmp.longreads.gtf.tmap 174075 > $stringtie_dir/stringtie.roc
fi

# Get PR-AUC
python $bin_dir/PR_auc.py $stringtie_dir/stringtie.roc > $stringtie_dir/pr_auc

# Get potential novel isoforms
grep "class_code \"j\"" $stringtie_dir/gffcmp.annotated.gtf > $stringtie_dir/class_code_j_transcripts

# Fetch transcripts sequences
gffread $stringtie_dir/longreads.gtf -g $local_ref_genome -w $stringtie_dir/longreads.fa

# Run Transrate
if [ ! -d $stringtie_dir/transrate_stringtie ]
then
    mkdir $stringtie_dir/transrate_stringtie
fi

echo "Running Transrate..."
transrate --assembly $stringtie_dir/longreads.fa --reference $ref_transcriptome --threads 8 --output $stringtie_dir/transrate_stringtie > $stringtie_dir/transrate_stringtie/transrate.log
echo "Done Transrate"
cd $curr_dir
