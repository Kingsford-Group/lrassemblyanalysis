#!/bin/bash
#
# Laura Tung
#
# Perform post-analysis after the isoseq full analysis for a BioSample is completed.
#
# Usage: post_isoseq_analysis.sh <BioSample_ID> <Organism> <Install_dir>
#
# <BioSample_ID>: BioSample ID.
# <Organism>: human or mouse
# <Install_dir>: Installation directory where this repo resides
#
# Run this script under the directory for this BioSample ID, after the isoseq full analysis for this BioSample ID is successfully completed.

if [ "$#" != "3" ]; then
        echo "Usage: post_isoseq_analysis.sh <BioSample_ID> <Organism> <Install_dir>"
        echo "<BioSample_ID>: BioSample ID"
        echo "<Organism>: human or mouse"
        echo "<Install_dir>: Installation directory where this repo resides"
        echo "Run this script under the directory for this BioSample ID, after the isoseq full analysis for this BioSample ID is successfully completed."
        exit
fi

run_id=$1
organism=$2
curr_dir="$PWD"

analysis_dir=${run_id}_full_analysis

bin_dir=$3/lrassemblyanalysis/src

if [ $organism == 'mouse' ]
then
    echo "mouse"
    ref_annotation=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/mouse/gtf/Mus_musculus.GRCm38.92.gtf
    ref_transcriptome=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/mouse/cDNA/Mus_musculus.GRCm38.cdna.all.fa
else
    echo "human"
    ref_annotation=/home/mingfus/data/transcriptomics/ensembl/human/gtf/Homo_sapiens.GRCh38.90.gtf
    ref_transcriptome=/mnt/disk27/user/ltung/longreadscallop/data/longreads/mashmap/Homo_sapiens.GRCh38.cdna.all.fa
fi

# take care of ccs_fastq
mkdir ccs_fastq
ln -s $curr_dir/$analysis_dir/tasks/pbcoretools.tasks.bam2fastq_ccs-0/ccs.fastq $curr_dir/ccs_fastq/ccs.fastq

# take care of ccs_full_length and ccs_non_full_length
mkdir ccs_full_length
ln -s $curr_dir/$analysis_dir/tasks/pbcoretools.tasks.gather_contigset-2/file.contigset.fasta $curr_dir/ccs_full_length/isoseq_flnc.fasta

mkdir ccs_non_full_length
ln -s $curr_dir/$analysis_dir/tasks/pbcoretools.tasks.gather_contigset-3/file.contigset.fasta $curr_dir/ccs_non_full_length/isoseq_nfl.fasta

mkdir ccs_flnc_and_nfl
ln -s $curr_dir/$analysis_dir/tasks/pbcoretools.tasks.gather_contigset-1/file.contigset.fasta $curr_dir/ccs_flnc_and_nfl/flnc_and_nfl.fasta


# take care of final_collapsed_isoforms dir
isoforms_dir=$curr_dir/final_collapsed_isoforms

mkdir $isoforms_dir

cd $isoforms_dir

ln -s $curr_dir/$analysis_dir/tasks/pbtranscript2tools.tasks.post_mapping_to_genome-0/output_mapped.fastq .

ln -s $curr_dir/$analysis_dir/tasks/pbtranscript2tools.tasks.post_mapping_to_genome-0/output_mapped.gff .

# convert fastq to fasta
bioawk -c fastx '{print ">"$name"\n"$seq}' output_mapped.fastq > output_mapped.fasta

# Run gffcompare
/home/mingfus/data/tools/bin/gffcompare -M -N -r $ref_annotation $isoforms_dir/output_mapped.gff

# Get ROC
if [ $organism == 'mouse' ]
then
    echo "mouse"
    gtfcuff roc $isoforms_dir/gffcmp.output_mapped.gff.tmap 108855 > $isoforms_dir/output_mapped.roc
else
    echo "human"
    gtfcuff roc $isoforms_dir/gffcmp.output_mapped.gff.tmap 174075 > $isoforms_dir/output_mapped.roc
fi

# Get PR-AUC
python $bin_dir/PR_auc.py $isoforms_dir/output_mapped.roc > $isoforms_dir/pr_auc

# Get potential novel isoforms
grep "class_code \"j\"" $isoforms_dir/gffcmp.annotated.gtf > $isoforms_dir/class_code_j_transcripts

# Run Transrate
mkdir $isoforms_dir/transrate_isoforms

echo "Running Transrate..."
transrate --assembly $isoforms_dir/output_mapped.fasta --reference $ref_transcriptome --threads 8 --output $isoforms_dir/transrate_isoforms > $isoforms_dir/transrate_isoforms/transrate.log
echo "Done Transrate"
cd $curr_dir
