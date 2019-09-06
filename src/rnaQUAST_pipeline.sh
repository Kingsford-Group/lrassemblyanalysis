#!/bin/bash
#
# Laura Tung
#
# Pipeline of rnaQUAST.
#
# Usage: rnaQUAST_pipeline.sh <Full_path_run_dir> <Organism> <Merge_dir>
#
# <Full_path_run_dir>: The directory of a dataset (full-path all the way to BioSample ID)
# <Organism>: human or mouse
# <Merge_dir>: ccs_flnc_and_nfl or ccs_merge

if [ "$#" != "3" ]; then
        echo "Usage: rnaQUAST_pipeline.sh <Full_path_run_dir> <Organism> <Merge_dir>"
        echo "<Full_path_run_dir>: The directory of a dataset (full-path all the way to BioSample ID)"
        echo "<Organism>: human or mouse"
        echo "<Merge_dir>: ccs_flnc_and_nfl or ccs_merge"
        echo "Run this script under the directory for a dataset."
        exit
fi

dir=$1
organism=$2
merge_dir=$3
curr_dir="$PWD"

bin_dir=/mnt/disk27/user/ltung/longreadscallop/bin

if [ $organism == 'mouse' ]
then
    echo "mouse"
    # REPLACE WITH YOUR ACTUAL PATH TO THE REFERENCE DATA (for the following 4 lines)
    ref_genome=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/mouse/GRCm38/GRCm38.fa
    ref_annotation=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/mouse/gtf/Mus_musculus.GRCm38.92.gtf
    gmap_index=/mnt/disk27/user/ltung/longreadscallop/data/gmapdb/GRCm38
    gene_db=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/mouse/gene_db/Mus_musculus.GRCm38.92.db
else
    echo "human"
    # REPLACE WITH YOUR ACTUAL PATH TO THE REFERENCE DATA (for the following 4 lines)
    ref_genome=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/human/GRCh38/GRCh38.fa
    ref_annotation=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/human/gtf/Homo_sapiens.GRCh38.90.gtf
    gmap_index=/mnt/disk27/user/ltung/longreadscallop/data/gmapdb/GRCh38
    gene_db=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/human/gene_db/Homo_sapiens.GRCh38.90.db
fi

minimap2_dir=$dir/$merge_dir/minimap2

scallop_dir=$minimap2_dir/scallop_isoseq
stringtie_dir=$minimap2_dir/stringtie_1
isoforms_dir=$dir/final_collapsed_isoforms

cd $minimap2_dir
if [ ! -d $minimap2_dir/rnaQUAST_output ]
then
    mkdir $minimap2_dir/rnaQUAST_output
fi

# Run rnaQUAST round 1
echo "Running rnaQUAST round 1..."
python $bin_dir/rnaQUAST.py --transcripts $scallop_dir/longreads.fa $stringtie_dir/longreads.fa $isoforms_dir/output_mapped.fasta --reference $ref_genome --gene_db $gene_db --gmap_index $gmap_index --output_dir $minimap2_dir/rnaQUAST_output --threads 40 --labels Scallop StringTie IsoSeq --no_plots --disable_infer_genes --disable_infer_transcripts --lower_threshold 0.75 --upper_threshold 0.95 > $minimap2_dir/rnaQUAST_output/rnaQUAST.log
echo "Done rnaQUAST round 1"

if [ ! -d $minimap2_dir/rnaQUAST_output_1 ]
then
    mkdir $minimap2_dir/rnaQUAST_output_1
fi

# Run rnaQUAST round 2
echo "Running rnaQUAST round 2..."
python $bin_dir/rnaQUAST.py --transcripts $scallop_dir/longreads.fa $stringtie_dir/longreads.fa $isoforms_dir/output_mapped.fasta --reference $ref_genome --gene_db $gene_db --gmap_index $gmap_index --output_dir $minimap2_dir/rnaQUAST_output_1 --threads 40 --labels Scallop StringTie IsoSeq --no_plots --disable_infer_genes --disable_infer_transcripts --lower_threshold 0.0 --upper_threshold 0.5 > $minimap2_dir/rnaQUAST_output_1/rnaQUAST.log
echo "Done rnaQUAST round 2"

cd $curr_dir
