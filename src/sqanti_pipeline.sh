#!/bin/bash
#
# Laura Tung
#
# Pipeline of SQANTI.
#
# Usage: sqanti_pipeline.sh <Full_path_run_dir> <Organism> <Merge_dir>
#
# <Full_path_run_dir>: The directory of a dataset (full-path all the way to BioSample ID)
# <Organism>: human or mouse
# <Merge_dir>: ccs_flnc_and_nfl or ccs_merge

if [ "$#" != "3" ]; then
        echo "Usage: sqanti_pipeline.sh <Full_path_run_dir> <Organism> <Merge_dir>"
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
    # REPLACE WITH YOUR ACTUAL PATH TO THE REFERENCE DATA (for the following 3 lines)
    ref_genome=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/mouse/GRCm38/GRCm38.fa
    ref_annotation=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/mouse/gtf/Mus_musculus.GRCm38.92.gtf
    gmap_index=/mnt/disk27/user/ltung/longreadscallop/data/gmapdb/GRCm38
else
    echo "human"
    # REPLACE WITH YOUR ACTUAL PATH TO THE REFERENCE DATA (for the following 3 lines)
    ref_genome=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/human/GRCh38/GRCh38.fa
    ref_annotation=/mnt/disk27/user/ltung/longreadscallop/data/genomes/ensembl/human/gtf/Homo_sapiens.GRCh38.90.gtf
    gmap_index=/mnt/disk27/user/ltung/longreadscallop/data/gmapdb/GRCh38
fi

minimap2_dir=$dir/$merge_dir/minimap2

scallop_dir=$minimap2_dir/scallop_isoseq
stringtie_dir=$minimap2_dir/stringtie_1
isoforms_dir=$dir/final_collapsed_isoforms

cd $isoforms_dir
if [ ! -d $isoforms_dir/sqanti ]
then
    mkdir $isoforms_dir/sqanti
fi

# Run SQANTI for Iso-Seq
echo "Running SQANTI for Iso-Seq..."
sqanti_qc.py -g -o isoseq -d $isoforms_dir/sqanti $isoforms_dir/output_mapped.gff $ref_annotation $ref_genome > $isoforms_dir/sqanti/sqanti.log
echo "Done SQANTI for Iso-Seq"


cd $scallop_dir
if [ ! -d $scallop_dir/sqanti ]
then
    mkdir $scallop_dir/sqanti
fi

# Run SQANTI for Scallop-LR
echo "Running SQANTI for Scallop-LR..."
sqanti_qc.py -g -o scallop -d $scallop_dir/sqanti $scallop_dir/longreads.gtf $ref_annotation $ref_genome > $scallop_dir/sqanti/sqanti.log
echo "Done SQANTI for Scallop-LR"


### For now StringTie gtf not working with SQANTI, so comment out the following code for now
#cd $stringtie_dir
#if [ ! -d $stringtie_dir/sqanti ]
#then
#    mkdir $stringtie_dir/sqanti
#fi

## Run SQANTI for StringTie
#echo "Running SQANTI for StringTie..."
#sqanti_qc.py -g -o stringtie -d $stringtie_dir/sqanti $stringtie_dir/longreads.gtf $ref_annotation $ref_genome > $stringtie_dir/sqanti/sqanti.log
#echo "Done SQANTI for StringTie"

cd $curr_dir
echo "Done SQANTI"

