# Overview: Long-Read Transcript Assembly Analysis

This repository analyzes PacBio transcriptomic long-read datasets extracted from [SRA](https://www.ncbi.nlm.nih.gov/sra) using
[**Scallop-LR**](https://github.com/Kingsford-Group/scallop/releases/tag/isoseq-v0.9.1), [Iso-Seq Analysis](https://www.pacb.com/documentation/smrt-link-software-installation-v5-1-0/), and [StringTie](https://ccb.jhu.edu/software/stringtie/). Scallop-LR is our released long-read transcript assembler, and StringTie is a 
leading short-read transcript assembler which can also assemble long reads. Iso-Seq Analysis is a software 
system developed by PacBio that takes subreads as input and outputs polished consensus isoforms (transcripts). 
The predicted transcripts from Iso-Seq Analysis, Scallop-LR, and StringTie are evaluated using multiple 
evaluation methods [Gffcompare](http://ccb.jhu.edu/software/stringtie/gff.shtml), [SQANTI](https://bitbucket.org/ConesaLab/sqanti), [rnaQUAST](http://cab.spbu.ru/software/rnaquast/), and [Transrate](http://hibberdlab.com/transrate/). 
   
# Datasets and the Directory Structure for Datasets

In most of the PacBio datasets in SRA, one BioSample has multiple SRA Runs because the experimenters used 
multiple "movies" to increase the coverage so that low-abundance, long isoforms can be captured in analysis. 
The experimenters also used a "size selection" sequencing strategy, and thus different SRA Runs are designated 
for different size ranges. Therefore, we use one BioSample instead of one SRA Run to represent one dataset 
in our analysis, and we merge multiple SRA Runs that belong to the same BioSample into that dataset.

The following are the 26 datasets used in the analysis with their corresponding SRA Study ID's and BioSample ID's.
Each dataset corresponds to one BioSample and named by the BioSample ID (except that datasets 15-18 are four replicates for one BioSample).
The data can be extracted from [SRA](https://www.ncbi.nlm.nih.gov/sra), preprocessed and merged into a BioSample-based dataset
using the scripts in this repository. 

Dataset	| BioSample	| SRA Study	| Organism
--------|-----------|-----------|---------
1	| SAMN00001694	| ERP015321	| Homo sapiens
2	| SAMN00001695	| ERP015321	| Homo sapiens
3	| SAMN00001696	| ERP015321	| Homo sapiens
4	| SAMN00006465	| ERP015321	| Homo sapiens
5	| SAMN00006466	| ERP015321	| Homo sapiens
6	| SAMN00006467	| ERP015321	| Homo sapiens
7	| SAMN00006579	| ERP015321	| Homo sapiens
8	| SAMN00006580	| ERP015321	| Homo sapiens
9	| SAMN00006581	| ERP015321	| Homo sapiens
10	| SAMN08182059	| SRP126849	| Homo sapiens
11	| SAMN08182060	| SRP126849	| Homo sapiens
12	| SAMN04563763	| SRP071928	| Homo sapiens
13	| SAMN07611993	| SRP098984	| Homo sapiens
14	| SAMN04169050	| SRP068953	| Homo sapiens
15	| SAMN04251426.1	| SRP065930	| Homo sapiens
16	| SAMN04251426.2	| SRP065930	| Homo sapiens
17	| SAMN04251426.3	| SRP065930	| Homo sapiens
18	| SAMN04251426.4	| SRP065930	| Homo sapiens
19	| SAMEA3374575	| ERP010189	| Mus musculus
20	| SAMEA3374576	| ERP010189	| Mus musculus
21	| SAMEA3374577	| ERP010189	| Mus musculus
22	| SAMEA3374578	| ERP010189	| Mus musculus
23	| SAMEA3374579	| ERP010189	| Mus musculus
24	| SAMEA3374580	| ERP010189	| Mus musculus
25	| SAMEA3374581	| ERP010189	| Mus musculus
26	| SAMEA3374582	| ERP010189	| Mus musculus


We have both human and mouse datasets, and they are grouped under `human/` and `mouse/` directories. Under them, 
there is a directory for each SRA Study (named by the SRA Study ID). Under each SRA Study directory, there are 
a set of directories for all the SRA Runs that belong to this SRA Study, each is named by an SRA Run ID. There 
is also a directory called `BioSamples/`. Under `BioSamples/`, there are a set of directories for all the 
BioSamples that belong to this SRA Study, each is named by a BioSample ID. Each BioSample directory is dedicated 
to a BioSample-based dataset. Each BioSample directory should have a file called `SRA_Runs` that contains all 
the SRA Run ID's that belong to this BioSample, and each line of this file is one SRA Run ID.
   
Please use the above directory structure for datasets, for running the scripts in this repository for analysis.

# Tools/Programs and References Used in Analysis

The following are the tools/programs that are used in the analysis and their corresponding versions:

  Tool/Program     |  Version
  -----------------|---------------------------------
  [Iso-Seq Analysis](https://www.pacb.com/documentation/smrt-link-software-installation-v5-1-0/) |  Iso-Seq2 from SMRT Link v5.1.0.
  [Minimap2](https://github.com/lh3/minimap2)         |  v2.2.
  [StringTie](https://ccb.jhu.edu/software/stringtie/)        |  v1.3.2d.
  [Scallop-LR](https://github.com/Kingsford-Group/scallop/releases/tag/isoseq-v0.9.1)       |  v0.9.1.
  [Gffcompare](http://ccb.jhu.edu/software/stringtie/gff.shtml)       |  v0.9.9c.
  [SQANTI](https://bitbucket.org/ConesaLab/sqanti)           |  v1.2.
  [rnaQUAST](http://cab.spbu.ru/software/rnaquast/)         |  v.1.5.1.
  [Transrate](http://hibberdlab.com/transrate/)        |  v1.0.3.
  [GMAP](http://research-pub.gene.com/gmap/)             |  version 2017-09-30.
  [gffread](http://ccb.jhu.edu/software/stringtie/gff.shtml)          |  
  [gtfcuff](https://github.com/Kingsford-Group/rnaseqtools)          |  
  [bamkit](https://github.com/Shao-Group/bamkit)           |  
   
You need to download and compile these tools/programs. Please replace the locations of these executables
in the scripts of this repository by your actual locations of your executables.

The following are the reference genomes, reference annotations, and reference transcriptomes (from [Ensembl](https://uswest.ensembl.org/index.html)) 
that are used in the analysis:

  Reference                |  Human                           |  Mouse 
  -------------------------|----------------------------------|---------------------------------
  reference genome         |  GRCh38                          |  GRCm38 
  reference annotation     |  Homo_sapiens.GRCh38.90.gtf      |  Mus_musculus.GRCm38.92.gtf
  reference transcriptome  |  Homo_sapiens.GRCh38.cdna.all.fa |  Mus_musculus.GRCm38.cdna.all.fa

Please replace the locations of the reference genomes, reference annotations, reference transcriptomes, 
gene database, gmap db, and gmap reference sets (for Iso-Seq) in the scripts of this repository by your 
actual locations of them.

# Analyze a BioSample-based Dataset with Iso-Seq Analysis, Scallop-LR, and StringTie

1. Pre-Process each SRA Run that belongs to this BioSample:

    1. Under the SRA Study directory, create a directory named by this SRA Run ID.

        In this SRA Run directory, download the hdf5 of this SRA Run from SRA and untar the hdf5.


    2. Copy `run_bax2bam.sh` from this repository to the current SRA Run directory. 

        Edit `run_bax2bam.sh` according to the current SRA Run ID and movie (bax.h5 files).

        Run `run_bax2bam.sh` to convert bax.h5 files to subreads bam files:
        ```
        run_bax2bam.sh &
        ```


    3. Create bai index files for subreads bam files:
        ```
        create_bai.sh &
        ```

2. Run Iso-Seq Analysis for this BioSample:

    (1) Under the `BioSamples/` directory of the SRA Study directory, create a directory named by this BioSample ID.
        In this BioSample directory, create the `SRA_Runs` file. The `SRA_Runs` file should contain all the SRA 
        Run ID's of this BioSample, and each line is one SRA Run ID.


    (2) Create the merged dataset from all the SRA Runs of this BioSample and perform Iso-Seq full-analysis:
        ```
        biosample_isoseq.sh <BioSample_ID> <Top_dir> <Organism> &
        ```


    (3) Post analysis, including Gffcompare and Transrate on final isoforms 
        (after the Iso-Seq full-analysis is successfully completed):
        ```
        post_isoseq_analysis.sh <BioSample_ID> <Organism> <Install_dir> &
        ```

3. Run Minimap2 + Scallop-LR for this BioSample, and run Gffcompare and Transrate on the predicted transcripts:

    ```
    minimap2_scallop_isoseq_allreads_pipeline.sh <Full_path_run_dir> <Organism> <Install_dir> &
    ```

4. Run Minimap2 + StringTie for this BioSample, and run Gffcompare and Transrate on the predicted transcripts:

    ```
    minimap2_stringtie_1_allreads_pipeline.sh <Full_path_run_dir> <Organism> <Install_dir> &
    ```

5. Run rnaQUAST on Scallop-LR transcripts, StringTie transcripts, and Iso-Seq isoforms for this BioSample:

    ```
    rnaQUAST_pipeline.sh <Full_path_run_dir> <Organism> <Merge_dir> &
    ```

6. Run SQANTI on Scallop-LR transcripts and Iso-Seq isoforms for this BioSample:

    ```
    sqanti_pipeline.sh <Full_path_run_dir> <Organism> <Merge_dir> &
    ```

The descriptions for the command-line arguments of these scripts can be found at the beginning of each script.
The results of Scallop-LR and StringTie are located under the auto-generated `ccs_flnc_and_nfl/minimap2/` directory. 
The results of Iso-Seq Analysis are inside the auto-generated `final_collapsed_isoforms/` directory.

