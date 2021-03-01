#### Benchmarking 

- #### [Plant Genomics & Gene Editing Congress: USA](https://www.global-engage.com/event/plant-genomics/)

- #### [A high-quality chromosomal genome assembly of Diospyros oleifera Cheng](https://academic.oup.com/gigascience/article/9/1/giz164/5707454)

- #### GALA: gap-free chromosome-scale assembly with long reads [Preprint](https://www.biorxiv.org/content/10.1101/2020.05.15.097428v2.full.pdf) | [Github](https://github.com/ganlab/gala)


- #### [Improvements in the Sequencing and Assembly of Plant Genomes](https://www.biorxiv.org/content/10.1101/2021.01.22.427724v1.full)


They just used available datasets. Such a nice idea!


- #### [Benchmarking Oxford Nanopore read assemblers for high-quality molluscan genomes](https://www.biorxiv.org/content/10.1101/2020.12.31.424979v1.full)



- ##### Assemblers used

> The following assemblers were used for the benchmarking, including the long-read only assemblers (Canu [12], Flye [13], Wtdbg2 [14], Miniasm [15], NextDenovo (https://github.com/Nextomics/NextDenovo), NECAT [6], Raven [16], and Shasta [7]) and hybrid assemblers (MaSuRCA [17] and QuickMerge [18]). 


- ##### Canu long time computing issue

> Canu was not tested on the M. coruscus genome due to the extremely intensive computing time required for this large genome. 

- ##### Error correction point

> Previous analyses have suggested using corrected ONT reads could improve the genome assembly [19]. To check the effect of this has on the assemblies, the ONT reads that were corrected and / or trimmed by Canu and NECAT were also tested. 

- ##### Subsampling effect or coverage effect 

> To check whether including the shorter ONT reads could affect the assembly. The ONT reads were also sub-sampled with different cut-off lengths (see Table 1 for the lengths used). 

![sub](https://www.biorxiv.org/content/biorxiv/early/2021/01/02/2020.12.31.424979/T1/graphic-1.large.jpg?width=800&height=600&carousel=1)


- ##### Computing aspect

> CPU hour was calculated in the Slurm workload manager system by recording the program start and end time point. However, since the hardware configuration in each node varied, the CPU hour presented is only an indicator of the relative trend among different assemblers.



- ##### Polishing point

The assembled contigs were polished at least three times with Flye, and heterozygous contigs were removed with the purge_dup pipeline [20]. The resultant genomes were polished twice using Pilon version 1.23 [21] with Illumina reads. 


- ##### Genome quality assessment 


The genome completeness of each assembly was thoroughly monitored at each step using BUSCO v4.0.6 with odb10 metazoan dataset [22]. The genome quality of the Scaly-foot Snail assemblies were assessed by QUAST version 5.0.2 [23] comparing against the formerly published version of the genome as a reference [8]. QUAST calculates genome assembly characteristics such as N50 and total size, but also assesses mis-assemblies with minimap2. 


- ##### Command lines used:
 
 
The detailed commands and settings used for all analyses can be found in the Supplementary Information.




```Python


The commands used in this study.
## Remove Bacteria contaminated Illumina reads by Kraken2
/home/share/kraken2-2.0.8-beta/kraken2 -db /home/share/kraken2-2.0.8-beta/minikraken2_v1_8GB --threads 40 --paired SFG1.fq SFG2.fq --classified-out SFG_Bac#.fq --unclassified-out SFG_Host#.fq

#Remove the symbiont contaminated ONT reads
minimap2 -ax map-ont -t 40 ~/Data/SFG/Symbiont_genome/Bac.fasta SFG_HAC_3Kb.fq | samtools fastq -n -f 4 - > SFG_HAC_3Kb_BacFree.fq

#wtdbg2 assembly
~/App/wtdbg2/wtdbg2 -i ../Flye_5Kb_canu_trimm_ass/SFG_HA.trimmedReads.fasta.gz -t 40 -o SFG_HAC_canuTrim
~/App/wtdbg2/wtpoa-cns -t 40 -i SFG_HAC_canuTrim.ctg.lay.gz -fo SFG_HAC_canuTrim.ctg.lay.fa

# Flye assembly and polishing 
flye --nano-raw ../SFG_HAC_3Kb_BacFree.fq -g 0.3566g -o flye_ONT_HAC_3kb -t 40
flye --polish-target ../SFG_10Kb_HAC.fa --nanopore-raw ../SFG_HAC_3Kb_BacFree.fq --iterations 3 --out-dir flye_po_R3 --threads 40

# minimap + miniasm assembly
minimap2 -X -t 20 -x ava-ont ../SFGB1_clean_3Kb.fq ../SFGB1_clean_3Kb.fq > reads.paf
miniasm -f ../SFGB1_clean_3Kb.fq reads.paf > SFG_default.gfa
awk '/^S/{print ">"$2"\n"$3}' SFG_default.gfa | fold > SFG_default.fa

# canu assembly 
~/App/canu-2.0/Linux-amd64/bin/canu -fast -p SFG_HA -d SFG_HA  genomeSize=0.3566g -nanopore-raw SFG_HAC_3Kb_BacFree.fq corOutCoverage=200 corMhapSensitivity=normal correctedErrorRate=0.105 minReadLength=5000  useGrid=true gridOptions=--partition=oces

#NECAT assembly
necat.pl config SFG_HAC_config.txt
necat.pl correct SFG_HAC_config.txt
necat.pl assemble SFG_HAC_config.txt
necat.pl bridge SFG_HAC_config.txt
# the configure file
PROJECT=SFG_HAC
ONT_READ_LIST=readlist.txt
GENOME_SIZE=356600000
THREADS=40
MIN_READ_LENGTH=5000
PREP_OUTPUT_COVERAGE=70
OVLP_FAST_OPTIONS=-n 500 -z 20 -b 2000 -e 0.5 -j 0 -u 1 -a 1000
OVLP_SENSITIVE_OPTIONS=-n 500 -z 10 -e 0.5 -j 0 -u 1 -a 1000
CNS_FAST_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
CNS_SENSITIVE_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
TRIM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 1 -a 400
ASM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400
NUM_ITER=2
CNS_OUTPUT_COVERAGE=50
CLEANUP=1
USE_GRID=false
GRID_NODE=0
GRID_OPTIONS=
SMALL_MEMORY=0
FSA_OL_FILTER_OPTIONS=
FSA_ASSEMBLE_OPTIONS=
FSA_CTG_BRIDGE_OPTIONS=
POLISH_CONTIGS=true

#NextDenovo configure file
[General]
job_type = local
job_prefix = nextDenovo
task = all # 'all', 'correct', 'assemble'
rewrite = yes # yes/no
deltmp = yes
rerun = 3
parallel_jobs = 4
input_type = raw
input_fofn = ./input.fofn
workdir = ./01_rundir

[correct_option]
read_cutoff = 1k
seed_cutoff = 11653
blocksize = 1g
pa_correction = 2
seed_cutfiles = 2
sort_options = -m 1g -t 2 -k 50
minimap2_options_raw = -x ava-ont -t 10
correction_options = -p 15

[assemble_option]
minimap2_options_cns = -x ava-ont -t 10 -k17 -w17
nextgraph_options = -a 1

#Raven assembly
~/App/raven/build/bin/raven ../SFG_HAC_10Kb_BacFree.fa.gz -t 40 -p 1 >Raven_10Kb.fasta 2> log.txt

#Shasta assembly
sudo ~/App/shata_v0.4.0/shasta-Linux-0.4.0 --input ../Mcor_ONT_10kb.fa --memoryMode filesystem --memoryBacking 2M

#QuickMerge
delta-filter -r -q -l 10000 SFG_hybrid.delta > SFG_hybrid.rq.delta
quickmerge -d SFG_hybrid.rq.delta -q MaSu_HAC35x_FlyeP3_PD_Pilon2.fasta -r SFG_Flye_HAC10kb_P2_PD_pilon2.fasta -hco 5.0 -c 1.5 -l 2164900 -ml 10000 -p SFG_QM

# Purge_dup version 1.2.3
minimap2 -x map-ont SFG_FlyepolishedR3.fasta SFG_ONT_3Kb_BacFree.fq -t 40 > reads.paf
pbcstat *.paf
calcuts PB.stat -l13 -m51 -u144 > cutoffs 2>calcults.log
split_fa SFG_FlyepolishedR3.fasta > polished_3.fasta.split
minimap2 -x asm5 -DP polished_3.fasta.split polished_3.fasta.split -t 40 >split.paf
purge_dups -2 -T cutoffs -c PB.base.cov split.paf >dups.bed 2>purge_dups.log
/get_seqs dups.bed SFG_FlyepolishedR3.fasta
#or with Illumina reads
bowtie2-build -f SFG_canu_flyeP3.fasta SFG_canu --threads 40
bowtie2 -p 40 --maxins 800 -x SFG_canu -1 SFG_trim_1.fq -2 SFG_trim_2.fq 1>PE.sam 2>SE.err
samtools view -bS PE.sam >PE.bam -@ 20
~/App/purge_dups/bin/ngscstat PE.bam
~/App/purge_dups/bin/calcuts TX.stat -l13 -m51 -u144 > cutoffs 2>calcults.log
~/App/purge_dups/bin/split_fa SFG_canu_flyeP3.fasta > polished_3.fasta.split
minimap2 -x asm5 -DP polished_3.fasta.split polished_3.fasta.split -t 40 > polished_3.fasta.split.paf
~/App/purge_dups/bin/purge_dups -2 -T cutoffs -c TX.base.cov polished_3.fasta.split.paf > dups.bed 2> purge_dups.log

#Pilon error correction
bowtie2-build -f canu_FlyeP3_PD2_pilon1.fa SFS --threads 40
bowtie2 -p 40 -D 20 -R 2 -N 1 -L 18 -i S,1,0.50 --maxins 1200 -x SFS -1 SFG_clean_1.fq -2 
SFG_clean_2.fq 1>SFSPE500.sam 2> SFSPE500.err
grep -E "@|NM:" SFSPE500.sam | grep -v "XS:" > SFSPE500_uniq.sam
samtools view -bS SFSPE500_uniq.sam > SFSPE500_uniq.bam -@ 40
samtools sort SFSPE500_uniq.bam -m 5G -@ 10 -o SFSPE500_uniq_sorted.bam
java -jar ~/App/picard/picard.jar MarkDuplicates I= SFSPE500_uniq_sorted.bam O= SFSPE500_uniq_sorted_dedupe.bam METRICS_FILE=metrics.txt
samtools index SFSPE500_uniq_sorted_dedupe.bam
java -Xmx180G -jar ~/App/pilon-1.23/pilon-1.23.jar --genome canu_FlyeP3_PD2_pilon1.fa --frags SFSPE500_uniq_sorted_dedupe.bam --diploid --threads 40

#QUAST analysis
quast.py ./canu/Canu.fasta ./flye/Flye.fasta ./masurca/MaSuRCA.fasta ./miniasm/Miniasm.fasta ./necat/NECAT.fasta ./nextdenovo/NextDenovo.fasta
 ./raven/Raven.fasta ./shasta/Shasta.fasta ./wtdbg2/Wtdbg2.fasta ./quickmerger/QuickMerger.fasta -r ./Csqv1.1/Csq_v1.1.fa -t 40

#repeatmodeler and repeatmasker
BuildDatabase -name Mcor_v2.0 ../Mcor_v2.0.fasta
~/App/RepeatModeler-2.0.1/RepeatModeler -database Mcor_v2.0 -pa 40 
RepeatMasker -species all -pa 8 -div 30 Mcor_v2.0.fasta
RepeatMasker -lib /home/sunj/Data/Mcoru/Annotations/RepeatModeler/Mcor_v2.0-families.fa -pa 10 -div 30 Mcor_v2.0.fasta.masked

# Braker for Augustus training
~/App/BRAKER-2.1.5/scripts/braker.pl --genome=SFG_Flye_HAC10kb_P2_PD_pilon2.fasta.masked.masked --species=Chrysomallon --bam=all.sorted.bam --cores 40

# maker configure file
#-----Genome (these are always required)
genome=SFG_Flye_HAC10kb_P2_PD_pilon2.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----EST Evidence (for best results provide a file for at least one)
est=Trinity_all_0.97.fasta #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=Mollu_Prot_50AA_0.95.fa  #protein sequence file in fasta format (i.e. from mutiple organisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=all #select a model organism for RepBase masking in RepeatMasker
rmlib=./SFG_flye-families.fa #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=~/App/maker-3.01.03/data/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species=Chrysomallon #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
run_evm=1 #run EvidenceModeler, 1 = yes, 0 = no
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
allow_overlap= #allowed gene overlap fraction (value from 0 to 1, blank for default)

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=20 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=10000 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=50 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=1 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=20000 #length for the splitting of hits (expected max intron size for evidence alignments)
min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=1 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=300 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP=/home/sunj/temp #specify a directory other than the system default temporary directory for temporary files



```

- ##### Table of the overview of the benchmarking project


| Assembler                   | Input Reads                          | settings                                                | N50/NG50\* (assuming genome size =1593Mb) | No. of Scaffolds | Max     | No. of Scaf. Over 1Mb | Total\_length | BUSCO (metazoan, ODB10)                         | cpu hour |
| --------------------------- | ------------------------------------ | ------------------------------------------------------- | ----------------------------------------- | ---------------- | ------- | --------------------- | ------------- | ----------------------------------------------- | -------- |
| shasta v0.4.0               | Raw ONT, 15Kb cutoff (65X)           | default                                                 | 218.9Kb/317.6Kb                           | 22499            | 4.03Mb  | 77                    | 2.14Gb        | C:69.1%\[S:65.1%,D:4.0%\],F:10.7%,M:20.2%,n:954 | 98       |
|                             |                                      |                                                         |                                           |                  |         |                       |               |                                                 |          |
| wtdbg2 v2.5                 | Raw ONT, 15Kb cutoff (65X)           | \-x preset3 -L 15000                                    | 721.7Kb/1.01Mb                            | 16590            | 5.96Mb  | 448                   | 2.02Gb        | C:56.1%\[S:55.6%,D:0.5%\],F:11.4%,M:32.5%,n:954 | 820      |
| wtdbg2 v2.5                 | NECAT corr, 15Kb cutoff (50X)        | \-x preset4                                             | 772.9Kb/1.11Mb                            | 16451            | 9.37Mb  | 444                   | 2.01Gb        | C:82.8%\[S:78.9%,D:3.9%\],F:4.7%,M:12.5%,n:954  | 9218     |
| wtdbg2 v2.5                 | NECAT corr & trim, 15Kb cutoff (50X) | \-x preset4                                             | 756.3Kb/1.12Mb                            | 15960            | 9.34Mb  | 455                   | 1.99Gb        | C:83.4%\[S:79.9%,D:3.5%\],F:5.1%,M:11.5%,n:954  | 21370    |
|                             |                                      |                                                         |                                           |                  |         |                       |               |                                                 |          |
| flye v2.7.1-b1590           | NECAT corr, 15Kb cutoff (50X)        | \--nano-corr                                            | 250.7Kb/403.8Kb                           | 28333            | 3.07Mb  | 116                   | 2.52Gb        | \-                                              | 11428    |
| flye v2.7.1-b1590           | Raw ONT, 23Kb cutoff (46X)           | \--nano-raw                                             | 274.8Kb/456.2Kb                           | 27310            | 2.77Mb  | 160                   | 2.57Gb        | C:87.3%\[S:75.1%,D:12.2%\],F:5.8%,M:6.9%,n:954  | 1440     |
|                             |                                      |                                                         |                                           |                  |         |                       |               |                                                 |          |
| Raven v1.1.10               | Raw ONT, 15Kb cutoff (65X)           | default                                                 | 325.2Kb/475.9Kb                           | 10382            | 1.73Mb  | 95                    | 2.55Gb        | C:87.5%\[S:78.0%,D:9.5%\],F:3.8%,M:8.7%,n:954   | 1240     |
|                             |                                      |                                                         |                                           |                  |         |                       |               |                                                 |          |
| NECAT 20200119\_Linux-amd64 | Raw ONT, 15Kb cutoff (65X)           | PREP\_OUTPUT\_COVERAGE=70, CNS\_OUTPUT\_COVERAGE=50     | 1.24Mb/1.87Mb                             | 5153             | 6.56Mb  | 887                   | 2.62Gb        | C:86.2%\[S:68.9%,D:17.3%\],F:4.7%,M:9.1%,n:954  | 23738    |
|                             |                                      | Flye Polish three times                                 | 1.25Mb/1.88Mb                             | 5076             | 6.57Mb  | 895                   | 2.63Gb        | C:90.3%\[S:66.0%,D:24.3%\],F:3.1%,M:6.6%,n:954  |          |
|                             |                                      | Flye Polish three times, purge duplicate                | 1.46Mb/1.76Mb                             | 2841             | 6.57Mb  | 711                   | 1.96Gb        | C:89.2%\[S:87.6%,D:1.6%\],F:3.8%,M:7.0%,n:954   |          |
|                             |                                      | Flye Polish three times, purge duplicate 2nd round      | 1.46Mb/1.74Mb                             | 2646             | 6.57Mb  | 703                   | 1.92Gb        | C:89.4%\[S:87.9%,D:1.5%\],F:3.8%,M:6.8%,n:954   |          |
|                             |                                      |                                                         |                                           |                  |         |                       |               |                                                 |          |
|                             |                                      |                                                         |                                           |                  |         |                       |               |                                                 |          |
| NextDenovo v2.3.0           | Raw ONT                              | seed\_cutoff = 23987, 45X                               | 2.54Mb/3.40Mb                             | 1839             | 12.50Mb | 575                   | 2.07Gb        | C:88.6%\[S:86.6%,D:2.0%\],F:4.1%,M:7.3%,n:954   | 6700     |
|                             |                                      | Flye Polish three times                                 | 2.51Mb/3.37Mb                             | 1839             | 12.40Mb | 573                   | 2.06Gb        | C:89.8%\[S:86.7%,D:3.1%\],F:3.9%,M:6.3%,n:954   |          |
|                             |                                      | Flye Polish three times, purge duplicate                | 2.70Mb/3.36Mb                             | 1591             | 12.40Mb | 548                   | 1.93Gb        | C:89.9%\[S:88.4%,D:1.5%\],F:3.6%,M:6.5%,n:954   |          |
|                             |                                      | Flye Polish three times, purge duplicate, Pilon Round 2 | 2.70Mb/3.37Mb                             | 1591             | 12.4Mb  | 547                   | 1.93Gb        | C:95.8%\[S:94.1%,D:1.7%\],F:1.5%,M:2.7%,n:954   |          |
|                             |                                      |                                                         |                                           |                  |         |                       |               |                                                 |          |
| MaSuRCA v3.4.1              | Illumina PE150, Clean ONT            | LHE\_COVERAGE=25, FLYE\_ASSEMBLY=0                      | 813.6Kb/1.70Mb                            | 10541            | 11.02Mb | 561                   | 2.67Gb        | C:96.0%\[S:70.9%,D:25.1%\],F:1.5%,M:2.5%,n:954  | 46100    |
|                             |                                      | purge duplicate                                         | 1.17Mb/1.61Mb                             | 4570             | 11.02Mb | 515                   | 1.98Gb        | C:95.6%\[S:93.4%,D:2.2%\],F:1.8%,M:2.6%,n:954   |          |
|                             |                                      | purge duplicate, POLCA polishing                        | 1.17Mb/1.61Mb                             | 4570             | 11.02Mb | 515                   | 1.98Gb        | C:95.8%\[S:93.6%,D:2.2%\],F:1.7%,M:2.5%,n:954   |          |






- #### [Benchmarking of long-read correction methods](https://academic.oup.com/nargab/article/2/2/lqaa037/5843804)
- #### [A comprehensive evaluation of long read error correction methods](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07227-0)

- Nanopore only benchmarking Nextdenovo mentionned [Benchmarking Oxford Nanopore read assemblers for high-quality molluscan genomes](https://www.biorxiv.org/content/10.1101/2020.12.31.424979v1.full)



- Good for writing material and methods | [Ultracontinuous single haplotype genome assemblies for the domestic cat (Felis
catus) and Asian leopard cat (Prionailurus bengalensis)](https://pubmed.ncbi.nlm.nih.gov/33305796/)


- Very nice paper of benchmarking with nice gggplot2 graph on Adineta vaga paper [Overcoming uncollapsed haplotypes in long-read
2 assemblies of non-model organisms](https://www.biorxiv.org/content/10.1101/2020.03.16.993428v2.full.pdf)




Fig1: Statistics of raw assemblies obtained from the full PacBio dataset (raw assemblies), with a preliminary read-filtering step (keeping only reads larger than 15 kb), or a subsequent removal of uncollapsed haplotypes with HaploMerger2, purge_dups, or purge_haplotigs. a) Assembly scores for size, N50, completeness and haploidy. b) Long-read coverage distribution over the contigs.

![fig1](https://www.biorxiv.org/content/biorxiv/early/2020/12/15/2020.03.16.993428/F1.large.jpg?width=800&height=600&carousel=1)




Fig3:  Statistics of the PacBio and Nanopore assemblies depending on sequencing depth, with a) assembly size, b) N50, c) complete single-copy BUSCOs and d) haploidy. The assemblies were ran on five random subsamplings of the long-read datasets.


![fig2](https://www.biorxiv.org/content/biorxiv/early/2020/12/15/2020.03.16.993428/F3.large.jpg?width=800&height=600&carousel=1)





Fig4: Computational resources (RAM and CPU time) used by the assemblers. a) Maximum RAM usage and b) mean CPU time depending on sequencing depth. Canu and NextDenovo were not included in this comparison as they were run on different machines.


![fig3](https://www.biorxiv.org/content/biorxiv/early/2020/12/15/2020.03.16.993428/F4.large.jpg?width=800&height=600&carousel=1)







- Awesome review of long reads sequencing assembly [Piercing the dark matter: bioinformatics of long-range sequencing and mapping](https://www.nature.com/articles/s41576-018-0003-4)












- #### [Haplotype-resolved de novo assembly with phased assembly graphs](https://arxiv.org/abs/2008.01237)  | [nature methods redcube file](https://www.nature.com/articles/s41592-020-01056-5.epdf?sharing_token=gOs_Vf3Mp87PwxwOWBR9TdRgN0jAjWel9jnR3ZoTv0PEptI9_4gtBx6ljxr0whf0cYzRf6jNFncOK-h9I2pNj7zgjAjzTvz5DZ6OR5woNG7_ZnL517PlPbK8h-g9oucxb3hTlS62DXkNKpafRKz0oqAvl8bMAcGgXyVIDKdNBME%3D&s=03)



[Heng Li lab](https://hlilab.github.io/)


Supplementary Material for
“Haplotype-resolved de novo assembly with phased assembly graphs”
S1 Software commands
S1.1 Hifiasm
To produce primary assemblies of homozygous samples (M. musculus, Z. mays and CHM13), hifiasm (version 0.7) was run with the following command which does not purge haplotig duplications:

```bash

hifiasm -o <outputPrefix> -t <nThreads> -l0 <HiFi-reads.fasta>

```



For heterozygous samples, hifiasm was run with the following command:

```bash

hifiasm -o <outputPrefix> -t <nThreads> <HiFi-reads.fasta>

```

We added ‘-D10’ for the octoploid F. × ananassa because the default k-mer cutoff seems too low:


```bash
hifiasm -o <outputPrefix> -t <nThreads> -D10 <HiFi-reads.fasta>
```


For trio-binning assembly, we first built the paternal trio index and the maternal trio index by yak (version
r55) with the following commands:

```bash

yak count -b37 -t <nThreads> -o <pat.yak> <paternal-short-reads.fastq>
yak count -b37 -t <nThreads> -o <mat.yak> <maternal-short-reads.fastq>

```

and then we produced the paternal assembly and the maternal assembly with the following command:


```bash

hifiasm -o <outputPrefix> -t <nThreads> -1 <pat.yak> -2 <mat.yak> <HiFi-reads.fasta>

```


S1.2 Falcon-Unzip
Falcon-kit (version 1.8.1) was run with the following HiFi-specific options:

```bash
length cutoff pr = 8000
ovlp daligner option = -k24 -h1024 -e.98 -l1500 -s100
ovlp HPCdaligner option = -v -B128 -M24
ovlp DBsplit option = -s400
overlap filtering setting = --max-diff 200 --max-cov 200 --min-cov 2 --n-core 24 --min-idt
98 --ignore-indels

```
Falcon-unzip-kit (version 1.3.7) was run with default options.


S1.3 HiCanu
For primary assembly, HiCanu (version 2.0) was run with the following command line:

```bash
canu -p asm -d <outDir> genomeSize=<GSize> useGrid=false maxThreads=<nThreads> \
-pacbio-hifi <HiFi-reads.fasta>



```






The contigs labeled by ‘suggestedBubbles=yes’ were removed from the primary assembly. For triobinning assembly, we ran HiCanu in two steps as recommended. We partitioned the HiFi reads by parental
short reads with the following command:




```bash
canu -haplotype -p asm -d <outDir> genomeSize=<GSize> useGrid=false \
maxThreads=<nThreads> -haplotypePat <pat-reads.fq> -haplotypeMat <mat-reads.fq> \
-pacbio-raw <HiFi-reads.fasta>

```


Note that ‘-pacbio-raw’ was used to partition HiFi reads followed the document of HiCanu. We then
perform HiCanu assemblies on partitioned reads.



S1.4 Peregrine
For primary assembly, Peregrine (version 0.1.6.1) was run with the following command, where 48 is the
number of threads in use:

```bash

docker run -it -v <workDir>:/wd --user $(id -u):$(id -g) cschin/peregrine:0.1.6.1 asm \
/wd/Input.fnfo 48 48 48 48 48 48 48 48 48 --with-consensus --with-alt --shimmer-r 3 \
--best n ovlp 8 --output <outDir>

```


For trio-binning assembly, we first used HiCanu to partition HiFi reads by parental short reads, and then
assembled the each haplotype individually by Peregrine.



S1.5 Purge dups
Purge dups (version 1.2.3) was used to postprocess the output primary assemblies of HiCanu for all heterozygous samples. The commands are as follows:

```bash



minimap2 -I6G -xmap-pb <contigs.fa> <HiFi-reads.fasta> -t <nThreads> > <read-aln.paf>
bin/pbcstat <read-aln.paf>
bin/calcuts PB.stat > cutoffs
bin/split fa <contigs.fa> > <split.fa>
minimap2 -I6G -xasm5 -DP <split.fa> <split.fa> -t <nThreads> > <ctg-aln.paf>
bin/purge dups -2 -T cutoffs -c PB.base.cov <ctg-aln.paf> > <dups.bed>
bin/get seqs <dups.bed> <contigs.fa>

```

Since running Purge Dups in default cannot produce primary assembly of HiCanu with right size for HG002,
we manually adjusted the cutoffs thresholds of Purge Dups as follows “5 7 11 30 22 42”.


S1.6 Running asmgene
We aligned the cDNAs to the reference genome and contigs by minimap2 r974 and evaluated the gene
completeness with paftools.js from the minimap2 package:

```bash


minimap2 -cxsplice:hq -t <nThreads> <contigs.fa> <cDNAs.fa> > <aln.paf>
paftools.js asmgene -i.97 <ref.paf> <asm.paf>

```


We set the sequence identity threshold to be 97% with ‘-i.97’ to tolerate low per-base accuracy of ONT
assemblies. For trio binning assemblies, we added option ‘-a’ to evaluate genes mapped to the autosomes
only. When evaluating multi-copy genes retained in an assembly, we replaced ‘-i.97’ to ‘-i.99’ to increase the resolution.



S1.7 Computing NGA50
We used minigraph (version 0.10-dirty-r361) and paftools (version 2.17-r974-dirty) to calculate the NGA50
of each asssembly:

```bash
minigraph -xasm -K1.9g --show-unmap=yes -t <nThreads> <ref.fa> <asm.fa> > <asm.paf>
paftools.js asmstat <ref.fa.fai> <asm.paf>

```

In comparison to minimap2, minigraph tends to generate longer alignments and is more robust to highly
variable regions.



S1.8 BUSCO
BUSCO (version 3.0.2) was used with the following command:

```bash

python3 run BUSCO.py -i <asm.fa> -m genome -o <outDir> -c <nThreads> -l <lineage dataset>


```

where ‘lineage dataset’ is set to tetrapoda for R. muscosa and set to embryophyta for F. × ananassa
and S. sempervirens.

S1.9 Determining resolved BACs
The resolution of BAC for different assemblies was evaluated using the pipeline at: https://github.
com/skoren/bacValidation, except that we added option ‘-I6g’ to minimap2. The BAC libraries for
CHM13 and HG00733 can be found at https://www.ncbi.nlm.nih.gov/nuccore/?term=VMRC59+
and+complete and https://www.ncbi.nlm.nih.gov/nuccore/?term=VMRC62+and+complete, respectively.





S1.10 Running yak evaluation

We used yak (version r55) to measure the per-base consensus accuracy (QV), the switch error rate and the
hamming error rate. For QV evaluation, we first built the index for the short reads coming from the same
sample:

```bash
yak count -b37 -t <nThreads> -o <sr.yak> <short-reads.fastq>
yak qv -t <nThreads> <sr.yak> <contigs.fa>


```



To evaluate the switch error rate and the hamming error rate, we first built the indexes from the paternal
short reads as in section S1.1 and then estimate k-mer based error rates as follows:


```bash
yak trioeval -t <nThreads> <pak.yak> <mat.yak> <contigs.fa>

```   
    
    
S1.11 Dipcall
For the male sample HG002, we ran dipcall (version 0.1) as follows:

```bash
dipcall.kit/run-dipcall -x dipcall.kit/hs37d5.PAR.bed <prefix> hs37d5.fa \
<pat-asm.fa> <mat-asm.fa> > <prefix.mak>
make -j2 -f <prefix.mak>



```



For the female sample HG00733, we removed option ‘-x’. We used the GRCh37 variant of ‘hs37d5.fa’
here because GIAB works best with hs37d5.


S1.12 Evaluating collapsed misassemblies for inbred samples

We used scripts at: https://github.com/lh3/CHM-eval/blob/master/misc/clustreg.js, and https:
//github.com/lh3/CHM-eval/blob/master/misc/select-collapse-het.js. The commands are
as follows:

```bash
minimap2 -axasm20 -t <contigs.fa> <HiFi-reads.fasta> | samtools sort -o <aln.bam> -
htsbox pileup -vcf <contigs.fa> -q20 -Q20 -l5000 -S5000 -s5 <aln.bam> > <var.vcf>

./select-collapse-het.js -c <readCoverage> <var.vcf> | ./clustreg.js -n10
where ‘-l’ and ‘-S’ filter out alignments short than 5kb.


```








#### comparison hi-fi and pacbio using rice [example](https://academic.oup.com/gigascience/article/9/12/giaa123/6034784)


`Contiguity of the ONT and PB assemblies. (a) Treemaps for contig length difference between the ONT` How to do [tree maps](https://www.r-graph-gallery.com/234-a-very-basic-treemap.html)





#### prokaryote [example](https://f1000research.com/articles/8-2138)

#### plant genome Macadamia jansenii [example](https://www.biorxiv.org/content/10.1101/2020.03.16.992933v1)  | [video](https://youtu.be/FuyMgjROOIk?t=150)

#### metazoa Adineta vaga [example](https://www.biorxiv.org/content/10.1101/2020.03.16.993428v2.full.pdf)

####  rice [example](https://www.biorxiv.org/content/10.1101/2020.02.13.948489v1.full.pdf)

#### Example [mixte](https://academic.oup.com/bib/article/20/3/866/4590140) | Publish my work [here?](https://academic.oup.com/bib/pages/About)

#### Example [mixte](https://pubs.acs.org/doi/10.1021/acs.jafc.0c01647)

#### Angel Dominguez [notes](https://www.google.co.kr/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwj0i7-bo67tAhUKq5QKHddcAfM4ChAWMAd6BAgKEAI&url=https%3A%2F%2Fzenodo.org%2Frecord%2F345098%2Ffiles%2Fscientific_reports_assembly_long_reads%25282%2529.pdf%3Fdownload%3D1&usg=AOvVaw1lEIcsWnug_oOWnophIZWA)

#### Usefull for introduction | [Opportunities and challenges in long-read sequencing data analysis](https://link.springer.com/article/10.1186/s13059-020-1935-5)

#### Ten step for assembly genome [project](https://f1000research.com/articles/7-148/v1)

#### Twelve steps for assembly genome [project](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008325)


#### [High-Quality Genome Assembly of Fusarium oxysporum f. sp. lini](https://www.frontiersin.org/articles/10.3389/fgene.2020.00959/full) | [The Genome Sequence of Five Highly Pathogenic Isolates of Fusarium oxysporum f. sp. lini](https://apsjournals.apsnet.org/doi/pdf/10.1094/MPMI-05-20-0130-SC)




# I have a hard time to install NECAT due to the updating of the perl 



First i followwed the instruction from the developper



```

$ wget https://github.com/xiaochuanle/NECAT/releases/download/v0.0.1_update20200803/necat_20200803_Linux-amd64.tar.gz
$ tar xzvf necat_20200803_Linux-amd64.tar.gz
$ cd NECAT/Linux-amd64/bin
$ export PATH=$PATH:$(pwd)

```


But when 



I got this error




```

syntax error at /home/kplee/programs/NECAT/Linux-amd64/bin/Plgd/Project.pm line 51, near "$cfg{"
Global symbol "%env" requires explicit package name at /home/kplee/programs/NECAT/Linux-amd64/bin/Plgd/Project.pm line 52.
Global symbol "%env" requires explicit package name at /home/kplee/programs/NECAT/Linux-amd64/bin/Plgd/Project.pm line 55.
Global symbol "%env" requires explicit package name at /home/kplee/programs/NECAT/Linux-amd64/bin/Plgd/Project.pm line 56.
syntax error at /home/kplee/programs/NECAT/Linux-amd64/bin/Plgd/Project.pm line 57, near "}"
syntax error at /home/kplee/programs/NECAT/Linux-amd64/bin/Plgd/Project.pm line 84, near "$env{"
Global symbol "$cfg" requires explicit package name at /home/kplee/programs/NECAT/Linux-amd64/bin/Plgd/Project.pm line 84.
Global symbol "@jobs" requires explicit package name at /home/kplee/programs/NECAT/Linux-amd64/bin/Plgd/Project.pm line 89.
Global symbol "$prjDir" requires explicit package name at /home/kplee/programs/NECAT/Linux-amd64/bin/Plgd/Project.pm line 95.
Global symbol "$env" requires explicit package name at /home/kplee/programs/NECAT/Linux-amd64/bin/Plgd/Project.pm line 101.
Global symbol "$env" requires explicit package name at /home/kplee/programs/NECAT/Linux-amd64/bin/Plgd/Project.pm line 115.
Global symbol "$cfg" requires explicit package name at /home/kplee/programs/NECAT/Linux-amd64/bin/Plgd/Project.pm line 115.
Global symbol "$cfg" requires explicit package name at /home/kplee/programs/NECAT/Linux-amd64/bin/Plgd/Project.pm line 121.
syntax error at /home/kplee/programs/NECAT/Linux-amd64/bin/Plgd/Project.pm line 121, near "$cfg{"
/home/kplee/programs/NECAT/Linux-amd64/bin/Plgd/Project.pm has too many errors.
Compilation failed in require at ./necat.pl line 13.
BEGIN failed--compilation aborted at ./necat.pl line 13.



```




first I install perl using


```bash

sudo yum install perl

```

perl -v show that we have 5.16


```
This is perl 5, version 16, subversion 3 (v5.16.3) built for x86_64-linux-thread-multi
(with 41 registered patches, see perl -V for more detail)

Copyright 1987-2012, Larry Wall

Perl may be copied only under the terms of either the Artistic License or the
GNU General Public License, which may be found in the Perl 5 source kit.

Complete documentation for Perl, including FAQ lists, should be found on
this system using "man perl" or "perldoc perl".  If you have access to the
Internet, point your browser at http://www.perl.org/, the Perl Home Page.





```


to update to 5.26 I found this [page](https://www.softwarecollections.org/en/scls/rhscl/rh-perl526/)


```bash

# 1. Install a package with repository for your system:
# On CentOS, install package centos-release-scl available in CentOS repository:
$ sudo yum install centos-release-scl

# On RHEL, enable RHSCL repository for you system:
# $ sudo yum-config-manager --enable rhel-server-rhscl-7-rpms

# 2. Install the collection:
$ sudo yum install rh-perl526

# 3. Start using the software collection:
$ scl enable rh-perl526 bash


```

and then I checked

perl -v and bingo I got


```

This is perl 5, version 26, subversion 3 (v5.26.3) built for x86_64-linux-thread-multi
(with 27 registered patches, see perl -V for more detail)

Copyright 1987-2018, Larry Wall

Perl may be copied only under the terms of either the Artistic License or the
GNU General Public License, which may be found in the Perl 5 source kit.

Complete documentation for Perl, including FAQ lists, should be found on
this system using "man perl" or "perldoc perl".  If you have access to the
Internet, point your browser at http://www.perl.org/, the Perl Home Page.

```


so I ran 


./necat.pl


and I got this

```

Smartmatch is experimental at /home/kplee/programs/NECAT/Linux-amd64/bin/Plgd/Project.pm line 263.
Usage: necat.pl correct|assemble|bridge|config cfg_fname
    correct:     correct rawreads
    assemble:    generate contigs
    bridge:      bridge contigs
    config:      generate default config file

```













Out.

[Insert ref in markdown](https://blog.sakuragawa.moe/adding-footnotes-to-github-flavored-markdown/)








