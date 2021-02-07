#### Benchmarking 



- Very nice paper of benchmarking with nice gggplot2 graph on Adineta vaga paper [biorvx](https://www.biorxiv.org/content/10.1101/2020.03.16.993428v2.full.pdf)


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








