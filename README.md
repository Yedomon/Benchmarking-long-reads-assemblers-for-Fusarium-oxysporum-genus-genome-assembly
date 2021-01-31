#### Benchmarking 

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








