#### An attempt at avoiding sequencing artifacts when doing comparative genomics on draft genome assemblies

####Get exon-level coding annotations (BUSCO)

[BUSCO](busco.ezlab.org/) is the CEGMA successor for curated sets of single-copy orthologs for specific clades/lineages.  
Requires Python 3 -- [how to keep separate Python3 installation](http://askubuntu.com/questions/244544/how-do-i-install-python-3-3/290283)  

`py3  BUSCO_v1.1b1.py -g genomeassemby.fasta -m all -l yourlineage -sp closest_species_for_augustus -c 20 -o output_folder`  
* py3 -- symlink for Python3  
* -g -- genome assembly FASTA  
* -m -- perform all analysis steps  
* -l -- BUSCO lineage  
* -sp -- [AUGUSTUS species](http://augustus.gobics.de/binaries/README.TXT)  
* -c -- number of cores
 
In the BUSCO output directory, retrieve all of the exon annotations (other options include gene, transcript, intron, start_codon, stop_codon) 

`#!/bin/bash`  
`grep -r 'CDS' gffs/ | sed 's/gffs.*://g' | awk '$3 ~ /CDS/'  > cds`  
`awk 'NR==FNR{a[NR]=$0;next}{for (i in a){split(a[i],x," ");if ($3==x[1]&&x[4]>=$4&&x[5]<=$5)print x[1],x[2],$1,x[4],x[5],x[6],x[7],x[8],x[9] >> "cds.gff"}}' OFS='\t' cds full_table_*`

#### Create orthology map with [Mercator](https://github.com/hyphaltip/cndtools/tree/master/apps/mercator) by Colin Dewey

Rename exon GFF annotation to be the same as genome file (i.e. genome.fasta & genome.gff) and follow steps in link above to receive a directory of contiguous sequences for alignment. This program is particularly useful since it assembles scaffolds into larger ones if ordering of exon anchors is continuous across all genomes. Even better if one genome is of reference quality.

#### Genome alignment

In the Mercator output directory, you will find a series of numbered folders 

Optional step for reformatting each HAL alignment to reorder tree with a different root (species name given to --refGenome argument)

parallel "hal2maf --inMemory --maxBlockLen 1000000 --global --refGenome tcas {} {.}.maf" ::: *.hal

parallel "maf2hal --inMemory --deflate 0 --refGenome tcas {} {.}.HAL" ::: *.maf

in directory with HAL alignments and directories created for each species:

find . -type d | sed '1d' | sed 's/\.\///g' | while read -r line; do parallel "halBranchMutations --refFile "$line"/{.}.rearrangements --parentFile "$line"/{.}.dups --snpFile "$line"/{.}.snp --delBreakFile "$line"/{.}.del {} "$line"" ::: *.HAL; done;
