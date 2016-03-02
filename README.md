####Alignment of continuous coding regions for multiple genomes

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
 
In the BUSCO output directory, retrieve all of the exon annotations (other options include gene, transcript, intron, start_codon, stop_codon). Since these orthologs are under single-copy control, it is prudent to remove any Duplicated entries in the full_table_genomename file to avoid inferences drawn from sequencing artifacts.

`#!/bin/bash`  
`grep -v "Duplicated" full_table_* > good_table`  
`grep -r 'CDS' gffs/ | sed 's/gffs.*://g' | awk '$3 ~ /CDS/'  > cds`  
`awk 'NR==FNR{a[NR]=$0;next}{for (i in a){split(a[i],x," ");if ($3==x[1]&&x[4]>=$4&&x[5]<=$5)print x[1],x[2],$1,x[4],x[5],x[6],x[7],x[8],x[9] >> "cds.gff"}}' OFS='\t' cds good_table`

#### Create orthology map with [Mercator](https://github.com/hyphaltip/cndtools/tree/master/apps/mercator) by Colin Dewey

Rename exon GFF annotation to be the same as genome file (i.e. genome.fasta & genome.gff) and follow steps in link above to receive a directory of contiguous sequences for alignment. This program is particularly useful since it assembles scaffolds into larger ones if ordering of exon anchors is continuous across all genomes. Even better if one genome is of reference quality.

#### Genome alignment

In the Mercator output directory, you will find a series of numbered folders containing the orthologous sequences. These represent a much smaller percentage of the original genome, but the burden of sequencing artifacts will be further reduced by removing everything except conserved coding regions and the regions between them, validated by cross-species comparison. 

I prefer the [progressiveCactus](https://github.com/glennhickey/progressiveCactus) aligner for this task, built on the existing work of LastZ and Pecan. Even though it is version 0.0, it has already been used in [publications](https://scholar.google.com/scholar?hl=en&q=progressivecactus&btnG=&as_sdt=1%2C44&as_sdtp=) and has not thrown any errors with any of the programs in its suite. It will give you a [comprehensive list of mutations](https://github.com/glennhickey/hal/blob/master/README.md) at base pair resolution (inversions, duplications, transpositions, snps, deletions) and already has support for detecting constrained elements. 

Make sure current working directory is set to top level of Mercator sequence output directory (i.e. directory containing all numbered sequence directories). This code steps into each directory and performs the alignment, one at a time. Here, I have the Newick tree file of my species at one directory level higher.

` --out-dir=

Optional step for reformatting each HAL alignment to reorder tree with a different root (species name given to --refGenome argument)

parallel "hal2maf --inMemory --maxBlockLen 1000000 --global --refGenome tcas {} {.}.maf" ::: *.hal

parallel "maf2hal --inMemory --deflate 0 --refGenome tcas {} {.}.HAL" ::: *.maf

in directory with HAL alignments and directories created for each species:

find . -type d | sed '1d' | sed 's/\.\///g' | while read -r line; do parallel "halBranchMutations --refFile "$line"/{.}.rearrangements --parentFile "$line"/{.}.dups --snpFile "$line"/{.}.snp --delBreakFile "$line"/{.}.del {} "$line"" ::: *.HAL; done;
