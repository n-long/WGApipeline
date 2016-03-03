###Alignment of continuous coding regions across multiple genomes

####Get exon-level coding annotations (BUSCO)

[BUSCO](busco.ezlab.org/) is the [CEGMA](http://korflab.ucdavis.edu/datasets/cegma/) successor for curated sets of single-copy orthologs for specific clades/lineages.  
Requires Python 3 -- [how to keep separate Python3 installation](http://askubuntu.com/questions/244544/how-do-i-install-python-3-3/290283)  

`py3  BUSCO_v1.1b1.py -g genomeassemby.fasta -m all -l yourlineage -sp closest_species_for_augustus -c 20 -o output_folder`  
* py3 -- symlink for Python3  
* -g -- genome assembly FASTA  
* -m -- perform all analysis steps  
* -l -- BUSCO lineage  
* -sp -- [AUGUSTUS species](http://augustus.gobics.de/binaries/README.TXT)  
* -c -- number of cores  
* -o -- output directory
 
In the BUSCO output directory, retrieve all of the exon annotations with the following (other options include gene, transcript, intron, start_codon, stop_codon). Since these orthologs are under single-copy control, it is prudent to remove any Duplicated entries in the full_table_*  file to avoid inferences drawn from sequencing artifacts.

`#!/bin/bash`  
`cd output_folder`  
`grep -v "Duplicated" full_table_* > good_table`  
`grep -r 'CDS' gffs/ | sed 's/gffs.*://g' | awk '$3 ~ /CDS/'  > cds`  
`awk 'NR==FNR{a[NR]=$0;next}{for (i in a){split(a[i],x," ");if ($3==x[1]&&x[4]>=$4&&x[5]<=$5)print x[1],x[2],$1,x[4],x[5],x[6],x[7],x[8],x[9] >> "cds.gff"}}' OFS='\t' cds good_table`

#### Create orthology map with [Mercator](https://github.com/hyphaltip/cndtools/tree/master/apps/mercator)

Rename exon GFF annotation to be the same as genome file (i.e. genome.fasta & genome.gff) and follow steps in link above to receive a directory of contiguous sequences for alignment. This program is particularly useful since it assembles scaffolds into larger ones if ordering of exon anchors is continuous across all genomes. Even better if one genome is of reference quality.

#### Genome alignment

In the Mercator output directory, you will find a series of numbered folders containing the orthologous sequences. These represent a much smaller percentage of the original genome, but the burden of sequencing artifacts will be further reduced by removing everything except conserved coding regions and the regions between them, validated by cross-species comparison. 

I prefer the [progressiveCactus](https://github.com/glennhickey/progressiveCactus) aligner for this task, built on the existing work of LastZ and Pecan. Even though it is version 0.0, it has already been used in [publications](https://scholar.google.com/scholar?hl=en&q=progressivecactus&btnG=&as_sdt=1%2C44&as_sdtp=) and has not thrown any errors with any of the programs in its suite. It will give you a [comprehensive list of mutations](https://github.com/glennhickey/hal/blob/master/README.md) at base pair resolution and includes support for detecting constrained elements (PhyloP). The rest of this document assumes you are using progressiveCactus and  [HALtools](https://github.com/glennhickey/hal/blob/master/README.md) and have both included in your shell path. 

Make sure the current working directory is set to the top level of Mercator sequence output directory (i.e. directory containing all numbered sequence subdirectories). This code steps into each subdirectory and performs the alignment, one at a time. I keep the Newick tree file in the same directory that we start in, which contains all of the sequence information for the aligner. `workdir` is used for temp files, and `alignment.hal` will hold the alignment output.

`find "$PWD" -type d | sed 1d | while read -r line; do cd "$line" && echo "$line" && sh /home/hdd/4/progressiveCactus/bin/runProgressiveCactus.sh --overwrite --maxThreads=22 ../treefile workdir/ alignment.hal; done`

Having each alignment sequestered away in a subdirectory will not be convenient for the repeated analyses you will undoubtedly perform, so let's pull them all out into a new directory while maintaining the numbered folder name (this is important for retrieving the sequence coordinates).

`for subdir in *; do cp $subdir/alignment.hal ../HALtime/$subdir.hal;done`

Here are some optional commands for reformatting each HAL alignment to reorder the species tree with a different root (species name given to --refGenome argument)

`parallel "hal2maf --inMemory --maxBlockLen 1000000 --global --refGenome tcas {} {.}.maf" ::: *.hal`

`parallel "maf2hal --inMemory --deflate 0 --refGenome tcas {} {.}.HAL" ::: *.maf`

Let's create a subdirectory for each species inside the directory containing all of the .hal alignment files. Make sure they match the species name specified in the Newick tree file given to progressiveCactus. This command assumes you don't have any other subdirectories, otherwise you will need to modify the sed command or include another piped `grep -v` so that only the species subdirectories show up before `while read -r line` is executed. Make sure the `::: *.HAL` section matches the extension of your HAL alignment files. 

`mkdir species1/ species2/ species3/ etc...`

`find . -type d | sed '1d' | sed 's/\.\///g' | while read -r line; do parallel "halBranchMutations --refFile "$line"/{.}.rearrangements --parentFile "$line"/{.}.dups --snpFile "$line"/{.}.snp --delBreakFile "$line"/{.}.del {} "$line"" ::: *.HAL; done;`

Now each species subdirectory is populated with files containing inversions, duplications, transpositions, snps, and deletions that are numbered by the matching HAL alignment. NOTE: duplications (.dups) are only given with respect to the root of the species tree. 

The HAL mutation coordinates are given only in relation to the start and end of each Mercator segment. Therefore, in your Mercator directory, you will need to extract the coordinates of each numbered Mercator segment out of the `map` file by using the `genomes` file as reference 

My `genomes` file looks like 

`tcas    tmad    tconf   tfree   agla    dendro`

So I can extract the coordinates for each species as follows:

`cut -f1-5 map > tcas.coord; cut -f1,6-9 map > tmad.coord; cut -f1,10-13 map > tconf.coord; cut -f1,14-17 map > tfree.coord; cut -f1,18-21 map > agla.coord; cut -f1,22-25 map > dendro.coord`

We can intersect the HAL mutation coordinates with Mercator coordinates as follows (for greater ease you could first concatenate all genome-specific mutations, except for duplications for reason listed above):

awk 'NR==FNR{a[NR]=$0;next}{for (i in a){split(a[i],x," ");if ($1==x[1]) print x[2],x[3]+$3,x[3]+$4,$5,x[5] > "tmad.rearrangements.bed"}}' tmad.coord tmad.mutations

The logic operator `($1==x[1])` requires column 1 of file2 ($1 for tmad.mutations) and column 1 of file1 (x[1] for tmad.coord) to match, followed by printing of the specified columns. x[] is used for columns of file1 while $ specifies columns of file2. 

With the transformed coordinates, you can easily intersect these files with annotations from other sources using [`bedtools intersect`](http://bedtools.readthedocs.org/en/latest/content/tools/intersect.html)
