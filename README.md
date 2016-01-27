# WGApipeline

Optional step for reformatting each HAL alignment to reorder alignments with a new root (species name given to --refGenome argument)

parallel "hal2maf --inMemory --maxBlockLen 1000000 --global --refGenome tcas {} {.}.maf" ::: *.hal

parallel "maf2hal --inMemory --deflate 0 --refGenome tcas {} {.}.HAL" ::: *.maf

in directory with HAL alignments and directories created for each species:

find . -type d | sed '1d' | sed 's/\.\///g' | while read -r line; do parallel "halBranchMutations --refFile "$line"/{.}.rearrangements --parentFile "$line"/{.}.dups --snpFile "$line"/{.}.snp --delBreakFile "$line"/{.}.del {} "$line"" ::: *.HAL; done;
