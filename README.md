# WGApipeline

in directory with HAL alignments and directories created for each species:
find . -type d | sed '1d' | sed 's/\.\///g' | while read -r line; do parallel "halBranchMutations --refFile "$line"/{.}.rearrangements --parentFile "$line"/{.}.dups --snpFile "$line"/{.}.snp --delBreakFile "$line"/{.}.del {} "$line"" ::: *.HAL; done;
