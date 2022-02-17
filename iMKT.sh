#!/bin/bash -l

## bash script for running iMKT analysis

dataPath="$1"
multiFASTA="$2"

fastaPrefix4fold="$(cut -d'.' -f1 <<<"$multiFASTA")"_'4fold'
fastaPrefix2fold="$(cut -d'.' -f1 <<<"$multiFASTA")"_'2fold'

# echo $fastaPrefix

conda activate imktData
cd /bscb/bscb10/rs2474/messer/softwares/iMKTData-master/src

echo "MKT on 0fold + 4fold sites: "
./sfsFromFasta_runxi_4fold.py $dataPath $multiFASTA $fastaPrefix4fold.daf $fastaPrefix4fold.div
Rscript --vanilla run_iMKT.R $dataPath $fastaPrefix4fold.daf $fastaPrefix4fold.div

echo "MKT on 0fold + 2fold + 4fold sites all together: "
./sfsFromFasta_runxi_2fold.py $dataPath $multiFASTA $fastaPrefix2fold.daf $fastaPrefix2fold.div
Rscript --vanilla run_iMKT.R $dataPath $fastaPrefix2fold.daf $fastaPrefix2fold.div
