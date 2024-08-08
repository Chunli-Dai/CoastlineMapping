#!/bin/sh

echo running compile.sh
echo Start compiling files
rm -f run_Tilemain.sh mccExcludedFiles.log Tilemain

module unload matlab
#module load matlab/2022b # projcrs, only works for matlab Since R2020b

#which mcc

#mcc -m CoastTileMonoMain.m -a ~/codec2/ -a constant.m 
#/home4/cdai/workpfe/software/matlab2020bbin/bin/mcc -m CoastTileMonoMain.m -a ~/codec2/ -a constant.m

/swbuild/cdai/software/matlab2017bbin/R2017b/bin/mcc -m CoastTileMonoMain.m -a ~/codec2/ -a constant.m

mv CoastTileMonoMain Tilemain
cp ~/template/run_Tilemain.sh .

echo Compiling files end.
