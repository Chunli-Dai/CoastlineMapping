#/bin/sh -f

#mapformats was moved out of codec2

tilelist=('arcticdem_01_iceland' 'arcticdem_02_greenland_southeast' 'arcticdem_03_greenland_southwest' 'arcticdem_04_greenland_central' 'arcticdem_05_greenland_northeast' 'arcticdem_06_greenland_northwest' 'arcticdem_07_canada_ellesmere' 'arcticdem_08_canada_baffin' 'arcticdem_09_canada_victoria' 'arcticdem_10_canada_north_mainland' 'arcticdem_11_canada_north_hudson' 'arcticdem_12_canada_south_nwt' 'arcticdem_14_svalbard' 'arcticdem_15_russia_novaya_zemlya' 'arcticdem_18_russia_cherskly' 'arcticdem_19_russia_magadanskaya' 'arcticdem_20_russia_kamchatka' 'arcticdem_21_russia_yakutiya_east' 'arcticdem_22_russia_central_east' 'arcticdem_23_russia_yakutiya_west' 'arcticdem_24_russia_central_west' 'arcticdem_25_russia_norilsk' 'arcticdem_26_russia_petersburg' 'arcticdem_27_russia_murmansk' 'arcticdem_28_scandinavia' 'arcticdem_29_russia_franz_josef' 'arcticdem_30_russia_siberian_islands' 'arcticdem_31_alaska_south' 'arcticdem_34_alaska_north')

nlist=${#tilelist[@]}
echo Total number of tiles: $nlist.

echo ${tilelist[*]}

#check ~/codec2/readme.txt for steps before this.

for (( i=14; i<=$nlist; i++ ))
do
dir=${tilelist[$i-1]}
echo $i $dir
cd /u/cdai/work/runalltiles/
mkdir $dir

#dir=arcticdem_12_canada_south_nwt
cd $dir
cp /u/cdai/template/* .   # warning: need to change the directory
#cp ../rema_07_qml_ne/*.m . ; sed -i 's/rema_07_qml_ne/'$dir'/g' constant.m
cp ../arcticdem_34_alaska_north/*.m .
sed -i 's/arcticdem_34_alaska_north/'$dir'/g' constant.m #change directory
#grep $dir ademtiles_join-ademreg.csv > df
grep $dir ademtiles_join-ademreg_arcticdem2.csv > df
grep $dir ademtiles_join-ademreg_rema1.csv > df
grep $dir ademtiles_join-ademreg_rema2.csv > df #01 05 
grep $dir ademtiles_join-ademreg_earthdem.csv > df
sed -i 's/'$dir',/''/g' df 
mv df tilelist  #get tilelist
#Get mat0.mat
#/home4/cdai/workpfe/software/matlab2017bbin/R2017b/bin/mcc -m Tilemain_nov.m -a ~/codec2/ -a constant.m
#cp ~/template/run_Tilemain_nov.sh .

vi Tilemain_nov.m (update two lines in Tilemain_nov.m); update projgdal projstrin in constant.m
compile the lines within job_nov.pbs
qsub job_nov.pbs #run Tilemain_nov.m to update the arcticdem_nov.tif and get mat0.mat.

# for arcticdem or rema do this ; skip if earthdem
./compile.sh #to compile all matlab codes. Note: always load matlab/2020a version. Has to re-compile it for each region (not sure why)! If too slow (>10 minutes), try nohup ./compile.sh > outcompile &
#cp ~/template/run_Tilemain.sh . #add LD_LIBRARY_PATH to run_Tilemain.sh 

#bundle all jobs
nohup ./run_change_group_pgc_par.sh > outrun1 &

cd ../
#exit
done

