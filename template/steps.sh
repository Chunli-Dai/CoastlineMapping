#/bin/sh -f

dir=arcticdem_12_canada_south_nwt
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


