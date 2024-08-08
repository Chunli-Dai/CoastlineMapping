#!/bin/sh


echo work directory `pwd`

coderundir=`pwd`; #'../'; # Code with updated parameters for a specific job.
#arcticdem_08_canada_baffin/21_30_2_2/
str1=$(basename $coderundir) #21_30_2_2
workdir1=$str1
coderundir=$(dirname $coderundir)

echo main directory [Code with updated parameters for a specific job] for this region: $coderundir

#strtile=${str1:0:9}; # "$yidc"_"$xidc"_"$xidsc"_"$yidsc"
workdirtile=`pwd` # e.g. /mnt/c/scratch/sciteam/GS_bazu/user/chunli/arcticdem_08_canada_baffin/21_30_2_2
pwd

#55_16_2_1_coast_tide_thre50_v1.0.shp

ln -fs $coderundir/mat0.mat .
cp  $coderundir/*.m .

#Update projection in constant.m 
#Update projection information in constant.m based on tile names (e.g., utm10n_47_04_2_2).
subtile=$workdir1;
projgdal=$(awk -F"'" '/projgdal/{print $2}' constant.m)
str1=${subtile:0:3} # e.g.,'utm'
if [ "$str1" = "utm" ]; then
    echo "UTM projection"

    cp  $coderundir/compile.sh .
    
    utmid=${subtile:3:3}; #e.g.,10n
    #sed -i "s|37N|$utmid|g" constant.m
    nsstr=${subtile:5:1}; # 'n' or 's'
    #Northern hemisphere, e.g., UTM Zone 3N: EPSG 32603;
    if [[ "${nsstr,,}" == 'n' ]]; then
    sed -i "s|epsg:32637|epsg:326${utmid:0:2}|g" constant.m
    sed -i "s|UTM zone 37 north|UTM zone ${utmid:0:2} north|g" constant.m
     #Southern hemisphere, e.g., UTM Zone 2S: EPSG 32702;
    elif [[ "${nsstr,,}" == 's' ]]; then
    sed -i "s|epsg:32637|epsg:327${utmid:0:2}|g" constant.m
    sed -i "s|UTM zone 37 north|UTM zone ${utmid:0:2} south|g" constant.m
    fi

    #matlab license is blocked in PBS
    #compile
./compile.sh

elif [ "$projgdal" = "epsg:3413" ]; then  
    echo "Polar Stereographic North projection"

ln -fs $coderundir/run_Tilemain.sh     # Copy compiled Matlab code.  Jason Li @ 4/13/2021
ln -fs $coderundir/Tilemain            # Copy compiled Matlab code.  Jason Li @ 4/13/2021

#    sed -i "s|projgdal='epsg:32637';|projgdal='epsg:3413;'|g" constant.m
    # projstrin update
#   sed -i "s|UTM zone 37 north|polar stereo north|g" constant.m

elif [[ "${projgdal}" == 'epsg:3031' ]]; then
    echo "Antarctic Polar Stereographic projection, epsg:3031"

ln -fs $coderundir/run_Tilemain.sh     # Copy compiled Matlab code.  Jason Li @ 4/13/2021
ln -fs $coderundir/Tilemain            # Copy compiled Matlab code.  Jason Li @ 4/13/2021

#   sed -i "s|UTM zone 37 north|polar stereo south|g" constant.m

else
        echo Projection is not as expected. Double check projgdal= $projgdal projstrin in constant.m
fi

#ln -fs $coderundir/run_Tilemain.sh     # Copy compiled Matlab code.  Jason Li @ 4/13/2021
#ln -fs $coderundir/Tilemain            # Copy compiled Matlab code.  Jason Li @ 4/13/2021

echo 1 > input.txt
echo $workdir1 >> input.txt
#matlab -nodisplay -r "addpath('/home/howat.4/mosaicfunctions'),batch_initializeMosaic(firstTile,lastTile), exit"


# your logic ends here
END=$(date +%s)
echo End time: $END ; # test if job still stalls.

#DIFF=$(( $END - $START ))  #seconds
#DIFFhr=`awk "BEGIN {printf \"%.2f\n\", $DIFF/3600}"`
#echo job.pbs total run time for 1 tile within $str1 is $DIFFhr hours.


