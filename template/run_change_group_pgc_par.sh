#/bin/sh -f
# 
# Modified from chunli@qbc.loni.org:~/chunliproject/arcticdemapp/landslide/template/run_change_group_pgc_par.sh
# Nov 16, 2022

# get tilelist from a filelist
i=0
for line in `cat tilelist`
do
i=$i+1;
tilelist[$i-1]=$line;
done

nlist=${#tilelist[@]}
echo Total number of tiles: $nlist.

echo ${tilelist[*]}

shrundir=`pwd`; #e.g. /u/sciteam/chunli/scratch/chunli/arcticdem_08_canada_baffin/
rm joblist[0-9]* qsubj*.pbs
count=0
countj=1;

for (( i=1; i<=$nlist; i++ ))
do
tilename=${tilelist[$i-1]} #39_17_2_2 "$yidc"_"$xidc"_"$xidsc"_"$yidsc"
xidc=${tilename:3:2}
yidc=${tilename:0:2}

#workdirtile1="$yidc"_"$xidc"_"$xidsc"_"$yidsc"
workdirtile1=$tilename
echo workdirtile1 $workdirtile1

unlink $workdirtile1
if [ ! -d $workdirtile1 ] 
then
    mkdir $workdirtile1
fi

# 39_17_2_2
cd $workdirtile1
pwd
#pwdsv=`pwd`;
dir_tile50km=`pwd`;

cp $shrundir/prepcompile.sh .
cp $shrundir/job.pbs .

#add compiling files steps from job.pbs here;
./prepcompile.sh

#bundle every 36 jobs. #lou can have up to 36 cpus per node.
# 36 jobs cause error: walltime exceeded limit 72 hours.
let "count+=1"
echo $count
if [[ $count -lt 10 ]] ; then
ofile=$shrundir/joblist$countj # /u/sciteam/chunli/scratch/chunli/arcticdem_08_canada_baffin/joblist1
echo $ofile
strtmp=`pwd`
echo "$strtmp"'/ ' >> "$ofile"

else #==36 #when we collect 36 jobs, we submit them
ofile=$shrundir/joblist$countj
strtmp=`pwd`
echo "$strtmp"'/ ' >> $ofile

#submit job
pbsfile=$shrundir/qsubj$countj.pbs 
cp /u/cdai/template/parallel.sbatch $pbsfile
newtext=$ofile
oldtext="/home/chunli/chunliwork/work/landslide/testparallel/joblist"
echo $oldtext $newtext  
#use # as separator in sed.
sed -i 's#'$oldtext'#'$newtext'#g' $pbsfile
lastsubtdir=`pwd` #last subtile directory of a job list.
cd $shrundir  #so that output .err .out would be in the main work directory.
qsub $pbsfile 
cd $lastsubtdir #compatible with old logic.
#reset
let "countj+=1"
count=0

fi #36


#wait if the current number of jobs in queue is >=450
while true
do
sleep 1s #wait 5 seconds
njobs=`qstat -u cdai | wc -l`
#maxjobs=450;
maxjobs=100; #40; #20 ; #30 (crashed the server); #50;
let "njobs-=5"
if [[ $njobs -lt $maxjobs ]]
then
   echo The number of current jobs is $njobs ', < '$maxjobs . 
   break
else
   echo The number of current jobs is $njobs ', wait to submit new jobs until it is <' $maxjobs .
   sleep 10m #wait 1 minute
fi
done
# wait

cd $shrundir/ #../  #$shrundir

#echo list all jobs ${jobidg[@]} #list all jobs


done # tile list for (( i=1; i<=$nlist; i++ ))


#submit the last job
#submit job
pbsfile=$shrundir/qsubj$countj.pbs 
cp /u/cdai/template/parallel.sbatch $pbsfile
newtext=$ofile
oldtext="/home/chunli/chunliwork/work/landslide/testparallel/joblist"
echo $oldtext $newtext  
#use # as separator in sed.
sed -i 's#'$oldtext'#'$newtext'#g' $pbsfile
qsub $pbsfile 

