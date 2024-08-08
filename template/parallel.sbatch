#!/bin/bash

#PBS -A s5692
#PBS -W group_list=s5692 
#PBS -lselect=1:ncpus=6:model=rom_ait:mem=80GB
#PBS -l walltime=72:00:00
# #PBS -M chunlidai1@gmail.com
# #PBS -m abe
#PBS -P high

# User variables
joblistfile='/home/chunli/chunliwork/work/landslide/testparallel/joblist'
          # Joblist file (a list of subdirectories seperated by lines)
ppn=6
          # Job per node (you may not need 48. You can try and see if it makes a difference)

# Load GNU Parallel
#module load parallel/20190222/intel-19.0.5
module load pkgsrc

# Change work directory
echo PBS submit directory: $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

# Prepare joblist
joblist=(`cat $joblistfile`)

# Run parallelly
parallel --workdir $PBS_O_WORKDIR -j$ppn \
         'echo Work on subdirectory: {}; subdir1=`basename {}`; mkdir -p {//}/$subdir1; ln -sf $PWD/job.pbs {}job.pbs; cd {}; ./job.pbs > out1; ' ::: ${joblist[@]}

#Johnny's:  parallel --workdir $PBS_O_WORKDIR -j$ppn 'echo Work on subdirectory: {}; subdir1=`basename {}`; mkdir -p {}; ln -sf $PWD/job.pbs {}job.pbs; cd {}; ./job.pbs > out1;' ::: ${joblist[@]}

#to avoid error: "Segmentation fault      (core dumped) mv /u/cdai/work/runalltiles/arcticdem_34_alaska_north/$subdir1 /var/tmp/", 11/22/2022
# use cp and rm instead of mv

# Clear /var/tmp 
#cd /var/tmp
#find . -type d -user $USER -exec rm -rf "{}" +

