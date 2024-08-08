#!/bin/sh

# ./rerun.sh joblist60
file=$1

for dir in `cat $file`
do
	cd $dir
./prepcompile.sh
pwd
cd ../
done
