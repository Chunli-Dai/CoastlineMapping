Main routine: CoastTileMonoMain.m (for mono images)

For mono images:
Steps:
1\ (do this step one time only): add path (for SETSM, and imagery_utils) to your .bashrc file
export PATH=$PATH:/u/cdai/landslide/code1/SETSM:/nobackupp17/elarour/coastlinedata/code/imagery_utils/
#to avoid this error: error while loading shared libraries: libgeotiff.so.2: cannot open shared object file: No such file or directory
export LD_LIBRARY_PATH=/u/cdai/landslide/code1/SETSM/tools/tiff-4.0.3/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/u/cdai/landslide/code1/SETSM/tools/libgeotiff-1.4.2/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/u/cdai/landslide/code1/SETSM/tools/proj-5.1.0/lib:$LD_LIBRARY_PATH

2\ (do this step once per region): orthorectify imagery; Code by Erik at PGC (more details see notes_erik.txt).
	2.1\ Copy /u/cdai/template/job_ortho.pbs to workdir (e.g., /home4/cdai/workpfe/orthorectwork/), and Revise ncpus, parallel-processes, <input_imagery_dir> <output_ortho_dir>.
	     cd wordir
	     dir=arcticdem_03_greenland_southwest
	     cp ~/template/job_ortho.pbs .; sed -i 's/arcticdem_02_greenland_southeast/'$dir'/g' job_ortho.pbs  #change directory
	2.2\ Run command:
	     qsub job_ortho.pbs
	Which is based on:
python pgc_ortho.py <input_imagery_dir> <output_ortho_dir> --dem <dem_tif> --epsg 3413 --outtype UInt16 --stretch ns --resample near --no-pyramids
	e.g., python /nobackupp17/elarour/coastlinedata/code/imagery_utils/pgc_ortho.py /u/cdai/elarour/coastlinedata/imagery_by_region/arcticdem_31_alaska_south/ /u/cdai/work/orthorectwork/  --dem /nobackupp17/elarour/coastlinedata/EGM2008_Arctic.tif --epsg 3413 --outtype UInt16 --stretch ns --resample near --no-pyramids
	Computation time: 15042 images took 50 hours with 24 cpus, 2GB mem, 8GB virtual mem. 
	                   8594 images took 22 hours with 24 cpus. 
		          26812 images took 17 hours with 64 cpus.

#Carried out the following steps for each region:
3\ follow steps in /u/cdai/template/steps.sh

4\ check output for final results.
   The final output should be coastline shapefile (14_51_2_1_coast_v1.0.shp 14_51_2_1_coast_v1.0.shx), probability file (14_51_2_1_prob_v1.0.tif), and the number of repeated measurements file (14_51_2_1_nov_v1.0.tif).

Done.


#######################

Notice the computation time and memory usage:
Case 1: 1 tile (50 km by 50 km) with 63 images (around 3 repeats at each location):
Memory: it takes 38 GB.
Computation time: 3 hours.


