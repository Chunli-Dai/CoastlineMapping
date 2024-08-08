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
	2.1\ Copy /u/cdai/template/job_ortho.pbs to your workdir, and Revise ncpus, parallel-processes, <input_imagery_dir> <output_ortho_dir>.
	2.2\ Run command:
	     qsub job_ortho.pbs
	Which is based on:
python pgc_ortho.py <input_imagery_dir> <output_ortho_dir> --dem <dem_tif> --epsg 3413 --outtype UInt16 --stretch ns --resample near --no-pyramids
	e.g., python /nobackupp17/elarour/coastlinedata/code/imagery_utils/pgc_ortho.py /u/cdai/elarour/coastlinedata/imagery_by_region/arcticdem_31_alaska_south/ /u/cdai/work/orthorectwork/  --dem /nobackupp17/elarour/coastlinedata/EGM2008_Arctic.tif --epsg 3413 --outtype UInt16 --stretch ns --resample near --no-pyramids
	Computation time: 15042 images took 50 hours with 24 cpus, 2GB mem, 8GB virtual mem. 

3\ copy all files (tilelist CoastTileMonoMain.m  constant.m job.pbs) from template directory to your work directory, e.g. cp /u/cdai/template/* .
4\ Edit the file tilelist, to include only the tiles that you want to work on.
5\ In constant.m, set the folder directory for multidir, codedir,orthworkdir, tiledirnew.
   In CoastTileMonoMain.m change codedir (the directory of source code).

6\ In work directory, run qsub job.pbs

7\ check output for final results.
   The final output should be coastline shapefile (14_51_2_1_coast_v1.0.shp 14_51_2_1_coast_v1.0.shx), probability file (14_51_2_1_prob_v1.0.tif), and the number of repeated measurements file (14_51_2_1_nov_v1.0.tif).

Done.


#######################
Calculating coastlines using the water probability output, which should avoid the repeated computation and save computation resources.
1\ run matlab < CoastTileMonoProbMain.m
	codes: CoastTileMonoProbMain.m CoastTileMonoProb.m constant.m.

Notice the computation time and memory usage:
Case 1: 1 tile (50 km by 50 km) with 63 images (around 3 repeats at each location):
Memory: it takes 38 GB.
Computation time: 3 hours.


