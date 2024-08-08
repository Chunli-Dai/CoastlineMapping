Main routine: CoastTileMain.m (for stereo images), CoastTileMonoMain.m (for mono images)

For mono images:
Steps:(Same as above,except step 4)
0\ Edit the file tilelist, to include only the tiles that you want to work on.
1\ In constant.m, set the folder directory for multidir, codedir,orthworkdir, tiledirnew.
   In CoastTileMonoMain.m change codedir (the directory of source code).
2\ Skip thisif the imagery is already orthorectified: Revise multidir to the correct path for multispectral imagery in ./codec/orthorectf.sh
3\ Skip this if the imagery is already orthorectified: Orthorectificying all multispectral images: In work folder, run ./codec/orthorectf.sh
4\ In work directory, run sbatch job.slurm
	#matlab < ./codec/CoastTileMonoMain.m
5\ check output for final results
   The final output should be coastline shapefile (14_51_2_1_coast_v1.0.shp 14_51_2_1_coast_v1.0.shx), probability file (14_51_2_1_prob_v1.0.tif), and the number of repeated measurements file (14_51_2_1_nov_v1.0.tif).

Calculating coastlines using the water probability output, which should avoid the repeated computation and save computation resources.
1\ run matlab < CoastTileMonoProbMain.m
	codes: CoastTileMonoProbMain.m CoastTileMonoProb.m constant.m.

Notice the computation time and memory usage:
Case 1: 1 tile (50 km by 50 km) with 63 images (around 3 repeats at each location):
Memory: it takes 38 GB.
Computation time: 3 hours.


###################
For stereo images (not updated, not suggested due to the less amount of data)
Steps:
1\ In constant.m, assign the corresponding directory for tiledir, stripdir, multidir, codedir, orthworkdir, tiledirnew, which may contain asterisks.
   In CoastTileMain.m change codedir (the directory of source code).
2\ Revise multidir to the correct path for multispectral imagery in ./codec/orthorectf.sh
3\ Orthorectificying all multispectral images: In main folder, run ./codec/orthorectf.sh 
4\ In work directory, run matlab < ./codec/CoastTileMain.m
5\ check output for final results

