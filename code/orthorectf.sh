#/bin/sh -f
# Use job_ortho.pbs instead. This is not good for submitting jobs with more cpus. Slow in bash command.

python /nobackupp17/elarour/coastlinedata/code/imagery_utils/pgc_ortho.py /u/cdai/elarour/coastlinedata/imagery_by_region/arcticdem_34_alaska_north/ ~/work/orthorectwork/arcticdem_34_alaska_north  --dem /nobackupp17/elarour/coastlinedata/EGM2008_Arctic.tif --epsg 3413 --outtype UInt16 --stretch ns --resample near --no-pyramids
