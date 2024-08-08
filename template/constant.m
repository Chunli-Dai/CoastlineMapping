%directory
%General directory that contains the mosaic tile DEM files, such as /elev/dem/setsm/ArcticDEM/mosaic/v2.0/
macdir=[]; % Leave it blank for linux;Absolute main directory part.
tiledir=[macdir,'/u/cdai/data/ArcticDEMmosaics/']; %no longer used, no need to update.
tiledirnew=tiledir; %to store new downloaded dems
%stripdir='~/data*/ArcticDEM/region*/strips/2m/'; %directory of strip files; no longer used, no need to update.
%multidir='/u/cdai/data/lonidata//4904_2022apr19_NomeCoastlineOrthoImagery/ortho_imagery/';
multidir='/nobackupp27/cdai/orthorectwork/arcticdem_34_alaska_north/';

probdir='./output/'; % the output directory that contains the output of CoastTileMonoMain.m (probability results and image boundary results).

%control parameters
resr=2;
widm=10e3; %buffer widm of the a priori coastline, e.g., 10 km (best).
widmstat=2e3; %buffer widm of a priori coastline for calculating ndwi statistics, and the buffer of tile for calculation.
widm2=100.;   %get rid of a long thin (100 m widm) beach band for NDWI difference
cloudarea=100*100; % clouds over water area, minimum size 100m by 100m; if smaller, clouds get ignored. Also smallest island that would be kept.
lakearea=1000*500*4;%smallest water body (or lake) (m^2) that would be kept;
almt=2e3*2e3; %1e3*1e3;  %2e3*2e3;%(lost of data) %1e3*1e3 (suggested value); %minimum areas for each piece
cloudflag=1; %1 apply cloud detection; 0 do not apply cloud detection (four times faster).
novlmt=2; %if number repeats <= novlmt, set the area as edges/void. suggest value 3.
novmax=60; % if number of repeats > novmax, only select novmax of them. Suggest value 60
cntmin=25*25; %unit:pixels. The size of a priori land/ocean area should be big enough to ensure reliable statistical analysis of the histogram of the region. (Liu and Jezek, 2004)

%control parameters for multispec.m
threshold=0.5; %Suggest value 0.5 or 0.3; % general threshold of NDWI for water classification.
probthre=50; %[50., 5, 95];% threshold for water probability. %Best for migitating random coregistration offset: 50.
stdthres=0.5; % if NDWI STD > stdthres, discard the image. Suggest value 0.5
dmthres=0.6; % if mean_ocean - mean_land > dmthres, discard the image. Suggest value 0.6

% Use the following if utm; If Arctic/Antractica, update the following two manually.
projgdal='epsg:32637'; % epsg:3413 epsg:3031 epsg:32606
projstrin='UTM zone 37 north'; %'polar stereo north'; % 'polar stereo south'; 'UTM zone 45 north';

% Revise the parent directory for the following two lines.
orthworkdir=[macdir,'/u/cdai/work//orthorectwork/']; % no longer used; use job_ortho.pbs instead.

flagtide=1; %1 calculate tides; 0 no tides;
flagtidemethod=1; %1 accurate but slow; 2 fast but approximate

coregflag=0;% 0 no coregistration; Recommend 8.  8, MJ's setsm coregistration with 3 parameters.
flagrock=1; %1 use land surface as control. 0 do not use the control.

maxpxpy=15; maxpz=20; maxsigma=15;% Recommend: maxpxpy=15; maxpz=20;maxsigma=15;ArcticDEM with coregflag=8;

flagplot=0;
flagplotany=0; % Do not plot any final figures.
flagfilteradj=1;%Recommend 1. 1 apply filter in adjustOffsets.m, 0 do not apply filter (use all DEMs).

flaggage=0; %1 use gage water level, assume all tidal variations within 50 km is the same as the gage
gagefile='Nome2010to2022.csv';lat_gage=64+29.7/60;lon_gage=-(165+26.4/60)+360; %Nome gage 

currentdir=pwd;
imagesubdir=[currentdir,'/imagesubdir/'];  %directory of subset images

if ~exist(imagesubdir,'dir')
  mkdir(imagesubdir)
end

flagproj=2;
%flagproj 1, using projinv function which may produce wrong results using Matlab2017b;
            %         2, using gdaltransform

