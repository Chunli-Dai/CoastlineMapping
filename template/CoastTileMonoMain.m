% main program for getting coastline for each ArcticDEM tile
% Requirements: gdal software;ogr2ogr
% %%%% inputs needed
currentdir=pwd;
%addpath(genpath(currentdir));

% %%%% control parameters
codedir=['~/codec2/']; 
%addpath(genpath(codedir),'-end');
constant

widm0=widm+1e3; %buffer widm of the a priori coastline, e.g., 2km. need to be slightly larger than widm in Coastline.m

%Preparation: building folders
if ~exist('output','dir')
  mkdir('output')
end

if ~exist('mat0.mat','file') %readding boundary; time 1 hour for 90751 strip files and 10130 mono xml files

%shpname='./GSHHS/GSHHS_f_L1.shp';% a priori coastline shapefile
% https://www.soest.hawaii.edu/pwessel/gshhg/ 
if strcmp('projgdal','epsg:3031') % if 0; Antarctic
	% Version 2.3.7 Released June 15, 2017
	[status , cmdout ]=system(['find ',codedir,' -name GSHHS_f_L5.shp']);
else %ArcticDEM or EarthDEM or rema_01_subantarctic_islands rema_05_ellsworth_land
[status , cmdout ]=system(['find ',codedir,' -name GSHHS_f_L1.shp']);
end
shpname=deblank(cmdout);

% %%% Preparation: get the list of xml files which contain boundries
filename='monolist'; %'boundaries_reg31.dat';
if ~exist(filename,'file')
   str=sprintf('find  %s -name *[0-9].xml > %s',deblank(multidir),filename);
  [status, cmdout]=system(str);
end

fprintf ('\n Step 0: geting the boundary for all files in the region.')
%READ INPUT PARAMETERS; getting the boundaries for all files
fid = fopen(filename);
n = linecount(fid);
fid = fopen(filename);
% range=fscanf(fid, '%f', [4, n))';
f=cell(n,1);fdir=cell(n,1);
range=zeros(n,4);XYbg=cell(n,1); projgdalg=cell(n,1);
%exclude panchromatic bands; may be included later.
idd=[];
for i=1:n
   ifile=fgetl(fid);
   [demdir,name,ext] =fileparts([strtrim(ifile)]);
   f{i}=[name,ext];
   fdir{i}=[demdir,'/'];%working on two regions 
   satname=f{i}(1:4);

   % get the boundary from xml file
   [XYbi,rangei,projgdali]=imagebd(ifile);
   range(i,1:4)=rangei;XYbg{i}=XYbi;projgdalg{i}=projgdali;
    if strcmp(satname,'WV01')||strcmp(satname,'GE01')
        idd=[idd;i];
    end
end
range(idd,:)=[];f(idd)=[];fdir(idd)=[];XYbg(idd)=[];
projgdalg(idd)=[];
display(['demdir=',demdir])

save mat0.mat -v7.3

else 
load mat0.mat
constant %update projgdal
end


dx=100e3;x0=-4000e3;y0=-4000e3;
%xe=3400e3;ye=4000e3; %ArcticDEM Mosaic tiles coordinate reference;

% ArcticDEM mosaic tile grids
%32_33_2_2_5m_v2.0 ; 
%name convention; yid_xid_xids_yids
% inputtype=1;
%get the input from input.txt
filename='input.txt';
fid = fopen(filename);
inputtype=fscanf(fid, '%d', [1, 1])';
if inputtype ==1 %use this; ignore other options.
    str=fgetl(fid);tilefile=fgetl(fid);
    rang0=getbox(tilefile);
elseif inputtype==2 || inputtype==3
   latlon=fscanf(fid, '%f', [2, 1])';
   lateq=latlon(1);loneq=latlon(2);
elseif inputtype ==4
   rang0=fscanf(fid, '%f', [4, 1])';
end

%To be compatible with older versions. %55_16_2_1_2m_v3.0_reg_dem.tif
    %tilefile=[tilefile,'_2m_v3.0_reg_dem.tif'];

    % Find whether this tile contains any coastline.
    % Buffer the tile boundary by widm;
    %rang0=[x-widm0 x+dx/2+widm0 y-widm0 y+dx/2+widm0];
    rang0=[rang0(1)-widm0 rang0(2)+widm0 rang0(3)-widm0 rang0(4)+widm0];
    x0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];y0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];
    % x y to lat lon
    [lat0,lon0]=xy2latlon(x0,y0,projgdal);

    if strcmp('projgdal','epsg:3031') % if 0; Antarctic; check all points in coastline, such as GSHHS_f_L5.shp

	Sc=shaperead(shpname);
	
	cnt=0;
	for j=1:length(Sc) %330 polygons
		    latj=Sc(j).Y; lonj=Sc(j).X;
                    [sx,sy]=latlon2xy(latj,lonj,'epsg:3031');
                    isIntersecting = inpolygon(sx, sy, x0, y0);
        	    cnt=cnt+sum(isIntersecting(:));
        end
    else % arcticdem or earthdem; crop coastline shapefile, e.g., GSHHS_f_L1.shp

    bb = geoshape(lat0,lon0,'Geometry','Polygon');
    tileshape=sprintf('output/%s_tile.shp',tilefile);
    tilecoastname=sprintf('output/%s_tilegshhs.shp',tilefile);

%   tileshape=sprintf('output/%02d_%02d_%01d_%01d_tile.shp',yid,xid,xids,yids);
%   tilecoastname=sprintf('output/%02d_%02d_%01d_%01d_tilegshhs.shp',yid,xid,xids,yids);

    shapewrite(bb,tileshape);
    system(['rm ',tilecoastname])
    system(['time ogr2ogr -overwrite -clipsrc ',tileshape,' ',tilecoastname,' ',shpname]);
    %ogr2ogr -overwrite -clipsrc tile.shp tilegshhs.shp GSHHS/GSHHS_f_L1.shp
    S = shaperead(tilecoastname);
    cnt=length(S); %figure;mapshow(S);

    %incase the tile is inside the coastline polygon, instead of coastline crossing the tile.
    if cnt~=0
    Sbox=shaperead(tileshape);
    Sbox.area=polyarea( Sbox.X(~isnan(Sbox.X)) , Sbox.Y(~isnan(Sbox.X)));
    S(1).area2=polyarea( S(1).X(~isnan(S(1).X)) , S(1).Y(~isnan(S(1).X)));
    if abs(S(1).area2 - Sbox.area)<1e-8  %1e-4 degree (11m)
       cnt=0;
    end
    end %cnt

    end %projgdal

    if cnt==0; fprintf(['This tile covers no coastline ',tilefile]);exit; end
    
fprintf (['\n Working on tile:',tilefile,'; \n'])

	tic
    [Co]=CoastTileMono(tilefile, S,range,XYbg,projgdalg,f,fdir);
%   tileshape=sprintf('output/%02d_%02d_%01d_%01d_tile*',yid,xid,xids,yids);
	toc
	    system(['rm ',tileshape]);

return
