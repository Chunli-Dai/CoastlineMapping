% program for getting mat0.mat and arcticdem_nov.tif
% Requirements: gdal software;
% %%%% inputs needed
currentdir=pwd;
%addpath(genpath(currentdir));  %hi

% %%%% control parameters
codedir=['~/codec2/']; 
%addpath(genpath(codedir)); %hi
constant

widm0=widm+1e3; %buffer widm of the a priori coastline, e.g., 2km. need to be slightly larger than widm in Coastline.m

%Preparation: building folders
if ~exist('output','dir')
  mkdir('output')
end


%addpath(genpath(multidir));

%shpname='./GSHHS/GSHHS_f_L1.shp';% a priori coastline shapefile
if 0 % strcmp('projgdal','epsg:3031') % if 0; Antarctic %hi
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

    if strcmp(projgdal(1:7),'epsg:32') %earthdem tiles epsg:32606 
        %lat lon
	x0=-180; y0=-70;
	xe=180; ye=70;
	resrc=0.02; % 0.01 degree =1km;
    elseif strcmp(projgdal,'epsg:3413')  % arcticdem tiles, e.g., 51_08_2_2_01_02

       dx=100e3;x0=-4000e3;y0=-4000e3;%xe=3400e3;ye=4000e3; %ArcticDEM Mosaic tiles coordinate reference;
       xe=3400e3;ye=4000e3; %ArcticDEM Mosaic tiles coordinate reference;
       resrc=400;
    elseif strcmp(projgdal,'epsg:3031')  %Antarctica
        %   dx=100e3; x0=-4000e3;y0=-4000e3; %1000e3 larger than REMA tiles, e.g., 40_09_1_1
            dx=100e3; x0=-4000e3+1000e3;y0=-4000e3+1000e3;
	    xe=4000e3;ye=4000e3; 
            resrc=400;
    end

%ofile='../arcticdem_nov.tif';
ofile='../earthdem_nov.tif';
%ofile='../rema_nov.tif'; %hi
if exist(ofile,'file')
nov=readGeotiff(ofile);
nx=length(nov.x);ny=length(nov.y);
else %40m resolution
% integer type, 0 for nan.
%resrc=400;
nov.x=x0:resrc:xe;
nov.y=ye:(-resrc):y0;
nx=length(nov.x);ny=length(nov.y);
nov.z=uint16(zeros(ny,nx));
end

resrc=mean(diff(nov.x));
%M=logical(size(nov.z));

novt=nov; % for this region
novt.z=uint16(zeros(ny,nx));

for i=1:n
        XYbi=XYbg{i};
        Xb=XYbi(:,1);Yb=XYbi(:,2);

	if length(Xb)<=2|length(Xb(~isnan(Xb)))<=2|length(Yb(~isnan(Yb)))<=2
		fprintf(['\n Xb Yb bad for:',f{i},'\n'])
		continue;
	end

	if strcmp(projgdal(1:7),'epsg:32')
		%xy to lat lon
		projgdalj=projgdalg{i};
		xj=Xb;yj=Yb;
           	[latj,lonj]=xy2latlon(xj,yj,projgdalj);
		Xb=lonj;Yb=latj;
	end

        idx=round((Xb-nov.x(1))/resrc)+1;
        idy=round((Yb-nov.y(1))/(-resrc))+1;
        Mb = poly2mask(idx,idy, ny,nx); % build polygon mask       
        novt.z=novt.z + uint16(Mb);
end

%update nov with the new novt when nov is zero and novt is larger.
M=novt.z>nov.z;
nov.z(M)=novt.z(M);

%projstr='polar stereo north';
projstr=projstrin;
writeGeotiff(ofile,nov.x,nov.y,uint16(nov.z),12,255,projstr)

