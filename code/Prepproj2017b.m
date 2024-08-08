
% File preparation: To use old matlab 2017b projfwd, the new projcrs is
% converted to compatible struct and saved to files. Then the struct can be
% loaded and used for matlab2017b functions.

% Function run in after R2020b

% South zones 327, 01 to 60;
% North zones 326, 01 to 60;
% UTM Zone 2S: EPSG 32702; UTM Zone 3N: EPSG 32603;

codedir=['~/codec2/']; 
addpath(genpath(codedir),'-end');

%
% a simple geotiff file in polar stereo north
data.x=-2:2; data.y=-2:2; data.z=ones(length(data.y),length(data.x));

projstr='polar stereo north'; %'epsg:3413'
writeGeotiff('~/codec2/sample/sample.tif',data.x,data.y,double(data.z),5,0,projstr)

for icase=1:2
    if icase==1;c0=32600; %north
    elseif icase==2;c0=32700;  %south
    end
for zone=1:60
       crs_epsg=c0+zone; %e.g.32636;

%      crs_epsg=32636;%e.g.32636;
%        p1 = projcrs(crs_epsg); %projcrs, Since R2020b
       
       % projcrs to struct; too trivial
       
       % write sample geotiff files
       str=['gdalwarp -of GTiff -t_srs EPSG:',num2str(crs_epsg),' ~/codec2/sample/sample.tif ~/codec2/sample/sample',num2str(crs_epsg),'.tif'];
       [status , cmdout]=system(str);
       
       
end
       
end %icase
       

% test - correct
crs_epsg=32636;
samplefile=['~/codec2/sample/sample',num2str(crs_epsg),'.tif'];
proj = geotiffinfo(samplefile);
x=747657;y=4031478;
[lat, lon] = projinv(proj, x, y);
num2str([lon lat])
lon= 35.7614221046142; lat = 36.396568533103;
[x2, y2] = projfwd(proj, lat, lon);
num2str([x2, y2])

crs_epsg=32637;
samplefile=['~/codec2/sample/sample',num2str(crs_epsg),'.tif'];
proj = geotiffinfo(samplefile);
[x2, y2] = projfwd(proj, lat, lon);
num2str([x2, y2])

%Argentina
lat=-54.749747;lon= -67.241607;
crs_epsg=32719;
samplefile=['~/codec2/sample/sample',num2str(crs_epsg),'.tif'];
proj = geotiffinfo(samplefile);
[x2, y2] = projfwd(proj, lat, lon);
num2str([x2, y2])
        
%Turkey
% x=747657;y=4031478;
%[lat,lon]=xy2latlon(x,y,'epsg:32636'); 
%[lon lat]=35.761422104614205 36.396568533103043 
% echo 747657 4031478 | gdaltransform -s_srs EPSG:32636 -t_srs EPSG:4326
% lon lat = 35.7614221046142 36.396568533103 0
% [x2,y2]=latlon2xy(lat,lon,projgdalj);
%[x2, y2]= 747657 4031478
% convert to another zone 
%[x2,y2]=latlon2xy(lat,lon,'epsg:32637');
% x2 y2=209536.89578      4032809.1721;

%Argentina 
%lat=-54.749747;lon= -67.241607;
% [x2,y2]=latlon2xy(lat,lon,'epsg:32719'); 
% [x2, y2]='613174.89557      3931637.4999'; 
