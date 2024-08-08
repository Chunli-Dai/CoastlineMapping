function [Co]=CoastTileMono(tilefile,S,range,XYbg,projgdalg,f,fdir)
% coastline from orthoimages;mosaicking of coastlines
% Requirements: 
%               
% Versions:
% CoastTilev2.m: to save space: do not save all groups, just save the data that has the maximum repeats, and disgrad others.
%  	July 2018: collect the edge line of each images, and set the a priori coastline buffer zone as edge.
%		   and the edge of overlapping groups.
%		   Crop the area of coastal band to save space.
% 

% Chunli Dai chunlidai1@gmail.com
% March 2018
% 

% close all
% clear 
% clc


params.I = [1 0 0 0 0 0 0];
params.C = [50 0.001 0.001 0.05 0.0001];
%         params.G = [3000 20];
params.M = '3D';
params.V = [10 20 10];
params.G = [9000 20]; %Adjust the max height parameter for the 2012 Kamchatka Volcano

%load and set up the parameter
constant

%flagplot=1; %use the value in constant.m

res='2';
% demdir=dir(['/data2/ArcticDEM/region_',regionnum,'*']);
macdir='/Users/chunlidai/surge/';
% macdir='/Users/chunli/surge/';
macdir=[];

% addpath(genpath([macdir,'/../Box Sync/ISSMCoastline/']));
Co=[];

yr=365.25;
neq=1;
dsr=0.2;%0.04; %200m ; %0.2; % 8m to 40m
nsr=1./dsr;
%resr=2.; %str2double(res)/dsr;
% resr=40.; 
resrc=40.; %for coregisteration
% res=mt.info.map_info.dx;
dsr=res/resr;
% demdir=[macdir,'/data2/ArcticDEM/region_08_canada_baffin/tif_results/8m/'];
%coregflag=7;%1 parameter (vertical), 3 parameter (Ian's); 7 parameter (MJ's)

%Grids of ArcticDEM Tiles
% 54_06_2_2_5m_v2 yid_xid_xids_yids
% 1,2 ; 2,2
% 1,1 ; 2,1
%x=x0+(xid-1)*dx; %left edge of a box;
dx=100e3;x0=-4000e3;y0=-4000e3;%xe=3400e3;ye=4000e3; %ArcticDEM Mosaic tiles coordinate reference;

% %%% Preparation: get rang0 from tile name
rang0=getbox(tilefile); %no buffer
rang0b=rang0;
%tile boundary with buffer widm for selecting files and grouping; 
%need to be the buffer with of a priori incase of boundary being eroded.
rang0=[rang0(1)-widm rang0(2)+widm rang0(3)-widm rang0(4)+widm];
ranget=round(rang0/resr)*resr;rang0=ranget;

rang1=[rang0(1)-widmstat rang0(2)+widmstat rang0(3)-widmstat rang0(4)+widmstat];
ranget=round(rang1/resr)*resr;rang1=ranget;

%54_06_2_2_5m_v2.0_reg_dem.tif  54_06_2_2_coast_v1.0.tif 54_06_2_2_prob_v1.0.tif  54_06_2_2_lessdata_v1.0.tif 
% 43_51_1_1_coast_tide_thre50_v1.0.shp
ifilestr=[tilefile,'_']; % utm49s_42_07_2_1_  54_06_2_2_
%ofile1=['output/',ifile(1:10),'coast_v1.0.shp'];
ofile1=['output/',ifilestr,'coast_v1.0.shp'];
ofile2=['output/',ifilestr,'prob_v1.0.tif'];
%ofile3=[ifile(1:10),'lessdata_v1.0.tif'];
ofile3=['output/',ifilestr,'nov_v1.0.tif']; %number of overlapping coverages
ofile4=['output/',ifilestr,'bound_v1.0.tif']; %image boundaries, no data areas

% rang0 = [  -3450000    -3400000     1350000     1400000];

xeq=(rang0(1)+rang0(2))/2;yeq=(rang0(3)+rang0(4))/2;
[lateq,loneq]=xy2latlon(xeq,yeq,projgdal);
%[lateq,loneq]=polarstereo_inv(xeq,yeq,[], [],70,-45);
formatSpec = '%6.1f';

% rang0=[-3467 -3455 110 124 ]*1e3;
x0=[rang1(1) rang1(2) rang1(2) rang1(1) rang1(1) ];y0=[rang1(4) rang1(4) rang1(3) rang1(3) rang1(4) ];
[lat0,lon0]=xy2latlon(x0,y0,projgdal);
%[lat0,lon0]=polarstereo_inv(x0,y0,[], [],70,-45);
tx=rang1(1):resr:rang1(2);ty=rang1(4):-resr:rang1(3);
xout=tx;yout=ty; %needs a buffer zone of the box avoiding excluding small water body at boundaries.

% fast, using gdal for cropping shapefiles.
%load the coastline Hawii see coastlineArcticDEMv2.m
cnt=length(S);

% %%% Preparation: get water mask from a priori coastline
xw=rang0(1):resrc:rang0(2);yw=rang0(4):-resrc:rang0(3); % a priori water mask
nwx=length(xw);
% method 1 %poly2mask, fast 0.5 sec.
tic
smg=false(nwx,nwx);%water mask from a priori coastline shapefiles
for j=1:cnt
%   [sx,sy]=polarstereo_fwd(S(j).Y,S(j).X,[], [],70,-45);
    latj=S(j).Y; lonj=S(j).X;
    [sx,sy]=latlon2xy(latj,lonj,projgdal);
    dx=resrc;
    idx=round((sx-rang0(1))/dx)+1;idy=round((sy-rang0(4))/(-dx))+1;
%     idt=isnan(idx)|isnan(idy);
%     idx(idt)=[];idy(idt)=[];
    %Oct 10, 2018:fix bug 14; separate polygons using the NaN;
    M=[idx(:),idy(:)];
    idx = any(isnan(M),2);
    idy = 1+cumsum(idx);
    idz = 1:size(M,1);
    C = accumarray(idy(~idx),idz(~idx),[],@(r){M(r,:)});
    for k=1:length(C)
    idx=C{k}(:,1);idy=C{k}(:,2);   
    sm=poly2mask(idx,idy,nwx,nwx); % fast, apply to each polygon one by one.
    smg=smg|sm;
    end
end
wm=[];wm.x=xw;wm.y=yw;wm.z=smg;clear smg; %1 land, 0 water
if(sum(wm.z(:))==0||sum(wm.z(:))==nwx*nwx);fprintf('This tile contain no coastline (all land or all ocean).');return;end % if all land or all water, i.e. no coastline.
% get the buffer region of coastline: +- widm of the coastline;
Md1 = imdilate(wm.z, ones(widm/resrc));
Me1=imerode(wm.z, ones(widm/resrc)); % erode the box edge also; to fix: change widm0 in CoastTileMain to be larger than widm;
Mcb=Md1&~Me1; %coastal band
%to save space, void the buffer of the box
% figure;imagesc(wm.x,wm.y,Mcb);colorbar %plot the coastal band
toc

%create file 'imagelist.dat' which include lists of images only within the selected tile.
% takes 2 hours for 149129 images.
[range,XYbg,~,f,fdir]=creatsublist(rang0,range,XYbg,projgdalg,f,fdir); 


x=[range(:,1) range(:,2) range(:,2) range(:,1) range(:,1) ];y=[range(:,4) range(:,4) range(:,3) range(:,3) range(:,4) ];
id=find(range(:,2)>rang0(1) & range(:,1)<rang0(2) & range(:,4)>rang0(3) & range(:,3)<rang0(4));
% id=1:length(range(:,1));

fprintf ('\n Step 1: Preparing the grouping of common overlapping pieces.')
fprintf ('\n Step 1.1: getting the real Polygon boundary for all files over the output zone.')
%load all reg.txt and meta.txt files in the coverage
idregion=id;XYb=cell(size(idregion));dzxy=cell(size(idregion));count=0;
rmsreg=zeros(size(idregion));idd=[];dzxyd=zeros(length(idregion),3);
flagcb=zeros(size(idregion));
for j=1:length(idregion)
    i=idregion(j);
 % get the Polygon boundary for actual data
        demdir=fdir{i};
        infile= [demdir,'/',f{i}];
        %[XYbi,rangei]=imagebd(infile);
	XYbi=XYbg{i};
	Xb=XYbi(:,1);Yb=XYbi(:,2);
        XYb{j}=XYbi;
        
        % check whether this polygon intersect with the coastline band.
        idx=round((Xb-wm.x(1))/resrc)+1;
        idy=round((Yb-wm.y(1))/(-resrc))+1;
        Mb = poly2mask(idx,idy, nwx,nwx); % build polygon mask       
        overl=Mb&Mcb;
        if(sum(sum(overl))~=0);flagcb(j)=1;end
end
% only keep the strips that cover coastline.
idd=flagcb==0;
idregion(idd)=[];XYb(idd)=[];dzxy(idd)=[];dzxyd(idd,:)=[];rmsreg(idd)=[];% %hi
if(isempty(idregion));fprintf('No images along coastline.');return;end % 
fprintf(['\n \n',num2str(length(idregion)),' images in this region before coregistration. \n']);
%End of loading


if flagplot==1
%Plot the coverage
figure
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
hold all
% plot(x(id), y(id) ,'.-')
for i=1:length(id)
    plot(x(id(i),:), y(id(i),:),'b>-','linewidm',4)
%     hold off
end
hold on
plot(xeq,yeq,'r*','Markersize',12)
hold on;plot(rang0([1 2 2 1 1]), rang0([3 3 4 4 3]), 'r-')

%[lat,lon]=polarstereo_inv(x,y,[], [],70,-45);
[lat,lon]=xy2latlon(x,y,projgdal);
% lon(lon>=0)=lon(lon>=0)-360;
figure
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
hold all
% plot(lon', lat' ,'-')
for i=1:length(id)
    plot(lon(id(i),:), lat(id(i),:),'b>-','linewidm',4)
%     hold off
end
plot(loneq, lateq ,'r*','Markersize',12)
end %if flagplot

fprintf ('\n Step 1.2: Searching for overlappings of Polygons at regular big grids.')
fprintf (['\n In case of too many repeats, use only the summer and latest (in time) novmax*3=',num2str(novmax*3),' measurements. \n'])
tic
%sub-zones for searching overlappings
dx=800;%400;%3600;%7e3;%%17e3; %3600; %80; %3600; 
dxov=dx/4;%200;
rang0t=round(rang0/dxov)*dxov;
rangx=rang0t(1):dx:rang0t(2);nx=length(rangx)-1;
rangy=rang0t(3):dx:rang0t(4);ny=length(rangy)-1;
ns=nx*ny;rang0s=zeros(ns,4);
novlp=zeros(ny,nx);baovlp=zeros(ny,nx);idg=cell(ns,1);%idg{ns}=[];
enl=0;%;0.3;
for ix=1:nx
    for iy=1:ny
        ixy=iy+(ix-1)*ny;
        rang0s(ixy,:)=[rangx(ix) rangx(ix+1) rangy(iy) rangy(iy+1) ];
        id=find(range(:,1)<=rangx(ix)-enl*dx & range(:,2)>=rangx(ix+1)+enl*dx & range(:,3)<=rangy(iy)-enl*dx  & range(:,4)>=rangy(iy+1)+enl*dx);
%         novlp(iy,ix)=length(id);
        % Refinement of the overlapping count, checking the actual data
        % coverage overlapping with the subzone.
        x0si=[rang0s(ixy,1)-enl*dx rang0s(ixy,2)+enl*dx rang0s(ixy,2)+enl*dx rang0s(ixy,1)-enl*dx rang0s(ixy,1)-enl*dx ];
        y0si=[rang0s(ixy,4)+enl*dx rang0s(ixy,4)+enl*dx rang0s(ixy,3)-enl*dx rang0s(ixy,3)-enl*dx rang0s(ixy,4)+enl*dx ];
        idd=[]; str=cell(length(id),1);
        for j=1:length(id)
            % get the Polygon boundary for actual data
            i=id(j); 
	    str{j}=f{id(j)}(1:13);

            M=idregion==i;nt=sum(M); % only work on the strips that are selected (cover coastline (and successfully coregistered)).
            if nt==0 ;idd=[idd;j];continue;end
            Xb=XYbg{i}(:,1);
            Yb=XYbg{i}(:,2);

            % new method
            in = inpolygon(x0si,y0si,Xb,Yb); %whether subtiles are inside the polygon
            if any(in==0) %any subtile corners not in the boundary
               idd=[idd;j];
            end
        end % j=1:length(id)
        id(idd)=[];str(idd)=[];

        %get rid of the strips have the same date and same sensor
        [un idx_last idx] = unique(str(:));
        id1=1:length(id);idd=id1(~ismember(id1,idx_last));
	id(idd)=[];str(idd)=[];

	%Oct 9, 2018
	%in case of too many repeats, use only the summer and latest (in time) novmax*3 measurements for grouping.
	[idd,idkp]=selectdate(str,novmax*3);
	id(idd)=[];str(idd)=[];

        novlp(iy,ix)=length(id);
        idg{ixy}=sort(id);   %finding the DEMs at each zones.    
    end
end
display(['Counting overlapping...']);
toc

% [X,Y]=meshgrid(rangx(1:end-1),rangy(1:end-1));
[X,Y]=meshgrid((rangx(1:end-1)+rangx(2:end))/2,(rangy(1:end-1)+rangy(2:end))/2);
if flagplot==1
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
surf(X*1e-3,Y*1e-3,novlp);colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
hold on
plot3(xeq*1e-3, yeq*1e-3,1e2 ,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
saveas(gcf,'OverlappingCount_Alaska','fig')

%[LAT,LON]=polarstereo_inv(X,Y,[], [],70,-45);
[LAT,LON]=xy2latlon(X,Y,projgdal);
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
surf(LON,LAT,novlp);colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
hold on
plot3(loneq, lateq,1e2 ,'r*','Markersize',12)

saveas(gcf,'OverlappingCount','fig')
out=[LAT(:) LON(:),];
save -ascii grid.dat out
out=[LAT(:) LON(:) novlp(:)];
save -ascii latlonnov_barnesstrip8m.dat out
end

fprintf ('\n Step 1.3: grouping grids, grid in each piece has the same overlapping files.')
%identify the overlapping pieces based on the same files id
%https://www.mathworks.com/matlabcentral/newsreader/view_thread/304505
idgs=cell(ns,1); % ns: total number of big grids
for i=1:length(idg);idgs{i}=num2str(idg{i}(:)');end
[un idx_last idx] = unique(idgs(:));
unique_idx = accumarray(idx(:),(1:length(idx))',[],@(x) {(x)});
%Exclude small pieces;
suma=zeros(length(un),1);
for j=1:length(un)
suma(j)=length(unique_idx{j})*dx*dx; % total area of each piece.
end
%almt=1e3*1e3;% 2e3*2e3;%1e3*1e3; %Adjust this parameter; minimum areas for each piece %hi
idd=find(suma<max(dx*dx,almt));idnn=1:length(un);
idkp=idnn(~ismember(idnn,idd));idxn=zeros(size(idx));
for i=1:length(idx)
    idn=find(idkp==idx(i));
    if isempty(idn) idxn(i)=0;else;idxn(i)=idn;end
end
un(idd)=[];idx=idxn;
npc=length(un); %total number of groups
%for plotting
idun=zeros(ny,nx); %id for each piece (un)
for ixy=1:ns
    ix=floor((ixy-1)/nx)+1;
    iy=ixy-(ix-1)*nx;
    idun(iy,ix)=idx(ixy);
end

% For those deleted pieces, check if they form large connectted areas; 
% if so, restore some.
[un,idun]=restoregrid(un,idun,dx,idg,novlp);
npc=length(un); %total number of groups

if flagplot==1
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
surf(X*1e-3,Y*1e-3,idun);colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
hold on
plot3(xeq*1e-3, yeq*1e-3,1e2 ,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
% hold on;plot3(xriv(:,1)*1e-3,yriv(:)*1e-3,1e3*ones(size(riv(:,1))),'b-','linewidm',4)
saveas(gcf,'OverlappingPieces','fig')
end

% Initialize
nsel=npc;nsuby=length(yout);nsubx=length(xout);
% jump=zeros(nsuby,nsubx);%jumpstd=zeros(nsuby,nsubx);
%latout=zeros(nsuby,nsubx,3); %lat lon demref
%if flagplot==1 % need iselop in estimating tidal heights.
iselop=zeros(nsuby,nsubx,'int32'); %isel id of the chosen group
%end
%jumpc=cell(nsuby,nsubx); %collecting jump information
% lenfc=cell(nsuby,nsubx); %collecting jump information
% jumpcm=zeros(3,10,nsuby,nsubx);%collect in matrix for faster speed.
%Gather the final mask
jump=-1*ones(nsuby,nsubx,'int8'); %%-1,non value;1 water; 0 non-water
% nrpt=zeros(nsuby,nsubx);
novlpf=zeros(nsuby,nsubx,'int32');%maximum repeat of the chosen group at a pixel
novlpfc=zeros(nsuby,nsubx,'int32');%maximum repeat of the chosen group at a pixel before excluding clouds.
prob=255*ones(nsuby,nsubx,'uint8'); %probability of the final water mask at each pixel. 255 no data
%oflagc=255*ones(nsuby,nsubx,'int32');
Medgsib=false(nsuby,nsubx); %collect of all images edge lines.
%fid2 = fopen('trend.txt', 'w');
%fid4 = fopen('demTyndalbig.dat','w');

% Preallocate for genearting water mask files for each strip.
idg1=[];
for i=1:length(idg);idg1=[idg1;idg{i}(:)];end
isv=unique(idg1);
if isempty(isv);return;end
datarsv(length(isv))=struct('x',[],'y',[],'z',[]);
%load ../sv4.mat
% datamtrsv(length(isv))=struct('x',[],'y',[],'z',[]);
fprintf(['\n\n ',num2str(length(isv)),' images in this region after Step 1.3. \n']);

%Estimating tidal height for each image.
save test1tide.mat f range isv
[tideh]=gettide(f,range,isv);
save test1tide.mat f range isv tideh

%Coregistration after data reduction in Step 1.2.
fprintf ('\n Step 1.X: Coregister all orthoimages to each other. \n')
tic
idregion=isv; %index of f cells. isv is the selected index after the data reduction.
pg=zeros(length(idregion),3);dX4Sg=pg;
idd=[];
%Refer to river/rivergithub2/rivermainserial/ChangeRiver.m
%hi: worry that the orthorectification with a geoid may cause distortion to image on land, hence bad translational offsets.
if coregflag ~=0  
%tag.z 1 land area (used for coregistration), 0 ocean void area
tag=wm;tag.z=wm.z&Mcb; %within coastal band; 1 land area; 0 void area;
if ~exist('dX4Sg.mat','file')
  [dX4Sg,idregion,data0r]=coreg3image(rang0,idregion,XYbg,f,fdir,tag); %only load data within rang0;
 %[dX4Sg]=pxpymono(dX4Sg2,f2is,fdir2is,XYb2is,fis,fdiris,XYbis); %12 hours for 86 images;14 hours for 86 images
  save dX4Sg.mat dX4Sg idregion XYbg f fdir -v7.3
else
  load dX4Sg.mat %hi
  fprintf(['\n Loading existing dX4Sg.mat.'])
  [n1,m1]=size(dX4Sg);
  if n1~=length(idregion)
    warning('CoastTileMono.m: loaded dX4Sg not matching the ids!')
    return
  end
end
fprintf ('\n Done with coregistration.\n ')
toc
fprintf(['\n ',num2str(length(idregion)),' images in this region after coregistration.']);
clear data datar data0r
end %if coregflag


fprintf ('\n Step 1.4: loading of data, calculating NDWI, and assigning chosen output images for pixels in each group.')
% DEM collection and time series for each zone
%Find the pair closest to epicenter
% [iy,ix]=find((((X-xeq).^2+(Y-yeq).^2))<5e3.);
%Find the maximum overlapping pairs closest to epicenter
[iys,ixs]=find( novlp>=0);
% [iys,ixs]=find( baovlp==1);
clear Xt Yt; for i=1:length(ixs);Xt(i)=X(iys(i),ixs(i));Yt(i)=Y(iys(i),ixs(i));end
[dist,im]=min((Xt-xeq).^2+(Yt-yeq).^2);

display(['nsel=',num2str(nsel)])
tcpu1 = cputime;
ck0=clock;
output=[];
%tag=readGeotiff([macdir,'/home/chunli/scripts/Barnesrock.tif']);
%load([macdir,'/home/chunli/scripts/Barnesrock.mat']);
% load([macdir,'/home/chunli/scripts/BarnesrockRGI.mat']); %tag
ck1=clock;
%parpool(4)
%parfor isel=1:60
%for isel=1:npc  %nsel 
%display(['isel=',num2str(isel)]);
%id=strread(un{isel},'%d');

%probthre=[50., 5, 95];
probthrein=probthre; clear probthre
probthre=50; %probthrein(1);

save mat1.mat -v7.3

if exist('sv5.mat','file')
   load sv5.mat
   fprintf(['\n Loading existing sv5.mat.'])
else 

for isel=1:npc %[139,177];%[163,176];%1:npc %[97,108]%115 %1:npc  %[127, 134]
display(['isel=',num2str(isel)]);
id=strread(un{isel},'%d');%id of overlapping files


%August 2023; Oct 9, 2018
%in case of too many repeats, use only the novmax measurements
if length(id) >= novmax;
fprintf (['\n Too many repeats, use only the novmax = ',num2str(novmax),' measurements. \n'])
str=f(id)';
[idd,idkp]=selectdate(str,novmax);
id=id(idkp);
end

for j=1:length(id)%0%length(id)
% display(['Overlapping Files Date: ',f{id(j)}(6:13)])
display(['Overlapping Files Date: ',f{id(j)}])
end
if length(id)<1;continue;end

rang0sov=[max(range(id,1)) min(range(id,2)) max(range(id,3)) min(range(id,4))];
tx=rang0sov(1):resr:rang0sov(2);ty=rang0sov(4):-resr:rang0sov(3); %add the constraint of a priori coastline
[X,Y]=meshgrid(tx,ty);
M1=interp2(wm.x,wm.y,Mcb,tx,ty','*nearest',0); %Mcb 1  coastal band, 0 non coast
if sum(M1(:))==0; continue;end
rang0sov=[min(X(M1)) max(X(M1)) min(Y(M1)) max(Y(M1))];

ranget=[rang0sov/resr];
rang0sov=[ceil(ranget(1)) floor(ranget(2)) ceil(ranget(3)) floor(ranget(4))]*resr;
rang0st=[max(rang0sov(1),rang1(1)) min(rang0sov(2),rang1(2)) max(rang0sov(3),rang1(3)) min(rang0sov(4),rang1(4))]; %overlapping grid
x0st=[rang0st(1) rang0st(2) rang0st(2) rang0st(1) rang0st(1) ];y0st=[rang0st(4) rang0st(4) rang0st(3) rang0st(3) rang0st(4) ];

%plot the overlapping zone and two image boundaries
if 0
figure  %(1) 
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
plot(x0, y0,'g-','linewidm',4)
hold all
% plot(x0s', y0s','r-','linewidm',2)
plot(x0st, y0st,'r-','linewidm',6)
plot(x(id,:)', y(id,:)' ,'b>-') 
plot(xeq, yeq ,'r*','Markersize',12)
axis equal
plot(X(idun==isel),Y(idun==isel),'ro') %overlapping grids

hold all
plot3(x0st*1e-3, y0st*1e-3,100*ones(size(x0st)),'r-','linewidm',6)
end

%find the overlapping zone for this piece
mx0=find(xout>=rang0st(1) & xout<=rang0st(2) );
my0=find(yout>=rang0st(3) & yout<=rang0st(4) );
nyj=length(my0);nxj=length(mx0);
% demg=nan*ones(nyj,nxj,length(id));
demg=-1*ones(nyj,nxj,length(id),'int8');
[nyj,nxj]
if nyj<=1||nxj<=1;continue;end

%edge of the overlapping zone
Medgsib(my0(1),mx0(:))=1;Medgsib(my0(end),mx0(:))=1;Medgsib(my0(:),mx0(1))=1;Medgsib(my0(:),mx0(end))=1;
%the outer ring
oidx=max(mx0(1)-1,1):min(mx0(end)+1,nsubx);
oidy=max(my0(1)-1,1):min(my0(end)+1,nsuby);
Medgsib(oidy(1),oidx(:))=1;Medgsib(oidy(end),oidx(:))=1;Medgsib(oidy(:),oidx(1))=1;Medgsib(oidy(:),oidx(end))=1;
    
% get coastline/watermask for each file, and save it.
    t=zeros(length(id),1);
%   [~,idsort]=sort(t);id=id(idsort); %sort the id based on tim
    idd=[];
	for j=1:length(id)
	fprintf(['\n Working on strip ',num2str(j), ': ',f{id(j)}])
	ymd=f{id(j)}(6:13);i=id(j); 
	t(j)=datenum(ymd,'yyyymmdd');
	iisv=find(isv==i);
	if isempty(datarsv(iisv).x) % non exist
        demdir=fdir{i};
        infile= strrep([demdir,'/',f{i}],'.xml','.xml');
        clear data
	% find the overlapping range of data range and coastal band.
	rang0sov=range(i,:);
	tx=rang0sov(1):resr:rang0sov(2);ty=rang0sov(4):-resr:rang0sov(3); %add the constraint of a priori coastline
	[X,Y]=meshgrid(tx,ty);
	M1=interp2(wm.x,wm.y,Mcb,tx,ty','*nearest',0); %Mcb 1  coastal band, 0 non coast
	rangeov=[min(X(M1)) max(X(M1)) min(Y(M1)) max(Y(M1))];
	clear X Y
        tic
        data=multispecmono(infile,wm,rangeov); %given strip meta file, finding all image 1 multispectral imageries, orthrectiying, get water mask.
        display(['Getting water mask of this strip ',num2str(j)]);
        toc
	[n1,m1]=size(data.z);
        if isempty(data.z) || min([n1,m1])<=2 
        %Modification Sept. 25, 2018: save the empty results to avoid repeat calculation of NDWI.
            datar= struct();
            datar.x=1;datar.y=[];  datar.z=[]; 
            datarsv(iisv)=datar;

            idd=[idd;j];

	    %June 20, 2024;
	    if min([n1,m1])<=2; fprintf(['\n Mask size is <=2:',infile,'.\n']);end
            continue
        end  
    

        %coregister to a reference (e.g. ICESat). reg.txt file or feature tracking.
%         dzxyt=dzxyd(idregion==i,:);%dzxy{idregion==i};%lots of them missing along coast.
        idt1=find(idregion==i);
        if isempty(idt1)
            dzxyt=[0 0 0];
            fprintf(['Image coregistration not found (failed), remove this image: ', infile])
            
            %Modification Sept. 25, 2018: save the empty results to avoid repeat calculation of NDWI.
            datar= struct();
            datar.x=1;datar.y=[];  datar.z=[]; 
            datarsv(iisv)=datar;
            idd=[idd;j];
            continue
            
        else
            dzxyt=dX4Sg(idregion==i,:); 
            dX4S=dzxyt(2:3); 
        end

fprintf(['\n Loading its translational offsets pxy (',num2str(j),'/',num2str(length(id)),', isel ',num2str(isel),'): ',num2str(dX4S(:)')])
        
        %make coordinates exactly on (resr*n) xout yout grids.
        data.x=data.x-dX4S(1);data.y=data.y-dX4S(2);
        if ~isempty(data.coast)
        data.coast(:,1)=data.coast(:,1)-dX4S(1); data.coast(:,2)=data.coast(:,2)-dX4S(2);
        end
        
        ranget=[[min(data.x) max(data.x) min(data.y) max(data.y)]/resr];
        ranget=[ceil(ranget(1)) floor(ranget(2)) ceil(ranget(3)) floor(ranget(4))]*resr;
        tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);
	%July 2024
	n1=length(tx);m1=length(ty);
        if  min([n1,m1])<=2
            datar= struct();
            datar.x=1;datar.y=[];  datar.z=[];
            datarsv(iisv)=datar;
            idd=[idd;j];
            continue
        end

	try
        tz = interp2(data.x,data.y,data.z,tx,ty','*nearest',-1);
%       tz(isnan(tz))=-1;tz=int8(tz); %DOUBLE TO INT8
        datar= struct();
        datar.x=tx;datar.y=ty;  datar.z=tz; %datar.coast=data.coast;
        datarsv(iisv)=datar; %possible crash with parallel
	catch e %July, 2024
                fprintf('\n CoastTileMono.m There was an error! The message was:\n%s',e.message);
	    datar= struct();
            datar.x=1;datar.y=[];  datar.z=[];
            datarsv(iisv)=datar;
            idd=[idd;j];
            continue

	end

	%collect all image edge lines
	tzedg=datar.z==-1;
	tzd1=imdilate(tzedg, ones(3));
	tze1=imerode(tzedg,ones(3));
	tzl=tzd1&(~tze1); %edge lines
	%[X,Y]=meshgrid(datar.x,datar.y);
	try %
	tz = interp2(datar.x,datar.y,tzl,xout,yout','*nearest',0); %check time
	Medgsib=Medgsib|tz; %edge 1; non edge 0
	catch e %July, 2024
                fprintf('\n CoastTileMono.m There was an error! The message was:\n%s',e.message);
	end
      
	else %load water mask data for collecting demg
        datar=datarsv(iisv);

	%Modification Sept. 25, 2018:skip empty results.
	[n1,m1]=size(datar.z);
	if isempty(datar.z) || min([n1,m1])<=2
	    datar= struct();
            datar.x=1;datar.y=[];  datar.z=[];
            datarsv(iisv)=datar;
            idd=[idd;j];
            continue
        end

	end %if
    
    %assign DEM to the overlapping zone
    idx=find(datar.x>=rang0st(1) & datar.x<=rang0st(2));
    idy=find(datar.y>=rang0st(3) & datar.y<=rang0st(4));
    if (length(idy)~=nyj || length(idx)~=nxj);warning(['Size mismatch, interpolate instead, isel:',num2str(isel)]);
         tx=xout(mx0);ty=yout(my0);
	 try
         tz = interp2(datar.x,datar.y,datar.z,tx,ty','*nearest',-1);
         demg(1:nyj,1:nxj,j)=tz;
	 catch e %July, 2024
                fprintf('\n CoastTileMono.m There was an error! The message was:\n%s',e.message);
         end
        % idd=[idd;j];
    else
         demg(1:nyj,1:nxj,j)=datar.z(idy,idx);
    end
    
	end %for j
% % end of loading
	id(idd)=[];demg(:,:,idd)=[];t(idd)=[];
%save work1.mat demg
        demgc=demg;% total demg before excluding clouds.

        if 0 %method 1; if one image is mistakenly all water(e.g. GEOEYE image), all other images will be mistakenly deleted in clouds detection.
        % pick the larger values at the overlapping area. ...
        M1=max(demg,[],3); %change M1 from double to int8;
    % 	M1(isnan(M1))=-1;M1=int8(M1);

        else %method 2 ; Use the water probability, tend to reserve more data in clouds detection. Difference(Methods 1/2) not significant.
         Msum=sum((demg==1),3); %demg: 1, water, 0 non water, -1 void data.
        Msuma=sum((demg~=-1),3); %repeats of non void data.
        Msr=uint8(Msum./Msuma*100);Msr(Msuma==0)=255;
	M1=-1*ones(nyj,nxj,'int8');
        M1(Msr>=probthre&Msr~=255)=1;
        M1(Msr<probthre)=0;
        M1(Msuma==0)=-1;
        end %if 0

        %based on water coverage to exclude clouds on water. 
        % It's useful to let you know when no good data appears all.
        cloudth=cloudarea/resr/resr;%1e3*500/resr/resr % %1km by 500m clouds over water area 
        if cloudflag % cloud detection or not 
        nlb=widm2/resr; % widm2: long thin (100 m widm) beach band
        idd=[];
        for j=1:length(id)
            demgt=demg(:,:,j); %-1,non value;1 water; 0 non-water
            if length(id)==1 && 0 % use a priori information as referece; good on ocean, bad on coast due to poor resolution
                M1t=interp2(wm.x,wm.y,Md1,xout(mx0),yout(my0)','*nearest',1);% Md1:1 land ; 0 water
                M1t=int8(~M1t);%1 water
            else % use the stacking DEM as reference
                M1t=M1;
            end %
            M1t(demgt==-1)=-1;
            Mextra=abs(M1t-demgt); %extra non water area, which means clouds. 1 if clouds, or land (vs lake).
		%1 extra land, -1 extra water/could be shadows;trust the stacked mask
            
            Mextra3=icluster(Mextra,nlb,demgt,cloudth); %all clusters that are clouds over water
            
            cloudsum=sum(sum(Mextra3)); %pixels
            if cloudsum > cloudth 
                idd=[idd;j];
            end
        end       
        id(idd)=[]; demg(:,:,idd)=[];t(idd)=[];
        fprintf('\n Clouds detection. \n')
        end
        
        % Rerun the previous step to exclude the clouds that might have wrong translation parameters.
        % Recommend the following line, since the cloud shadow might introduce false water areas.
        M1=max(demg,[],3); % if not using this line, we have more coverage with more repeats but poor quality. 
        Msum=sum((demg==1),3); %1, water, 0 non water, -1 void data.
        Msuma=sum((demg~=-1),3); %repeats of non void data.
        Msr=uint8(Msum./Msuma*100);Msr(Msuma==0)=255;%Msr(~(Msr>=0&Msr<=100))=255;
        if isempty(M1) ;continue;end
        
        %store the final water mask for the tile.
        %matrix fast, 0.6sec vs 128 sec
        lenfg=sum((demg~=-1),3); % clouds excluded in demg.
        lenfgc=sum((demgc~=-1),3);%repeats include clouds
        lenfpre=novlpf(my0,mx0);novlpfsc=novlpfc(my0,mx0);
%       M=lenfg>lenfpre; 
        M=lenfg>lenfpre|(lenfg==lenfpre&lenfgc>novlpfsc); %hi ; choose the more data groups.
        %collect count of repeats of all data including clouds.
        novlpfsc(M)=lenfgc(M);novlpfc(my0,mx0)=novlpfsc;
        %collect count of repeats without clouds.
        novlpfs=novlpf(my0,mx0);
        novlpfs(M)=lenfg(M);novlpf(my0,mx0)=novlpfs;
        %collect water mask M1
        jumps=jump(my0,mx0);
        jumps(M)=M1(M);jump(my0,mx0)=jumps;
	%if flagplot==1 % need iselop in estimating tidal heights.
        %store group isel
        iselops=iselop(my0,mx0);
        iselops(M)=int32(isel);iselop(my0,mx0)=iselops;
	%end
        % %probability for this group isel
        probi=255*ones(size(M1),'uint8');  
        %0.1, only one image; 0.2, only one image after excluding clouds;
        if 1 % assgin for whole scene.
        probi(lenfgc==1)=uint8(0.1*100); 
        probi(lenfgc>1&lenfg==1&Msr~=0&Msr~=255)=uint8(0.2*100); %more reliable after cloud detection
	probi(lenfgc>1&lenfg==1&Msr==0)=uint8(0); %to fix bug 11;
        else  % assgin for water area only
        probi(lenfgc==1&Msr~=0&Msr~=255)=uint8(0.1*100); 
        probi(lenfgc>1&lenfg==1&Msr~=0&Msr~=255)=uint8(0.2*100);
        end
        probi(lenfg>1)=Msr(lenfg>1);
        probs=prob(my0,mx0); %select
        probs(M)=probi(M);prob(my0,mx0)=probs;
	%oflagc 0 no data; 1 one data; 2 one data after excluding clouds; 10 >= two good data; -> just use novlpfs

        display(['Collecting the final water mask for piece isel= ',num2str(isel),'...']);
%save work2.mat jumpc
end % for isel
%save sv4.mat datarsv isv -v7.3
save sv5.mat xout yout iselop novlpf novlpfc prob jump Medgsib datarsv isv -v7.3
if flagplot==1 
%save sv5.mat xout yout iselop novlpf novlpfc prob jump Medgsib -v7.3
else
%save sv5.mat xout yout novlpf novlpfc prob jump Medgsib -v7.3
end
whos
ck2=clock;
dt=ck2-ck1;tsec=dt(4)*3600+dt(5)*60+dt(6);
display(['Loading files and NDWI calculation take: ',num2str(tsec),' sec'])
end % if sv5.mat does not exist

%clear datarsv data 

%%%-1,non value;1 water; 0 non-water
% jump(novlpfc<=1)=-1; %hi rather missing data than wrong coastlines. %move after the fill based on a priori coastline to avoid wrong coastline along eroded a priori coast.
% Use probability, if probs<0.1, consider change the area to land. jump(probs<0.1)=0;

%write
probo=prob;
%probo(prob==0)=uint8(255); %make land transparent
%projstr='polar stereo north';
projstr=projstrin;
%OutName='prob.tif';
OutName=ofile2;
writeGeotiff(OutName,xout,yout,probo,1,255,projstr)
writeGeotiff(ofile3,xout,yout,uint8(novlpf),1,255,projstr)
clear probo novlpfc

%probthre=[50., 5, 95];
for probthre=probthrein

fprintf (['\n Coastline based on water probability threshold: ',num2str(probthre),' %. '])

%Based on probability only, prob=[0, 100]; if 255, it means no data.
jump=-1*ones(nsuby,nsubx,'int8'); %%-1,non value;1 water; 0 non-water
jump(prob>=probthre&prob~=255)=1;
jump(prob<probthre)=0;
Me1o=interp2(wm.x,wm.y,Me1,xout,yout','*nearest',1);%1 land, 0 water
jump(Me1o)=0;   %hi; fill inland with land;%for inland regions far from coast, set it as edges.
Md1o=interp2(wm.x,wm.y,Md1,xout,yout','*nearest',1);
jump(Md1o==0)=1;   %hi; fill the ocean with water
Medgsap=Me1o|~Md1o; %set the a priori buffer zone as edges.
clear Me1o Md1o
%jump(prob==255|novlpf<=1)=-1; %considear novlpf<=2 %no data/edge; Move after fill (Modfil) to mitigate the separation of water bodies (Bug 6).

%get coast from water mask.
%apply the coastline buffer band to the mask
%Medgs=isnan(jump);
Medgs=(jump==-1);%-1,non value;1 water; 0 non-water
Modj=jump;Modj(Medgs)=0;%edges as land,but would be filterd out again later.
Modj= bwareaopen(Modj, lakearea/resr/resr); %remove small clusters; smallest lake that would be kept;
%Around edge of image, add a thin band of water, so that the small area at edge be filled; then later remove it; 
Med=(imdilate(~Medgs,ones(4)))&Medgs;
Modj(Med)=1;%1 is water
%smallest island is 46 meters long by 16 meters (Bishop Rock).
Modfil = bwareaopen(~Modj, cloudarea/resr/resr); %fill small areas: 1e4*4m^2, smallest island that would be kept.
Modfil=~Modfil;
Modfil(Med)=0;%remove the artifical water around edge;
clear Med Modj %jump

%get the boundaries; 
%remove the edges of the Tile Box.
Medgstb=zeros(size(Modfil)); 
Medgstb(:,1)=1;Medgstb(:,end)=1;Medgstb(1,:)=1;Medgstb(end,:)=1;

%remove the buffer widm (widm0, 3km) of the tile box;
idx=xout<rang0b(1)|xout>rang0b(2);idy=yout<rang0b(3)|yout>rang0b(4);
Medgstb(:,idx)=1;Medgstb(idy,:)=1;

Medgsib(novlpf>=20)=0;%do not apply edges (especially image boundaries) if number of repeats > 20.
Medgs=Medgs|(prob==255|novlpf<=novlmt)|Medgsap; %add more edge based on repeats.
Medgs1=imdilate(Medgs,ones(3))|Medgstb|Medgsib; % combine all edges. Medgsib Medgsap

%save image boundary file
writeGeotiff(ofile4,xout,yout,uint8(Medgs1),1,255,projstr)

%B = bwboundaries(~Modfil,'noholes'); %^_^ keep the unwanted lake islands.
B = bwboundaries(~Modfil); % Delta area separated from ocean by road is treated as a like. Keep all large lakes.
n=length(B);xo=cell(n,1);yo=cell(n,1);
%clear prob Medgsib Medgstb
if flagplot~=1
%clear novlpf
end
tic
peri=sqrt(4*pi*cloudarea/resr/resr); % shortest perimeter of a given area.
for k=1:n %can be slow if poor data quanlity, lots of scatterred points;workstereoprob2wocoregcrop/55_06_2_1_coast_v1.0.shp includes island of lakes!
    xid=B{k}(:,2); yid=B{k}(:,1);zid=ones(size(xid));
    if (length(xid)<peri);continue;end %to be immune to bug 11.
    %get the boundaries in terms of mask
    BWb=zeros(size(Modfil));BWb((xid-1)*nsuby+yid)=1;
    Mt=BWb&Medgs1;%filter out edges
    ne=sum(Mt(:));
    if ne~=0 %find the id and replace them with NaN;
       [idy,idx]=find(Mt==1);
       
       for j=1:ne
          Ml=xid==idx(j)&yid==idy(j);
          zid(Ml)=0;%fall into edges
       end
       
    end
    x=xout(xid);y=yout(yid);    
    %[LAT,LON]=polarstereo_inv(x,y,[], [],70,-45);
    [LAT,LON]=xy2latlon(x,y,projgdal);
%     LAT(zid==0)=nan;LON(zid==0)=nan;%polygon but with only valid polylines displayed.
    zidi=find(zid==0);zidid=zidi(2:end)-zidi(1:end-1);
    %id=find(zidid==1);%delete sequential nans
%    idx=zidi(id+1);
    idb=find(zidid~=1);M=zidi;M(idb+1)=nan;
%   idxb=zidi(idb+1); LAT(idxb)=nan;LON(idxb)=nan;
    [idx,idn]=separatelines(M,10); %fix Bug 15
    %figure;plot(x,y,'go');hold on;plot(x(idx),y(idx),'ro-')
    LAT(zidi(idn))=nan;LON(zidi(idn))=nan;LON(idx)=[];LAT(idx)=[];
    
    xo{k}=LON;yo{k}=LAT;
end
fprintf('Retrieve boundary')
toc
%clear Mt BWb Medgs1
idd=find(cellfun(@isempty,xo)); %fix bug 12
xo(idd)=[];yo(idd)=[];
%save xoyo.mat xo yo -v7.3

if flagtide==1  %1 carry out the tide calculation; 0 no calculation
    %do nothing; coastlines with be saved with tidal heights. 
else
% shp1 = struct('Geometry', 'PolyGon', 'X', xo, 'Y', yo);
shp1 = struct('Geometry', 'PolyLine', 'X', xo, 'Y', yo);
%shapewrite(shp1, 'coastline.shp');
if ~isempty(shp1)
ofile1=['output/',ifilestr,'coast_thre',num2str(round(probthre)),'_v1.0.shp'];
shapewrite(shp1, ofile1);
end
end

if flagtide==1 % 1 carry out the tide calculation; 0 no calculation

%flagtidemethod=1; %1 accurate but slow; 2 fast but approximate

if flagtidemethod ==1  % accurate but slow;
%CoastTileMonobp2findimage.m
%Check the tidal height along the coastline for the selected image that matchs the final water mask.
%chose to estimate the tide at the final coast. time consuming.
tic
n=length(xo);zo=cell(n,1);zo_time=cell(n,1);
if flagplotany==1
figure;pointsize = 10;hold all;
end
for k=1:n 
    LON=xo{k};LAT=yo{k};Mv=isnan(LON)|isnan(LAT);zo1=nan(size(Mv)); zo_time1=nan(size(Mv));
    LON(Mv)=[];LAT(Mv)=[];
    %[x,y]=polarstereo_fwd(LAT,LON,[], [],70,-45);
    [x,y]=latlon2xy(LAT,LON,projgdal);
    xid=round((x-xout(1))/resr)+1;
    yid=round(-(y-yout(1))/resr)+1;

    iselopc=iselop((xid-1)*nsuby+yid); %iselop on coastlines
    z=zeros(size(LON));z_time=zeros(size(LON));%double
    n1=length(xid);

% point by point % time consuming: 4 days for tile 55_15_1_2_prob_v1.0.tif;
%consider working on every 3 points (6 m) to save time; 100 points
    if 0 % serial , step 1
    for ik=1:n1
	fprintf (['\n flagtidemethod 1, seg ',num2str(k),' , ',num2str(ik),' out of ',num2str(n1),'.'])
        [z_tideh]=matchmasksub(ik,iselopc,xid,yid,xout,yout,jump,un,isv,datarsv,tideh);
	z(ik)=z_tideh;
    end %ik
    
    else % parallel
    ikg=1:100:n1; 
    z_tideh_g=zeros(length(ikg),1); z_time_g=zeros(length(ikg),1);
%   poolobj=parpool(24);
    for j=1:length(ikg)    % point by point % time consuming: 4 days for tile 55_15_1_2_prob_v1.0.tif;
			%consider working on every 3 points (6 m) to save time; 100 points
	ik=ikg(j);
	fprintf (['\n flagtidemethod 1, seg ',num2str(k),' , ',num2str(ik),' out of ',num2str(n1),'.'])
	try
        [z_tideh,iisv]=matchmasksub(ik,iselopc,xid,yid,xout,yout,jump,un,isv,datarsv,tideh);
	%z(ik)=z_tideh;
	z_tideh_g(j)=z_tideh;
    i=isv(iisv); 
    ymdhhmmss=f{i}(6:19);
    z_time_g(j)=str2double(ymdhhmmss);
    	catch e %June 24, 2024
                fprintf('\n CoastTileMono.m There was an error! The message was:\n%s',e.message);
        end

    end %j
%   delete(poolobj)
	
    %assign values to z is run by every 10 pixels.
    for j=1:length(ikg)
        if j==length(ikg)
              id1=ikg(j):n1;
        else
          id1=ikg(j):(ikg(j+1)-1);
        end
        z(id1)=z_tideh_g(j); %use the same tidal height for 10 neighboring pixels.
        z_time(id1)=z_time_g(j);
    end %j
    end %if serial

    zo1(~Mv)=z;zo_time1(~Mv)=z_time;
    zo{k}=zo1;zo_time{k}=zo_time1;
if flagplotany==1
    scatter(LON, LAT, pointsize, z);colorbar
end
%     scatter(xo{k}, yo{k}, pointsize, zo1);colorbar
end %k
if flagplotany==1
colormap jet
saveas(gcf,'coastlinestide1','fig')
end

elseif flagtidemethod ==2 % 2 fast but approximate
%Check the median tidal height at the final coastline based on tides in all used images. 
%Based on assumption: Ideally if the stacking number is large, the tide height of the selected image would be approximately near mean tide.

tic
n=length(xo);zo=cell(n,1);
figure;pointsize = 10;hold all;
for k=1:n 
    BWb=zeros(size(Modfil));
    LON=xo{k};LAT=yo{k};Mv=isnan(LON)|isnan(LAT);zo1=nan(size(LON));
    LON(Mv)=[];LAT(Mv)=[];
    %[x,y]=polarstereo_fwd(LAT,LON,[], [],70,-45);
    [x,y]=latlon2xy(LAT,LON,projgdal);
	%x=xout(xid);
    xid=round((x-xout(1))/resr)+1;
    yid=round(-(y-yout(1))/resr)+1;
    %get the boundaries in terms of mask
    BWb((xid-1)*nsuby+yid)=1;

    iselopc=iselop((xid-1)*nsuby+yid); %iselop on coastlines
    iselopcu=unique(iselopc); %unique iselopc
    z=zeros(size(LON));
    for isel=iselopcu %group by group; or point by point 
	idopc=find(iselopc==isel);

    %get the tides of all used images, and find the median value
    id=strread(un{isel},'%d');%id of overlapping files
    Tideid=zeros(length(id),1);
    for j=1:length(id)
        i=id(j);iisv=find(isv==i);
        Tideid(j)=tideh(iisv);
    end
	z(idopc)=median(Tideid);
    end
    zo1(~Mv)=z;zo{k}=zo1;
    scatter(LON, LAT, pointsize, z);colorbar
%     scatter(xo{k}, yo{k}, pointsize, zo1);colorbar
end
colormap jet
saveas(gcf,'coastlinestide2','fig')

end %if flagtidemethod

fprintf(['Get tides along coastlines method ',num2str(flagtidemethod),' takes: '])
toc

save xoyo.mat xo yo zo zo_time -v7.3

%load xoyo.mat
% 43_51_1_1_coast_tide_thre50_v1.0.shp
ofile1b=['output/',ifilestr,'coast_tide_thre',num2str(round(probthre)),'_v1.0.shp']; 
flagshp=2; %1 points; 2 lines

if flagshp==1 %points
%file size is three times of the coastline file. 
xol=[];yol=[];zol=[];
for k=1:length(xo);xol=[xol;xo{k}(:)];yol=[yol;yo{k}(:)];zol=[zol;zo{k}(:)];end
nl=length(xol);
[shp1(1:nl).Geometry] = deal('Point');
%shp1 = rmfield(shp1, {'X', 'Y', 'Z'});
for k=1:nl;shp1(k).X = xol(k);shp1(k).Y = yol(k);shp1(k).Z = zol(k);;end
%shapewrite(shp1, 't1');

elseif flagshp==2 %lines %best
% Reduce file size: use Polyline, and one attribute (common tidal height) for one segment of line.
%divide line segments based on common tidal heights
clear xo2 yo2 zo2 zo2_time
count=0;
n=length(xo);
for k=1:n 
    
    zo1=zo{k};n1=length(zo1);
    df=diff(zo1);
    idj=find(df~=0); 

    %separate lines by non-zeros.
    M=zeros(n1,1);M(idj+1)=nan;
   [idx,idn]=separatelines(M,0); % M has to be size of n by 1;

    for j=1:length(idn)  %idn the start index of each segment
	count=count+1;
	if j==length(idn)
	  id1=idn(j):n1;
	else
	  id1=idn(j):idn(j+1)-1;
	end
	xo2{count}=xo{k}(id1); yo2{count}=yo{k}(id1); 
	zo1_t1=zo1(id1);zo_time_t1=zo_time{k}(id1);
	if (nanmax(zo1_t1)-nanmin(zo1_t1)) ~=0 
		fprintf(['\n Warning: tides along this segment is not the same for all points; [k j] = ',num2str([k j]),'.\n']);
		%caused by single point merged to previous segment. Not a bit issue.
	end
	zo2{count}=nanmedian(zo1_t1); %tidal height should be the same for all points on this segment.
    zo2_time{count}=nanmedian(zo_time_t1);
    end
end
%make sure zo has one value for each cell.
if count==0
	fprintf('\n xo2 is empty. No output. \n')
	shp1 = struct([]);
else
        shp1 = struct('Geometry', 'PolyLine', 'X', xo2, 'Y', yo2,'TidalHeight',zo2,'SRC_DATE',zo2_time);
end %if count

end %if flagshp

if ~isempty(shp1)
	shapewrite(shp1, ofile1b);
end % isempty

end %if flagtide==1 

% figure;imagesc(xout*1e-3,yout*1e-3,BWb);colorbar;view(0,-90);
% saveas(gcf,'coastlines','fig')

end % for probthre=probthrein

if flagplot==1
[X,Y]=meshgrid(xout,yout);
if 0
figure;imagesc(xout,yout,double(Modfil),'alphadata',~Medgs)
hold on;plot(X(M),Y(M),'ro')
xl=X(M);yl=Y(M);

hold on;plot(xl(probl>0.6)*1e-3,yl(probl>0.6)*1e-3,'r.');
hold on;plot(xl(probl<=0.6&probl>0.2)*1e-3,yl(probl<=0.6&probl>0.2)*1e-3,'m.');
hold on;plot(xl(probl<=0.2&probl>=0.1)*1e-3,yl(probl<=0.2&probl>=0.1)*1e-3,'k.');
hold on;plot(xl(probl<0.1)*1e-3,yl(probl<0.1)*1e-3,'c.');
print('-dpdf','-r300','Tilemask') 
% saveas(gcf,'Tilemask','fig')
savefig(gcf,'Tilemask','compact')

figure;hist(probl)
saveas(gcf,'probabilityw','fig')
end

% plot the chosen piece id for each output pixel
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]);
hold all
xt=xout*1e-3;yt=yout*1e-3;zt=iselop; %to check if the order is consistent with idun
imagesc(xt,yt,zt);
colormap jet;colorbar; shading interp; axis equal; view(0,90)
set(gca,'FontSize', 18);
plot(xeq*1e-3, yeq*1e-3,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
axis(rang0*1e-3)
title('Chosen piece isel')
print('-dpdf','-r300','Chosenisel') 
%saveas(gcf,'Chosenisel','fig')
%savefig(gcf,'Chosenisel','compact')

%Final count of repeats for the chosen piece 
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
imagesc(xt,yt,novlpf);
colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
hold on
plot(xeq*1e-3, yeq*1e-3,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
title('Final repeat, exclude clouds')
print('-dpdf','-r300','OverlappingCount_Final') 
%saveas(gcf,'OverlappingCount_Final','fig')
%savefig(gcf,'OverlappingCount_Final','compact')
close all
end %if flagplot==1

end
