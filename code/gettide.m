function [tideh]=gettide(f,range,isv)
%get tidal heights for each image. 
% isv the id of each images
% load test1tide.mat;

constant

Model='Model_tpxo9';

%get the center of image.
% tideh=zeros(size(isv));
n=length(isv);
SDtime=zeros(n,1);latm=zeros(n,1);lonm=zeros(n,1);
for iisv=1:length(isv)
%        iisv=find(isv==i);
    i=isv(iisv);
    %get time and location of image center.
    ymd=f{i}(6:19);
    SDtime(iisv)=datenum(ymd,'yyyymmddHHMMSS');
    rang0sov=range(i,:);
    xt=(rang0sov(1)+rang0sov(2))/2;yt=(rang0sov(3)+rang0sov(4))/2;
    %[latm(iisv),lonm(iisv)]=polarstereo_inv(xt,yt,[], [],70,-45);
    [latm(iisv),lonm(iisv)]=xy2latlon(xt,yt,projgdal);
end %iisv
M=lonm<0;lonm(M)=lonm(M)+360.;

if flaggage==0 % not use gage data
%Find the closet point of H grids to a target.
% Map the model grid
[lon,lat,H]=tmd_get_bathy(Model); %roughly 18 km grid; resolution 1/6 degree, see http://volkov.oce.orst.edu/tides/tpxo9_atlas.html
H(H<=0)=NaN;
M=~isnan(H);
[LON,LAT]=meshgrid(lon,lat);
X=[LON(M),LAT(M)]; %ocean grids with water column > 0.
XI=[lonm(:),latm(:)]; % image center locations.
tic;[k,d] = dsearchn(X,XI); toc %0.07sec

if flagplotany==1
figure;
imagesc(lon,lat,H,'alphadata',~isnan(H));colormap jet;colorbar;view(0,-90)
for i=1:length(k)
hold on;plot([XI(i,1),X(k(i),1)],[XI(i,2),X(k(i),2)],'g>-')
hold on;plot(XI(i,1),XI(i,2),'ro')
legend('nearest tide grid','image center')
end
saveas(gcf,'tidegrids','fig')
end

latp=X(k,2);lonp=X(k,1);
[tideh,conList]=tmd_tide_pred(Model,SDtime,latp,lonp,'z'); %output: tideh

elseif flaggage==1 %use gage water level, assume all tidal variations within 50 km is the same as the gage
    %given SDtime, lat_gage lon_gage; get tideh
%   gagefile='nomegage/Nome2010to2022.csv';lat_gage=64+29.7/60;lon_gage=-(165+26.4/60)+360; %Nome gage 
%     gagetable=readtable('nomegage/Nome2010to2022.csv'); % undesirably changed the format
%     [nrow,mcol]=size(gagetable);
fid = fopen(gagefile);
fgetl(fid) %reads line but does nothing with it; header     'Date,Time (GMT),Predicted (m),Preliminary (m),Verified (m)'
gagetable = textscan(fid, '%s %6.6d %f %s %f', 'Delimiter',','); % you will need to change the number   of values to match your file %f for numbers and %s for strings.
fclose (fid)

[nrow,~]=size(gagetable{1,1});
    
    for i=1:nrow
        hhmmss= num2str(gagetable{1,2}(i),'%6.6d');
        ymd=[gagetable{1,1}{i},hhmmss];
        gagetime(i)=datenum(ymd,'yyyymmddHHMMSS');
        gagepred(i)=gagetable{1,3}(i);
        gageobs(i)=gagetable{1,5}(i);
    end
%     figure;hold all; plot(gagetime,gageobs,'r.-',gagetime,gagepred,'b.');datetick('x','mm/yy')
%     xlabel('Date (Month/Year)');ylabel('Tidal height (m)')
%     set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]);
%     box on 
%       legend('NOAA observation','NOAA prediction')
%     
    
    %interpolate to image acquisitions; to do : gap in data
    tideh = interp1( gagetime , gageobs , SDtime,'linear') ;
    tideh_pred = interp1( gagetime , gagepred , SDtime,'linear') ;

    [SDtime_s,idsort]=sort(SDtime);
    tideh_s=tideh(idsort);
    tideh_pred_s=tideh_pred(idsort);

if flagplot==1
    figure;hold all; plot(SDtime_s,tideh_s,'r.-');datetick('x','mm/yy')
    plot(SDtime_s,tideh_pred_s,'b>-');datetick('x','mm/yy')
xlabel('Date (Month/Year)');ylabel('Tidal height (m)')
set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]);
box on
legend('NOAA observation','NOAA prediction')
end

end

%plot
if flagplotany==1
figure;plot(SDtime,tideh,'bo');datetick('x','mm/yy')
saveas(gcf,'tideallimages','fig')
end

%sort the date 
[SDtime_s,idsort]=sort(SDtime);
tideh_s=tideh(idsort);

if flagplotany==1
figure;hold all; plot(SDtime_s,tideh_s,'>-');datetick('x','mm/yy')
xlabel('Date (Month/Year)');ylabel('Tidal height (m)')
set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]);
saveas(gcf,'tideallimages','fig')
end

%save modeled tides 
listfile=['tides.txt'];
fidlist = fopen(listfile, 'w');
fprintf(fidlist,'%s \n','# Date (yyyymmddHHMMSS)   Modeled tidal heights (m)');
%sort isv SDtime tideh
for j=1:length(idsort) %[3,5,24,26] % 2012 Kamchatka Volcano 
        k=idsort(j);
        i=isv(k);
        ymd=f{i}(6:19);%datenum(ymd,'yyyymmddHHMMSS') same as SDtime(k), SDtime_s(j)
        fprintf(fidlist,'%s %12.2f \n',ymd, tideh_s(j));
end %j
fclose(fidlist);


% compare with gages
if 0
    latp1=64+29.7/60;lonp1=-(165+26.4/60)+360; %Nome gage
    [tideh_1,conList]=tmd_tide_pred(Model,SDtime,latp1,lonp1,'z');
    
%     lonp1=-164.655420912470+360;       latp1=63.0290610289717; %Yukon
% /Users/chunlidai/ArchivedBoxSync/SeaLevelST2019/runcoastline/yukoncheck5p.gmt 
    tideh_1_s=tideh_1(idsort);
    figure;hold all; plot(SDtime_s,tideh_1_s,'>-');datetick('x','mm/yy')
    xlabel('Date (Month/Year)');ylabel('Tidal height (m)')
    set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]);
  
    box on 
    legend('TPXO-9.1 Tidal height at Nome Gage');
  

end


if 0 %plot Fig. N8 
    
%plot the maximum range of tide heights at the tile region
minlon=min(XI(:,1));maxlon=max(XI(:,1));
minlat=min(XI(:,2));maxlat=max(XI(:,2));
idx=lon>minlon&lon<maxlon;idy=lat>minlat&lat<maxlat;  
SDtime3=max(SDtime)-366:(1/24):max(SDtime);
latp3=lat(idy);lonp3=lon(idx);M3=M(idy,idx);tiderange=zeros(size(M3));
LAT3=LAT(idy,idx);LON3=LON(idy,idx);
tic
for jj=1:length(LAT3(:))
[tideh3,conList]=tmd_tide_pred(Model,SDtime3,LAT3(jj),LON3(jj),'z');
% tiderange(jj)=max(tideh3)-min(tideh3);
tiderange(jj)=max(abs(tideh3));
end
toc
figure;
set(gcf,'Color','white')
set(gca,'FontSize', 8);
set(gcf, 'PaperPosition', [0 0 4 4]); 
set(gcf, 'PaperSize', [ 4 4]); 
hold all;
imagesc(lonp3-360,latp3,tiderange,'alphadata',M3)
colorbar;colormap jet
hl=xlabel('Longitude ($^{\circ}$)');
set(hl, 'Interpreter', 'latex');
hl=ylabel('Latitude ($^{\circ}$)');
set(hl, 'Interpreter', 'latex');
box on;
caxis([4 5])
    
%Time series at the most measured grid.
% k=[1456643];
ik=mode(k);
lonp5=-(156+57/60)+360; latp5=56+54/60+30/3600;%location of Fig. 9 
kid=find(abs(X(k,2)-latp5)<1/12&abs(X(k,1)-lonp5)<1/12);ik=mode(k(kid));
M=k==ik; %63 images near Nakalilok bay
latp1=X(ik,2);lonp1=X(ik,1); %-157.0000 56.8333
hold on;plot(lonm(M)-360,latm(M),'ko','markersize',6)
hold on;plot(lonp1-360,latp1,'w>','linewidm',4,'markersize',12)
axis square
axis([-157.4208 -156.3306 56.5509 57.1483])
%boundary of tile
bd=[
-157.0872   56.5730
 -157.3801   56.9665
 -156.6560   57.1261
 -156.3706   56.7303
 -157.0872   56.5730];
plot(bd(:,1),bd(:,2),'c-')
s10km=[-157.3333 56.5833
-157.1725 56.5833];
plot(s10km(:,1),s10km(:,2),'k-','linewidm',2)
text(-157.253-0.07, 56.610,'10 km')
saveas(gcf,'tiderange','fig')
print('-dpdf','-r300','tiderange') 

%Plot the estimated tidal height along the final coastlines.
load xoyo.mat
n=length(xo);pointsize = 10;hold all;
for k=1:n
    scatter(xo{k}, yo{k}, pointsize, zo{k});colorbar
end
caxis([-1 4])
    
%Fig. 9 b and c
% WV02_20130311223531_103001001FC82300_13MAR11223531-M1BS-500062155090_01_P001.xml
% WV02_20130414214424_10300100219D6B00_13APR14214424-M1BS-500071870100_01_P002.xml
timeb=datenum('20130414214424','yyyymmddHHMMSS'); %14 April 2013, 21:44 UTC%epoch of Fig. 9 b and c.
timec=datenum('20130311223531','yyyymmddHHMMSS'); %11 March 2013, 22:35 UTC
timea=datenum('20140925215150','yyyymmddHHMMSS');%WV02_20140925215150_103001003727E300_14SEP25215150-M1BS_R01C1-500106091130_02_P001.xml
% SDtime4=SDtime(M);
idb=find(abs(SDtime-timeb)<0.5/86400.); %151 %close to grid lonp5.
idc=find(abs(SDtime-timec)<1/86400.); %134 %closer to the grid left of lonp5.
ida=find(abs(SDtime-timea)<0.5/86400.);
hold on;plot(lonm([idb(:)])-360,latm([idb(:)]),'go','markersize',6)
hold on;plot(lonm([idc(:)])-360,latm([idc(:)]),'mo','markersize',6) 
fprintf(['\n Tidal height at 20130414214424 :',num2str(tideh(idb)),'\n']) %-1.051m
fprintf(['\n Tidal height at 20130311223531 :',num2str(tideh(idc)),'\n']) %1.6316m
tideh(ida) %0.93m

% SDtime1=min(SDtime):1/24.3:max(SDtime);%when interval is not an integer of hour.
SDtime1=floor(min(SDtime)):(1/24):max(SDtime);
[tideh1,conList]=tmd_tide_pred(Model,SDtime1,latp1,lonp1,'z');
figure;hold all;
plot(SDtime1,tideh1,'b.-');datetick('x','mm/yy')
plot(SDtime1(is:is+23),tideh1(is:is+23),'k>-');
hold on;plot(SDtime(M),tideh(M),'ro')
xlabel('Date');ylabel('Tidal height (m)')
saveas(gcf,'tideallimages1','fig')

%Mean at the hour
tm=0:23;tidem=zeros(length(tm),1);
x=(SDtime1-floor(SDtime1))*24;
for k=1:24
M2=abs(x-tm(k))<0.1;
tidem(k)=mean(tideh1(M2));
end

%Time series as hours of the day.
% sday=23+56/60+4/3600.;%a sidereal day lasts for 23 hours 56 minutes 4.091 seconds;
% SDtime1s=SDtime1*sday/24;
isb=find(abs(SDtime1-datenum('20130414000000','yyyymmddHHMMSS'))<0.5/24);
isc=find(abs(SDtime1-datenum('20130311000000','yyyymmddHHMMSS'))<0.5/24);
figure;
set(gcf,'Color','white')
set(gca,'FontSize', 8);
set(gcf, 'PaperPosition', [0 0 4 4]); 
set(gcf, 'PaperSize', [ 4 4]); 
hold all;
hold all;plot((SDtime1-floor(SDtime1))*24,tideh1,'b.');
plot(tm,tidem,'b>-','linewidm',1.5,'Markersize',6);
plot((SDtime1(isb:isb+23)-floor(SDtime1(isb:isb+23)))*24,tideh1(isb:isb+23),'m>-','linewidm',1.5,'Markersize',6);
plot((SDtime1(isc:isc+23)-floor(SDtime1(isc:isc+23)))*24,tideh1(isc:isc+23),'g>-','linewidm',1.5,'Markersize',6);
hold on;plot((SDtime(M)-floor(SDtime(M)))*24,tideh(M),'ko')
legend('Modeled time series of 13 years','Hourly mean tidal heights','14 April 2013','11 March 2013','All observations near this grid')
xlabel('Time (hours) of the day (GMT)');ylabel('Tidal height (m)')
box on
axis([0 24 -3 4])
print('-dpdf','-r300','tideallimages1hr') 
saveas(gcf,'tideallimages1hr','fig')

%Model at a certian time of the day.
x=(SDtime1-floor(SDtime1))*24;
M2=abs(x-12.71)<0.1;
figure;plot(x(M2),tideh1(M2),'b.');title('Tidal height (m) at 12.71 hour of all days')
figure;histogram(tideh1(M2),'Normalization','pdf','Facecolor','blue');
xlabel('Tidal height (m)');ylabel('PDF')
title('Tidal height (m) at 12.71 hour of all days')

figure;hold all;histogram(tideh1,'Normalization','pdf','Facecolor','blue');
histogram(tideh(M),'Normalization','pdf','Facecolor','red');
legend('Tidal model 1hr interval','Tidal model sampled by satellites')
xlabel('Tidal height (m)');ylabel('PDF')

figure;hold all;histogram(tideh1,'Normalization','cdf','Facecolor','blue');
histogram(tideh(M),'Normalization','cdf','Facecolor','red');
legend('Tidal model 1hr interval','Tidal model sampled by satellites')
xlabel('Tidal height (m)');ylabel('CDF')

end %if 0 %plot

end
