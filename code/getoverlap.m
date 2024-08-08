function [ratio1, ratio2, ratiom]=getoverlap(XYb1,XYb2,file2)
%Given two polygons,
% Calculate the overlapping ratio wrt. the polygon 1, polygon2, and mean
% ratio.
% XYbi=[Xb,Yb]; %n by 2
    
    resr10m=100;

    %remove NaNs
    Mn=any(isnan(XYb1),2);XYb1(Mn,:)=[];
    Mn=any(isnan(XYb2),2);XYb2(Mn,:)=[];
    [n1,m1]=size(XYb1);[n2,m2]=size(XYb2);
    if n1<3||n2<3
        fprintf('\n getoverlap.m input polygon has less than 3 vertices. \n')
        ratio1=0;ratio2=0;ratiom=0;return
    end

    %Matrix can be too big when two polygons are far apart. Use only one polygon
    %rang0=[min([XYb1(:,1);XYb2(:,1)]), max([XYb1(:,1);XYb2(:,1)]), min([XYb1(:,2);XYb2(:,2)]),max([XYb1(:,2);XYb2(:,2)])];
    rang0=[min([XYb1(:,1)]), max([XYb1(:,1)]), min([XYb1(:,2)]),max([XYb1(:,2)])]; %Only ratio1 is valid.
    
    ranget=round(rang0/resr10m)*resr10m;rang0=ranget;

    xout=rang0(1):resr10m:rang0(2);yout=rang0(4):-resr10m:rang0(3);
    nx0=length(xout);ny0=length(yout);
    %check the size of data0. 
    nmin=min([nx0,ny0]);nmax=max([nx0,ny0]);
    if nmin<=2||nmax>1000
	    %Image may be too big because of projection distortion when using a zone too far away.
        fprintf(['\n getoverlap.m input polygon size is too small or too big : [nx0 ny0]=',num2str([nx0 ny0]),'. \n'])
        ratio1=0;ratio2=0;ratiom=0;return
    end
    data0.x=xout;data0.y=yout;data0.z=nan(ny0,nx0); %initial matrix; 

    %ratio in percentage.

    idx=round((XYb1(:,1)-data0.x(1))/resr10m)+1;idy=round((XYb1(:,2)-data0.y(1))/(-resr10m))+1;
    sm1=poly2mask(idx,idy,ny0,nx0); % fast, apply to each polygon one by one.
    idx=round((XYb2(:,1)-data0.x(1))/resr10m)+1;idy=round((XYb2(:,2)-data0.y(1))/(-resr10m))+1;
    sm2=poly2mask(idx,idy,ny0,nx0); % fast, apply to each polygon one by one.
    %sm:  1 has data; 0 has no dem data.
    overlap=sm1&sm2;

    ratio1=sum(overlap(:))/sum(sm1(:))*100; %percentage of tile pixels covering the input box.

    %% ratio 2
    rang0=[min([XYb2(:,1)]), max([XYb2(:,1)]), min([XYb2(:,2)]),max([XYb2(:,2)])]; %Only ratio2 is valid.

    ranget=round(rang0/resr10m)*resr10m;rang0=ranget;

    xout=rang0(1):resr10m:rang0(2);yout=rang0(4):-resr10m:rang0(3);
    nx0=length(xout);ny0=length(yout);
    %check the size of data0. 
    nmin=min([nx0,ny0]);nmax=max([nx0,ny0]);
    if nmin<=2||nmax>1000
        fprintf(['\n getoverlap.m ',file2,'; Input polygon size is too small or too big: [nx0 ny0]=',num2str([nx0 ny0]),'. \n'])
        ratio1=0;ratio2=0;ratiom=0;return
    end
    data0.x=xout;data0.y=yout;data0.z=nan(ny0,nx0); %initial matrix; 

    %ratio in percentage.

    idx=round((XYb1(:,1)-data0.x(1))/resr10m)+1;idy=round((XYb1(:,2)-data0.y(1))/(-resr10m))+1;
    sm1=poly2mask(idx,idy,ny0,nx0); % fast, apply to each polygon one by one.
    idx=round((XYb2(:,1)-data0.x(1))/resr10m)+1;idy=round((XYb2(:,2)-data0.y(1))/(-resr10m))+1;
    sm2=poly2mask(idx,idy,ny0,nx0); % fast, apply to each polygon one by one.
    %sm:  1 has data; 0 has no dem data.
    overlap=sm1&sm2;
    ratio2=sum(overlap(:))/sum(sm2(:))*100; %percentage of tile pixels covering the input box.

    ratiom=(ratio1+ratio2)/2;

    % figure;imagesc(data0.x*1e-3,data0.y*1e-3,overlap);colorbar;hold on; plot(XYb1(:,1)*1e-3,XYb1(:,2)*1e-3,'r>-',XYb2(:,1)*1e-3,XYb2(:,2)*1e-3,'g>-')

return
end
