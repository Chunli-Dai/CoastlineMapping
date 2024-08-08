function rang0=getbox(tilefile); 
% get boundary from given tile name.
    constant
    width=0; %tile boundary with buffer width in meter.

    str1=tilefile(1:3);
    [dir,ifile0,ext] =fileparts(tilefile);
    ifile=ifile0;

%   dxs=subtiledx*1e3; %25 km;

    if strcmp(str1,'utm') %earthdem tiles e.g., utm10n_47_05_2_2_01_02
	dx=100e3; x0=150e3; % y0=0;
	ifile=ifile0(8:end); %update ifile
	nsstr=ifile0(6);
	
	%https://www.pgc.umn.edu/guides/stereo-derived-elevation-models/pgc-dem-products-arcticdem-rema-and-earthdem/
	if strcmp(nsstr,'n') % north
	        y0=0;
	elseif strcmp(nsstr,'s')
		y0=3300000;
	else
		fprintf(['\n getbox.m Projection is not as expected. ',tilefile,' ',projgdal,' ', nsstr,'\n'])
	end

    elseif strcmp(projgdal,'epsg:3413')  % arcticdem tiles, e.g., 51_08_2_2_01_02

    dx=100e3;x0=-4000e3;y0=-4000e3;%xe=3400e3;ye=4000e3; %ArcticDEM Mosaic tiles coordinate reference;
    %dxs=dx/2/10; %5km
    %	    dxs=dx/50; %2km
    elseif strcmp(projgdal,'epsg:3031')  %Antarctica
	%   dx=100e3; x0=-4000e3;y0=-4000e3; %1000e3 larger than REMA tiles, e.g., 40_09_1_1
	    dx=100e3; x0=-4000e3+1000e3;y0=-4000e3+1000e3; 
    else
	    fprintf(['\n getbox.m Projection is not as expected. ',tilefile,' ',projgdal,'\n'])
    end

    r=1;
    yid= sscanf(ifile(r:(r+1)), '%g', 1);
    xid= sscanf(ifile(r+3:(r+4)), '%g', 1);
    yids= sscanf(ifile(r+6), '%g', 1);
    xids= sscanf(ifile(r+8), '%g', 1);
%   yidss= sscanf(ifile(r+10:r+11), '%g', 1);
%   xidss= sscanf(ifile(r+13:r+14), '%g', 1);

%   x=x0+(xid-1)*dx+(xids-1)*dx/2+(xidss-1)*dxs;y=y0+(yid-1)*dx+(yids-1)*dx/2+(yidss-1)*dxs;
%   rang0=[x-width x+dxs+width y-width y+dxs+width]; %tile boundary with buffer width

    x=x0+(xid-1)*dx+(xids-1)*dx/2;y=y0+(yid-1)*dx+(yids-1)*dx/2;
    rang0=[x-width x+dx/2+width y-width y+dx/2+width]; %tile boundary with buffer width
    
    ranget=round(rang0/resr)*resr;rang0=ranget;

return
end
