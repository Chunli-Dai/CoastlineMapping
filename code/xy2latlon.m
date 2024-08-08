function [latj,lonj]=xy2latlon(xj,yj,projgdalj)
% Check
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
%

%Arctic
% [x2,y2]=latlon2xy(79.9641229,-99.7495626,'epsg:3413');
%num2str([x2, y2]); '-890000.0922     -629000.0661'
%[lat,lon]=xy2latlon(x2,y2,'epsg:3413');
%num2str([lat, lon],'%.15f'); '79.964122899998770 -99.749562600000004'
% Compare with: echo lon lat | gdaltransform -s_srs EPSG:4326 -t_srs EPSG:3413
%  xy= -890000.092229265 -629000.066181033  % difference 3e-5m 8e-5m

%Antarctic
%[x2,y2]=latlon2xy(-75,150,'epsg:3031');
%num2str([x2, y2]); '819391.61915     -1419227.9157'
%[lat,lon]=xy2latlon(x2,y2,'epsg:3031');
%num2str([lat, lon],'%.15f %.15f');  '-74.999999999997712 150.000000000000028'
% Compare with: echo lon lat | gdaltransform -s_srs EPSG:4326 -t_srs EPSG:3031
%  xy= 819391.619203618 -1419227.9157568  % difference 5e-5 meter
    
       constant

       if ~exist('flagproj','var');flagproj=1;end

       % 2000 points take 96 seconds for flagproj=2, gdal
       if length(xj)>50;fprintf(['\n Input points are too many; use flagproj=1.  \n']); flagproj=1;end

       %initialize
       latj=nan(size(xj));lonj=nan(size(xj));


	   if strcmp(projgdalj,'epsg:3413'); % Arctic
	      [latj,lonj]=polarstereo_inv(xj,yj,[],[],70,-45);
	   elseif strcmp(projgdalj,'epsg:3031') %Antarctica
	      [latj,lonj]=polarstereo_inv(xj,yj,[],[],-71,0);
	   else %utm
	   
	       %e.g., UTM Zone 2S: EPSG 32702; UTM Zone 3N: EPSG 32603;
	       %Using gdal, e.g., echo 747657 4031478 | gdaltransform -s_srs EPSG:32636 -t_srs EPSG:4326
	       % gdal too slow: % 5 points take 0.364243 seconds.

	       if isempty(projgdalj)|isnan(projgdalj)
		       latj=nan; lonj=nan;
		       return
	       end

	       crs_epsg=str2double(projgdalj(6:end));%e.g.32636;
	       if isnan(crs_epsg)
	          fprintf (['\n xy2latlon.m, crs_epsg= ',num2str(crs_epsg),', projgdalj=',projgdalj,'. \n'])
		  latj=nan; lonj=nan; return
       	       end

%       p1 = projcrs(crs_epsg); %projcrs, Since R2020b

		samplefile=['sample/sample',num2str(crs_epsg),'.tif'];
		try
            %flagproj 1, default using projinv function which may produce wrong results using Matlab2017b when image is too far from the map zone. 
            %         2, using gdaltransform 
            if flagproj==1
		    p1 = geotiffinfo(samplefile);
	        [latj, lonj] = projinv(p1, xj,yj);
            elseif flagproj==2
                for k=1:length(xj)
		    if isnan(xj(k))|isnan(yj(k));continue;end
                    %echo 747657 4031478 | gdaltransform -s_srs EPSG:32636 -t_srs EPSG:4326
                    clear str1
                    str1=sprintf('echo %f %f | gdaltransform -s_srs %s -t_srs EPSG:4326',xj(k),yj(k),projgdalj);
                    [status , cmdout ]=system(str1);
                    tmp= sscanf(cmdout,'%f');
                    latj(k)=tmp(2);lonj(k)=tmp(1);
                end
                
            end

	        catch e
%                   fprintf('\n xy2latlon.m There was an error! The message was:\n%s',e.message);
		    fprintf(['\n xy2latlon.m There was an error! Input is:',num2str([xj(:)' yj(:)']),projgdalj,'; The message was:',e.message]);

                    latj=nan; lonj=nan;
            end

   	   end

return
end
