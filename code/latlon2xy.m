function [xj,yj]=latlon2xy(latj,lonj,projgdalj)
       constant

       %initialize 
       xj=nan(size(latj));yj=nan(size(latj));

       if ~exist('flagproj','var');flagproj=1;end
       % 2000 points take 96 seconds for flagproj=2, gdal
       if length(latj)>50;fprintf(['\n Input points are too many; use flagproj=1.  \n']); flagproj=1;end

       if strcmp(projgdalj,'epsg:3413'); % Arctic
	      [xj,yj]=polarstereo_fwd(latj,lonj,[],[],70,-45);
	   elseif strcmp(projgdalj,'epsg:3031') %Antarctica
	      [xj,yj]=polarstereo_fwd(latj,lonj,[],[],-71,0);
	   else %utm
	   
	   % e.g., UTM Zone 2S: EPSG 32702; UTM Zone 3N: EPSG 32603;
	       %Using gdal, e.g., echo $lon $lat | gdaltransform -s_srs EPSG:4326 -t_srs EPSG:32719
	       crs_epsg=str2double(projgdalj(6:end));%e.g.32636;

	       %       p1 = projcrs(crs_epsg); %projcrs, Since R2020b
	       
	       samplefile=['sample/sample',num2str(crs_epsg),'.tif'];

           try
            %flagproj 1, default using projinv function which may produce wrong results using Matlab2017b; 
            %         2, using gdaltransform 
               if flagproj==1
	               p1 = geotiffinfo(samplefile);
	               [xj,yj] = projfwd(p1, latj,lonj); 
	           % wrong results of projfwd when using Matlab R2017b for earthdem_23_pacific/WV02_20160120222555_103001005082AD00_16JAN20222555-M1BS-500805500030_01_P010_u16ns32731.xml.  Correct using R2020b.
	           % This image is a rare case when image located in UTM zone 60s but writen to zone 31S.
               elseif flagproj==2
                for k=1:length(latj)
		    if isnan(lonj(k))|isnan(latj(k));continue;end
                    %echo $lon $lat | gdaltransform -s_srs EPSG:4326 -t_srs EPSG:32719
                    str1=sprintf('echo %f %f | gdaltransform -s_srs EPSG:4326 -t_srs %s',lonj(k),latj(k),projgdalj);
                    [status , cmdout ]=system(str1);
                    tmp= sscanf(cmdout,'%f');
                    yj(k)=tmp(2);xj(k)=tmp(1);
                end                   
               end

	       catch e
		    fprintf(['\n latlon2xy.m There was an error! Input is:',num2str([latj(:)' lonj(:)']),projgdalj,'; The message was:',e.message]);
		    % this error may happen when input coordinates is too far away from the epsg zone, which is good to filter out these images later.
                    xj=nan; yj=nan;
           end

   	   end

return
end
