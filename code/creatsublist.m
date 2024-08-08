function [range,XYbg,projgdalg,f,fdir]=creatsublist(rang0,range,XYbg,projgdalg,f,fdir); % need rang0 as the input

%create file 'imagelist.dat' which include lists of images only within the selected tile.
constant

fprintf('\n Select images within the given tile, which takes:');
tic
x0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];y0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];
XYb0=[x0(:),y0(:)];
overlap_ratio1=zeros(length(fdir),1);
for j=1:length(fdir)
	projgdalj=projgdalg{j}; % e.g.,'epsg:32637'
	XYbj=XYbg{j};
	%
	if ~strcmp(projgdalj,projgdal);
        %  fprintf(['\n Image has different projection than that in the constant.m: ',ifile,', ' , projgdalj,' \n']);
        %  fprintf(['\n Converting the projection for the boundary. \n']);
	   xj=XYbj(:,1);yj=XYbj(:,2);
	%  xy to lat lon
	   [latj,lonj]=xy2latlon(xj,yj,projgdalj);
	   if isnan(latj) % if all components are nans.
		   overlap_ratio1(j)=0;
		   continue;
	   end

	%  lat lon to xy
           [xj,yj]=latlon2xy(latj,lonj,projgdal); 
	   XYbj=[xj(:),yj(:)];

        end

	%
	try
	[overlap_ratio1(j),~,~]=getoverlap(XYb0,XYbj,f{j});
	catch e
                    fprintf('\n creatsublist.m,%s, There was an error! The message was:\n%s',f{j},e.message);
                    overlap_ratio1(j)=0;
        end

end

thres=1; %20; 0 include all; 20, discard images covering <20% of the box.

%M1=overlap_ratio1<thres;
M1=~(overlap_ratio1>=thres);

XYbg(M1)=[]; f(M1)=[];fdir(M1)=[];projgdalg(M1)=[];
range(M1,:)=[];

%filename=[currentdir,'/imagelist.dat'];
filename=['./imagelist.dat'];
fid = fopen(filename,'w');
for j=1:length(fdir)
   ifile0=[fdir{j},'/',f{j}]; %xml
   ifile= strrep(ifile0,'.xml','.tif');

   projgdalj=projgdalg{j}; % e.g.,'epsg:32637'

   %convert the projection to make sure all images have the same projection as in constant.m.
   if ~strcmp(projgdalj,projgdal);
      ofile=[imagesubdir,'/',f{j}];
      ofile=strrep(ofile,'.xml','.tif');
      fprintf(['\n Image has different projection than that in the constant.m: ',ifile,', ' , projgdalj,' \n']);
      fprintf(['\n Converting the projection for the boundary, but no cropping. \n']);
      %update 
      XYbj=XYbg{j};
      xj=XYbj(:,1);yj=XYbj(:,2);
      %  xy to lat lon
      [latj,lonj]=xy2latlon(xj,yj,projgdalj);

      %  lat lon to xy
      [xj,yj]=latlon2xy(latj,lonj,projgdal);
      XYbj=[xj(:),yj(:)];
      XYbg{j}=XYbj; range(j,:)=[min(xj) max(xj) min(yj) max(yj)];

      rangib=range(j,:);

      testrb=num2str([rangib(1) rangib(3) rangib(2)  rangib(4)]); %;[-te xmin ymin xmax ymax]

      str=['gdalwarp ',ifile,' ',ofile,' -t_srs ',projgdal,' -te ', testrb,' -r bilinear ',' -ot UInt16']
      [status, cmdout]=system(str)

      % March 25, 2024: copy xml file
      str=['cp ',ifile0,' ', imagesubdir] 
      [status, cmdout]=system(str);

      % update file names
      fdir{j}=imagesubdir; projgdalg{j}=projgdal;
   end

   fprintf(fid,'%s \n',ifile0);
end
fclose(fid);

toc

return
end
