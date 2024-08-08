function [XYbi,rangei,projgdali]=imagebd(ifile);
%get the Polygon boundary, and the rectangle range of a mono image from .xml file
% e.g. /data1/pgc_projects/dai_aleutians_multi_mono/imagery//WV02/WV02_20130923221601_1030010026BD6B00_13SEP23221601-M1BS-500127380110_01_P001.xml 

	constant 
	projgdali=''; % initialize; %need to update the retrieval of projgdal from xml or meta.txt mdf.txt.


	vstr={'<ULLON>','<ULLAT>','<URLON>','<URLAT>','<LRLON>','<LRLAT>','<LLLON>','<LLLAT>'};
    flagfmt1=1; %flag for format 1; default format
    
    vstr2={'<upperLeft>','<upperRight>','<lowerRight>','<lowerLeft>'};
    flagfmt2=1; %

	metafile=ifile;
	c=textread(metafile,'%s','delimiter','\n');
    
    %check the format 
    i=1;str=vstr(i);
    r=find(~cellfun(@isempty,strfind(c,str)));
    if isempty(r)
        flagfmt1=0;
    end
    i=1;str=vstr2(i);
    r=find(~cellfun(@isempty,strfind(c,str)));
    if isempty(r)
        flagfmt2=0;
    end

% Get the Footprint Vertices X, Y, close the loop
    if flagfmt1==1 
	n=length(vstr)/2;
	lon=zeros(n+1,1);lat=zeros(n+1,1);
	for i=1:n*2
	str=vstr(i);
	r=find(~cellfun(@isempty,strfind(c,str)));
	if isempty(r)        
		warning(['xml file is different as anticipated.',ifile])
        XYbi=[0 0]; 	rangei=[0 0 0 0];
		return;
	end %
	c2=c{r(1)};
	r1=strfind(c2,'>');r2=strfind(c2,'</');
    c2([1:r1(1),r2(1):end])='';
	z = sscanf(c2, '%g', 1);
	j=ceil(i/2);
	if mod(i,2)  %1 odd number, 0 even
           lon(j)=z;
	else
	   lat(j)=z;
	end
	end % if i

	%add projgdali
    str='<EPSG_CODE>';
    r=find(~cellfun(@isempty,strfind(c,str)));
    c2=c{r(1)};
	pattern = '<EPSG_CODE>(\d+)</EPSG_CODE>';
	matches = regexp(c2, pattern, 'tokens');
	% Check if there is a match
	if ~isempty(matches)
            % Extracted number as a string
            projgdali= matches{1}{1};
            projgdali=['epsg:',projgdali];
            
    else
            fprintf(['\n EPSG code not found in:',ifile,'\n']);
    end

    elseif flagfmt2 ==1
        %/data4/EarthDEM/alaska_2018oct26/qb_wv_alaska_metadata/QB02_20030408203848_1010010001C95401_03APR08203848-P1BS-000000073404_01_P001.xml 
    n=length(vstr2);
	lon=zeros(n+1,1);lat=zeros(n+1,1);
	for i=1:n
	str=vstr2(i);
	r=find(~cellfun(@isempty,strfind(c,str)));
	if isempty(r)        
		warning(['xml file is different as anticipated.',ifile])
        XYbi=[0 0]; 	rangei=[0 0 0 0];
		return;
	end %
	c2=c{r(1)+1};
	r1=strfind(c2,'>');r2=strfind(c2,'</');
    c2([1:r1(1),r2(1):end])='';
	zlat = sscanf(c2, '%g', 1);
    
    c2=c{r(1)+2};
	r1=strfind(c2,'>');r2=strfind(c2,'</');
    c2([1:r1(1),r2(1):end])='';
	zlon = sscanf(c2, '%g', 1);

       lon(i)=zlon;
	   lat(i)=zlat;
	end % for i        
        
    else 
        		warning(['xml file is different as anticipated.',ifile])
        XYbi=[0 0]; 	rangei=[0 0 0 0];
		return;
    end
    
	lon(n+1)=lon(1);lat(n+1)=lat(1); %close the loop
	%[Xb,Yb]=polarstereo_fwd(lat,lon,[], [],70,-45);
        [Xb,Yb]=latlon2xy(lat,lon,projgdali);

        XYbi=[Xb(:),Yb(:)]; %n by 2
	rangei=[min(Xb) max(Xb) min(Yb) max(Yb)];

	[demdir,filename,ext] =fileparts([strtrim(ifile)]);
        satname=filename(1);
	if strcmp(satname,'W')&& ( (max(Xb)-min(Xb))>50e3||(max(Yb)-min(Yb))>50e3)
	warning(['\n imagebd.m ',ifile,' file size is too big; Need to check. \n'])
	end
return
end
