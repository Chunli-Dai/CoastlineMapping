function writeGeotiff(OutFileName,x,y,z,fmt,nodata,projstr)

gdalpath =[];
if ismac
%     gdalpath = '/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/';
    gdalpath = './'; % chunli dai, use this line
end

%if ismac
%	tempfile =  [tempname(tempdir),'.envi'];
%else
%	 tempfile =  [tempname('/scratch/tmp'),'.envi'];
%end

outdir=OutFileName(1:find(OutFileName=='/',1,'last'));
tempfile =  [tempname(outdir),'.envi'];


if strcmp(projstr,'polar stereo north')
    enviwrite(tempfile,x,y(:),z,'format',fmt,'proj',projstr);
elseif strcmp(projstr,'polar stereo south')
    enviwrite(tempfile,x,y(:),z,'format',fmt,'proj',projstr);
    %enviwrite(tempfile,x,y(:),z,'format',fmt,'proj','polar stereo south');
elseif strcmp(projstr(1:3),'UTM') %projstr='UTM zone 45 north'
%     enviwrite(tempfile,x,y(:),z,'format',fmt,'proj','UTM','hemi','south','zone',23);
    newstr=split(projstr)
    enviwrite(tempfile,x,y(:),z,'format',fmt,'proj','UTM','hemi',newstr{4},'zone',str2num(newstr{3}));
end

% confirm file write
if ~(exist(tempfile,'file') && exist([tempfile,'.hdr'],'file'))
    disp('Warning, envi file not written in getImage, stopping')
    keyboard
end

system([gdalpath ,'gdal_translate -co bigtiff=if_safer -co compress=lzw -co tiled=yes -a_nodata ',...
    num2str(nodata),' ',tempfile,' ', OutFileName]);

delete([tempfile,'*'])
