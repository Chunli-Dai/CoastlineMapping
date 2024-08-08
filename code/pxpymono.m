function [dX4Sg]=pxpymono(dX4Sg2,f2is,fdir2is,XYb2is,fis,fdiris,XYbis)
%To collect the translational parameters of reference mono images that are components of strip files.
%And to coregsiter all mono image to the reference image
%dX4Sg2,f2is,fdir2is: Translational parameters, filenames, and directories of strip files.
% fis: filenames of mono images.
%output: Mrefmo : flag for fis that whether it is reference mono image or not.
%	psrefmo: the pzxy for the reference mono image.

%save test1.mat dX4Sg2 f2is fdir2is XYb2is fis fdiris XYbis -v7.3

n=length(fis);n2=length(f2is);
Mrefmo=false(n,1);
psrefmo=zeros(n,3);
flagref=zeros(n,1);

fprintf(['\n Collect all strip orthoimages as reference, total count:',num2str(n2),'\n'])
for j=1:n2
stripmetafile= [fdir2is{j},'/',f2is{j}];
[~,filename,~]=fileparts(stripmetafile);
satname=filename(1:4);

p3=dX4Sg2(j,:);

%copied from rivergithub2/multispecstrip.m
%get image 1 filenames, and dx2;
c=textread(stripmetafile,'%s','delimiter','\n');
rs=find(~cellfun(@isempty,strfind(c,'scene')));

if ~isempty(rs) %strip DEM
nsce=length(rs)-1; %number of scenes
rmsmeta=zeros(nsce,1);idd=[];dzxyd=zeros(nsce,3);mfile=cell(nsce,1);
for i=1:nsce
% c1=c{rs(1)+i};r1=strfind(c1,'dem.tif');c1(1:(r1+6))='';
%new release scene name is '_dem_smooth.tif'
c1=c{rs(1)+i};
matchStr = regexp(c1,'.tif','split');c1=matchStr{end};
tmp=sscanf(c1, '%g', 4);
rmsmeta(i)=tmp(1);dzxyd(i,1:3)=tmp(2:4);

c1=c{rs(1+i)+3};r1=strfind(c1,'Image 1=');if(isempty(r1)) Warning('Image 1 not found') ; end
r1=strfind(c1,'/');c1(1:(r1(end)))='';
satname=c1(1:4); %add this line, in case of strip file W1W1_20130916_1020010025278100_1020010026377B00_seg1_2m_meta.txt

if strcmp(satname,'WV01')
    mfile{i}=deblank(c1);
else
    mfile{i}=deblank(strrep(c1,'P1BS','M1BS'));%e.g.WV02_20160304214247_1030010052B75A00_16MAR04214247-M1BS-500641617080_01_P009.tif
end

end

else %rs is empty; scene DEM meta file
% compatible for scene DEM meta file.
check={'Image 1=','/','.tif'};
Mcheck=contains(c,check{1})&contains(c,check{2})&contains(c,check{3});
c1all=c(Mcheck);
nsce=sum(Mcheck(:));
rmsmeta=zeros(nsce,1);idd=[];dzxyd=zeros(nsce,3);mfile=cell(nsce,1);

for i=1:nsce
c1=c1all{i};
r1=strfind(c1,'/');c1(1:(r1(end)))='';

c11=deblank(c1);
satname=c11(1:4);%use this, since strip could be W1W2_20110825_1020010014850E00_103001000C5BD600_seg1_2m_meta.txt

if strcmp(satname,'WV01')
    mfile{i}=deblank(c1);
else
    mfile{i}=deblank(strrep(c1,'P1BS','M1BS'));%e.g.WV02_20160304214247_1030010052B75A00_16MAR04214247-M1BS-500641617080_01_P009.tif
end
end
end %if ~isempty(rs)

%find the unique images, delete the repeat file names, average the dzxy, and rms
[un idx_last idx] = unique(mfile,'stable');
unique_idx = accumarray(idx(:),(1:length(idx))',[],@(x) {(x)});
mfile=un;nsce=length(idx_last);
rmsmetau=zeros(nsce,1); dzxydu=zeros(nsce,3);
for i=1:nsce
rmsmetau(i)=mean(rmsmeta(unique_idx{i}));
dzxydu(i,1:3)=mean(dzxyd(unique_idx{i},:),1);
end
fprintf(['\n Strip file ',stripmetafile,' has ',num2str(nsce),'mono images.'])
fprintf(['\n They are:',mfile{1:end}])

for is=1:nsce
ntffile=strrep(mfile{is},'.tif','.ntf');
ntffile(end-3:end)='';%delete the extension.
ntffile=strrep(ntffile,'_temp',''); %new release
rs=find(~cellfun(@isempty,strfind(fis,ntffile)));
if ~isempty(rs) %if water mask file is found.
%   data=Msv(rs);
	p2=dzxydu(is,1:3);
	ps=p3+p2;
	Mrefmo(rs)=1;
	flagref(rs)=j;
%psrefmo(rs,1:3)=ps;%Not accurate, since mono image is not aligned with DEM.
	fprintf(['\n Strip id, mono id, ps (not accurate):',num2str([j, is, ps]),', ', ntffile])
else
    warning(['Image scene file ',ntffile,' not found!'])
end
end 
end % for j

%dX4Sg=psrefmo; %apply the collected mono image translational parameters.

fprintf(['\n Coregister all mono image to the reference images.'])
idtar=1:n;%find(Mrefmo~=1);
%idtar=[72,82]; %test two mono images (79 80) of strip WV02_20160827_103001005A735B00_103001005CAE7900_seg1_2m_dem.tif
fprintf(['\n Number of targe images:',num2str(length(idtar)),'\n'])
dX4Sg=zeros(length(idtar),3);
for j=1:length(idtar)
%fprintf(['\n hi1 j=',num2str(j)])
	is=idtar(j);
	tarimage=[fdiris{is},'/',fis{is}];

    Xb=XYbis{is}(:,1); Yb=XYbis{is}(:,2);
    rang0=[min(Xb) max(Xb) min(Yb) max(Yb)];
    dx=2;dy=dx;
    nwx=round((rang0(2)-rang0(1))/dx);nwy=round((rang0(4)-rang0(3))/dy);
    sx=Xb;sy=Yb;
    idx=round((sx-rang0(1))/dx)+1;idy=round((sy-rang0(3))/(dy))+1;
    sm0=poly2mask(idx,idy,nwy,nwx); % fast, apply to each polygon one by one.
%fprintf(['\n hi2 j=',num2str(j)])

	if Mrefmo(is)==1 %The mono image is part of one reference strip file, use the corresponding strip file.
	   isref=flagref(is);
	else
	%find the most overlapped reference image.
        nptall=zeros(n2,1);
        for jj=1:n2
          k=jj;
          sx=XYb2is{k}(:,1); sy=XYb2is{k}(:,2);
          idx=round((sx-rang0(1))/dx)+1;idy=round((sy-rang0(3))/(dy))+1;
          sm=poly2mask(idx,idy,nwy,nwx); % fast, apply to each polygon one by one.
          nptall(jj)=sum(sum(sm.*sm0))/sum(sum(sm0))*100; %
        end %k
        [maxratio,jjref]=max(nptall);
        isref=jjref;
        if (maxratio <= 2)
            warning(['For this target mono image the best overlapping reference strip covers only ',num2str(maxratio),'% of it']);
        end
	end %if
	refimage=[fdir2is{isref},'/',f2is{isref}];
    fprintf(['\n Mono image id, tar ref: ',num2str([is]),' ', tarimage,' ',refimage])

	odircoreg='./outcoreg/';

	[txy,flag]=coregmonostrip(refimage,tarimage,odircoreg);

	if flag==0 %coregistration failure
           fprintf(['\n Mono image coregistration failure id: ',num2str([is]),' ',tarimage,' ',filecpt,'\n'])
	else
	   ps=[0 txy]; %zxy
%	psref=psrefmo(isref,1:3);
	   psref=[dX4Sg2(isref,1:3)];psref(1)=0;
	   dX4Sg(is,1:3)=ps+psref;
	   fprintf(['\n Mono image id, pzxy, pzxytile: ',num2str([is, ps,ps+psref]),' ',tarimage,' ',refimage,'\n'])
	end  %  coregistration failure or not
end

return
end
