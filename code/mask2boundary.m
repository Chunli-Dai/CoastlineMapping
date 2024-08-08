function [Co]=mask2boundary(xout,yout,Modfil,Medgs1);
%Given water mask Modfil (logical; 1 water, 0 land) get its boundary.
%Modified from /home/dai.56/arcticdemapp/coastline/codec2/CoastTileMono.m
%input: Medgs1, a matrix (logical), 1 is the void area to be removed.
% xout, yout: polar stereographic coordinates of the matrix Modfil

constant  %get cloudarea

%boundary of box
Medgstb=false(size(Modfil));
Medgstb(:,1)=1;Medgstb(:,end)=1;Medgstb(1,:)=1;Medgstb(end,:)=1;%boundary of the box
Medgs1=Medgstb|Medgs1;

Co=[];
ofile1='boundary.shp';

resx=mean(xout(2:end)-xout(1:end-1));resy=mean(yout(2:end)-yout(1:end-1));
resr=mean([abs(resx),abs(resy)]);

nsuby=length(yout);nsubx=length(xout);

B = bwboundaries(~Modfil,'noholes'); %^_^ keep the unwanted lake islands.
% B = bwboundaries(~Modfil); % Delta area separated from ocean by road is treated as a like. Keep all large lakes.
n=length(B);xo=cell(n,1);yo=cell(n,1);
%clear prob Medgsib Medgstb
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
clear Mt BWb Medgs1
idd=find(cellfun(@isempty,xo)); %fix bug 12
xo(idd)=[];yo(idd)=[];

if 0
    out=[];
    for i=1:length(xo)
    out=[out;yo{i}(:), xo{i}(:)]; %lat lon
    end
    save -ascii outline.dat out

    %save xoyo.mat xo yo -v7.3
end

if 0
    figure;
    hold all
    for k=1:length(xo)
    plot(xo{k},yo{k},'.-')
    end
end

% shp1 = struct('Geometry', 'PolyGon', 'X', xo, 'Y', yo);
shp1 = struct('Geometry', 'PolyLine', 'X', xo, 'Y', yo);
%shapewrite(shp1, 'coastline.shp');
if ~isempty(shp1)
shapewrite(shp1, ofile1);
end

return
end
