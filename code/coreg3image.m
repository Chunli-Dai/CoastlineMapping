function [dX4Sg,idregion,data0r]=coreg3image(rang0,idregion,XYbg,f,fdir,tag)
%Modified from /home/dai.56/arcticdemapp/landslide/code1/coreg3.m
%Input:
%       rang0 : the boundary of target zone;
%       idregion: index of strip DEMs within the target zone;
%       XYbg:     the boundary of the actual data of each strip DEM;
%       f,fdir:   the filename and directory of all available strip DEMs;
%       tag: a priori land surfaces.
%
%Modification 

close all
constant

odircoreg='./outcoregortho/';
if ~exist(odircoreg,'dir')
    mkdir(odircoreg)
end

params.I = [1 0 0 0 0 0 0];
params.C = [50 0.001 0.001 0.05 0.0001];
%         params.G = [3000 20];
params.M = '3D';
params.V = [10 20 10];
params.G = [9000 20]; %Adjust the max height parameter for the 2012 Kamchatka Volcano
%coregflag=8;%1 parameter (vertical), 3 parameter (Ian's); 7 parameter (MJ's); 8, MJ's setsm coregistration with 3 parameters.
%flagplot=0;

% find the DEM mosaics within rang0.
%xr0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];yr0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];

resrc=40;
%resr=2;%10;%40; %2; % 2 m cause out of memory for 1.7GB DEM strip for 3 parameters (Nuth and Kabb) coregistration.


fprintf ('\n Step 1.1b: Coregister all pairs of orthoimages.')
tic

pg=zeros(length(idregion),3);rmsreg2=zeros(length(idregion),1);
%dX4Sg=pg;
dX4Sg=nan(length(idregion),3);
if coregflag==7
%   clear dX4Sg
end
idd=[];idd2=[];

ndem=length(idregion);
count=0; %record the count of coregistered pairs;
clear offsets
for idem1=1:ndem %index of DEM 1

    %read DEM 1
    i=idregion(idem1);
    XYbi=XYbg{i};
    metafile=[fdir{i},'/',f{i}]; metafileref=metafile;
    sat=filename2sat(metafile); ymd=filename2ymd(metafile);
    textref=[sat,'_',ymd];

    for idem2=(idem1+1):ndem
        j=idregion(idem2);  
        XYbj=XYbg{j};
    
        fprintf(['\n Coregistering pair ',num2str([idem1  idem2]),': ', f{i}, ' vs ',f{j},'.\n'])
    
    %Add constraint: for each pair, the overlapping area has to be > 1% of each scene DEM.
    %all space of two DEMs;
    range1=[min([XYbi(:,1);XYbj(:,1);]) max([XYbi(:,1);XYbj(:,1);]) min([XYbi(:,2);XYbj(:,2);]) max([XYbi(:,2);XYbj(:,2);])];
    ranget=round(range1/resrc)*resrc;
    tx=ranget(1):resrc:ranget(2);ty=ranget(4):-resrc:ranget(3);

    idx=round((XYbi(:,1)-tx(1))/resrc)+1;
    idy=round(-(XYbi(:,2)-ty(1))/resrc)+1;
    inref= poly2mask(idx,idy, length(ty),length(tx)); % faster than inpolygon
    arearef=(sum(inref(:)));%pixels

    idx=round((XYbj(:,1)-tx(1))/resrc)+1;
    idy=round(-(XYbj(:,2)-ty(1))/resrc)+1;
    intar= poly2mask(idx,idy, length(ty),length(tx)); % faster than inpolygon
    areatar=(sum(intar(:)));%pixels
    overlaparea=sum(sum(inref&intar));
    ratioref=overlaparea/arearef*100; 
    ratiotar=overlaparea/areatar*100; 
    %e.g. 2km length overlap for a 120 km long strip and 16 km long scene.
    if min([ratioref ratiotar])> 1 && mean([ratioref ratiotar])> 5 %good case
    else
        fprintf(['\n Overlapping area of this pair is too small:',num2str([ratioref ratiotar]),' %.',]);
        continue;
    end
    
    metafile=[fdir{j},'/',f{j}];
    sat=filename2sat(metafile); ymd=filename2ymd(metafile);
    texttar=[sat,'_',ymd];

    %For faster computation (compatible with coregflag8 when skipping the generation of datatarz )
    % Crop out non rock area
    Mrang0=false(size(inref));
    idx=find(tx>=rang0(1) & tx<=rang0(2));
    idy=find(ty>=rang0(3) & ty<=rang0(4));
    Mrang0(idy,idx)=true;
    if exist('tag','var')&&flagrock==1 %must use rock as control points
      rocktag1 = interp2(tag.x,tag.y,tag.z,tx,ty','*nearest',0); %logical
      M1=inref&intar&logical(rocktag1)&Mrang0;
    else
      M1=inref&intar&Mrang0;
    end
    [X,Y]=meshgrid(tx,ty);
    xmin=min(X(M1==1));xmax=max(X(M1==1));
    ymin=min(Y(M1==1));ymax=max(Y(M1==1));
    rangeov=[xmin xmax ymin ymax];rangeov=round(rangeov/resr)*resr;
    if length(rangeov) < 4
       fprintf(['\n No Overlapping. '])
       continue;
    end

      %get rocktag
    if exist('tag','var')
        if ~isempty(tag)
    %       rocktag = interp2(tag.x,tag.y,tag.z,refdem.x,refdem.y','*nearest',0); %logical
	    idx=find(tag.x>=rangeov(1) & tag.x<=rangeov(2));
	    idy=find(tag.y>=rangeov(3) & tag.y<=rangeov(4));
	    rocktag=tag.z(idy,idx);

            rocktag = logical(rocktag);
        else
            rocktag=[];
        end
    else
        rocktag=[]; %for mp1 in coregisterdems.m
    end
    % 

    %For faster computation (skip the following lines with readdem, but compatible with old version)
    if coregflag==8 && flagplot==0
	datatarz=[]; tardem=[];refdem=[];
    else %get datar datatarz tardem refdem

    end %if coregflag==8 && flagplot==0

    iter= 1;%initialize
    perrdef=[0.1 0.1 0.1]/2; %[4 4 4];%the default sigma0 is 4 m same as in adjustOffsets.m
    perr=perrdef;%[0.01 0.01 0.01];
    p=[NaN NaN NaN];
    %output:z2out, cpts,dz;p,perr;
        
    if coregflag==8
        
        odircoregi=[deblank(odircoreg),'/i',num2str(i),'j',num2str(j),'/'];
        if ~exist(odircoregi,'dir')
          mkdir(odircoregi)
        end
        
	%to comment these
        if strcmp(metafileref(end-7:end),'meta.txt')
        refimage=strrep(metafileref,'meta.txt',demext);
        elseif strcmp(metafileref(end-7:end),'_mdf.txt')
        refimage=strrep(metafileref,'mdf.txt',demext);
        end
        if strcmp(metafile(end-7:end),'meta.txt')
        tarimage=strrep(metafile,'meta.txt',demext);
        elseif strcmp(metafile(end-7:end),'_mdf.txt')
        tarimage=strrep(metafile,'mdf.txt',demext);
        end
	
	tarimage=metafile;refimage=metafileref;

    fprintf(['\n Mono image id, tar ref: ',' ', tarimage,' ',refimage])
        
        %coregfile=[odircoregi,'DEM_coreg_result.txt'];
        coregfile=[odircoregi,'coreg_result.txt'];

        if exist(coregfile,'file')
        %in case files are already genearated in the previous run
        %read coregistration results. 

        clear in
	in.flag=1; %1 ortho; not exist: dem
        in.idem2=idem2;   in.odircoregi=odircoregi;     in.tarimagep=tarimage;
        in.datatarz=datatarz; in.tardem=tardem;in.refdem=refdem;
        %[p,perr,dzmms,dz,cpts,iter]=readcoregsetsm(in);
        [p,~,~,~,~,iter]=readcoregsetsm(in);
        
        else %run the coregistration again;

        if isempty(rocktag)||sum(rocktag(:))<10
            if flagrock==1 %must use rock as control points
                warning(['There is no rock surface as control points:',metafile]); p=[NaN NaN NaN];perr=perrdef;
                iter=49; continue
            else %for change detection, it's okay if there is no rock surface as control points.
		try
%               [refimagep,tarimagep,~,~]=prepareMJ(refimage,tarimage,odircoregi,[],rangeov);

	        catch e
                fprintf('prepareMJ.m There was an error! The message was:\n%s',e.message);
	        iter=49; continue
	        end
            end
        else
	    try
%           [refimagep,tarimagep,~,~]=prepareMJ(refimage,tarimage,odircoregi,[],rangeov,tag);
	    catch e
            fprintf('prepareMJ.m There was an error! The message was:\n%s',e.message);
	    iter=49; continue
	    end
        end

        [txy,flag]=coregmonostrip(refimage,tarimage,odircoregi,tag);
        %	txy=txy_org.txy; rho=txy_org.rho;

        if flag==0 %coregistration failure
           fprintf(['\n Mono image coregistration failure id: ',tarimage,' \n'])
        else
           p=[0 txy]; %zxy
           fprintf(['\n Mono image id, pzxy: ',num2str([ p ]),' ',tarimage,' ',refimage,'\n'])
        end  %  coregistration failure or not
        
        end % if file exist; read coregistration results -done.

	%sigma=dzmms(3);
        %if  ~(sigma<=maxsigma) || isnan(sigma) || ~(max(abs(p(2:3))) < maxpxpy) % in case there are NaNs in p.
        if  ~(max(abs(p(2:3))) < maxpxpy) % 
            warning(['coregistration failure:',metafile]); 
	    fprintf(['\n pzxy= ',num2str([p(:)]'),'.\n'])
	    p=zeros(3,1);%perr=perrdef;
            iter=49;   
        end
    
	end % if coregflag
    	fprintf(['\n Coregistration method:',num2str(coregflag),'; pzxy=',num2str(p(:)'),'.\n'])

     %check the quality
    if iter==49
%         idd=[idd;j];
        continue
    end
    
    %perr is for the dispersion matrix; can't be zeros; otherwise; the result would be all NAN;
    % regression failure in coregisterdems.m
    if any(perr==0) ||any(isnan(perr))
        M1=(perr==0);
        perr(M1)=perrdef(M1); %use the default std
        continue %choose to skip this results;
    end
    
    %save results to offsets
	dzmms=perrdef;%hi

    count=count+1;
    offsets.i(count)=idem1; %dem 1 index
    offsets.j(count)=idem2; %dem 2 index
    offsets.dz(count)=-p(1); %[n×1 double] z offset (dem 2 - dem 1) -> should be dem1 -dem2; otherwise results are wrong!
    offsets.dx(count)=-p(2); % [n×1 double] x offset (dem 2 - dem 1) -> see testadjustoffsets.m.
    offsets.dy(count)=-p(3); % [n×1 double] y offset (dem 2 - dem 1)
    offsets.dze(count)=perr(1); % [n×1 double] z offset 1-sigma error %has to be non zero, since this is for the dispersion matrix;
    offsets.dxe(count)=perr(2); % [n×1 double] x offset 1-sigma error
    offsets.dye(count)=perr(3); % [n×1 double] y offset 1-sigma error 
    offsets.mean_dz_coreg(count)=dzmms(1); % [n×1 double] mean diff in z after corgestration
    offsets.median_dz_coreg(count)=dzmms(2); % [n×1 double] median diff in z after corgestration
    offsets.sigma_dz_coreg(count)=dzmms(3); % [n×1 double] std dev of diff in z after corgestration
		

    if 0 %flagplot==1

	close all

        % The case when coregistration is not applied.
        z2n = interp2(tardem.x' ,tardem.y,double(datatarz) ,refdem.x',refdem.y,'*linear');
        nsr=round(resrc/resr);%1; %plot low resolution data
        
        [X,Y]=meshgrid(refdem.x,refdem.y);
        [LATa,LONa]=polarstereo_inv(X,Y,[],[],70,-45);
        if ~isempty(cpts)
            [LATc,LONc]=polarstereo_inv(cpts(:,1),cpts(:,2),[],[],70,-45);
        else
            LATc=[]; LONc=[];
        end

        figure;set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 6 2.5]);
        % imagesc(data0r.x*1e-3,data0r.y*1e-3,z2n-data0r.z);caxis([-5 5]);colorbar
        surf(LONa(1:nsr:end,1:nsr:end),LATa(1:nsr:end,1:nsr:end),z2n(1:nsr:end,1:nsr:end)-double(datarefz(1:nsr:end,1:nsr:end))); shading interp;
        colorbar;colormap jet;view(0,90)
        hl=xlabel('Longitude ($^{\circ}$)');
        set(hl, 'Interpreter', 'latex');
        hl=ylabel('Latitude ($^{\circ}$)');
        set(hl, 'Interpreter', 'latex');
        caxis([-5 5])
        title([texttar,' - ',textref,'; NO coregistration'])
        ofile=[texttar,'m',textref,'Nocoreg'];
        print('-dpng','-r400',ofile) 
        try
        saveas(gcf,ofile,'fig')
        end

        figure;set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 6 2.5]);
%         surf(LONa(1:nsr:end,1:nsr:end),LATa(1:nsr:end,1:nsr:end),z2out(1:nsr:end,1:nsr:end)-double(datarefz(1:nsr:end,1:nsr:end))); shading interp;
        surf(LONa(1:nsr:end,1:nsr:end),LATa(1:nsr:end,1:nsr:end),dz(1:nsr:end,1:nsr:end)); shading interp;
        % imagesc(data0r.x*1e-3,data0r.y*1e-3,z2out-data0r.z);caxis([-5 5]);colorbar
        % title('2011/10/08-2013/05/26 DEM (m); After coregistration')
        title([texttar,' - ',textref,'; After coregistration'])
%         plot(X(rockfilter)*1e-3,Y(rockfilter)*1e-3,'k.')
        hold on;plot(LONc,LATc,'k.','Markersize',0.1,'Linewidm',2) %Plot control points 
        colorbar;colormap jet;view(0,90)
        hl=xlabel('Longitude ($^{\circ}$)');
        set(hl, 'Interpreter', 'latex');
        hl=ylabel('Latitude ($^{\circ}$)');
        set(hl, 'Interpreter', 'latex');
        caxis([-8 8])
        ofile=[texttar,'m',textref,'wcoreg'];
        print('-dpng','-r400',ofile) 
        try
        saveas(gcf,ofile,'fig')
        end
        
        %Plot control points %see  scripts/plotsrtmvs.m
        if coregflag==7
        
            if 0 %plot hillshade
        tmpy=reshape(out{1,1}.y,size(out{1,1}.z));tmpx=reshape(out{1,1}.x,size(out{1,1}.z));
        datadsx=tmpx(1,:)';datadsy=tmpy(:,1); % to check the space
        [LAT,LON]=polarstereo_inv(tmpx,tmpy,[],[],70,-45);
        %control points
%         [LATc,LONc]=polarstereo_inv(out{1,2}.x(~out{1,2}.z(:)),out{1,2}.y(~out{1,2}.z(:)),[],[],70,-45);
        [LATc,LONc]=polarstereo_inv(cpts(:,1),cpts(:,2),[],[],70,-45);
        %interpolate to regular lon lat mesh grids for PLOTTING FIGURES
        lonmesh=linspace(min(LON(:)),max(LON(:)),length(datadsx));
        latmesh=linspace(min(LAT(:)),max(LAT(:)),length(datadsy));
        [LONmesh,LATmesh]=meshgrid(lonmesh,latmesh);
        [xm,ym]=polarstereo_fwd(LATmesh,LONmesh,[],[],70,-45);
        tz=double(out{1,3}.z);tz(tz== -9999) = NaN; 
        demrefmesh= interp2(tmpx,tmpy,tz,xm,ym,'linear',nan);
        demrefmesh(isnan(demrefmesh))=-9999;

        % plot control points in reference DEM hillshade
        hills=hillshade(double(out{1,3}.z),datadsx,datadsy,'plotit');
        set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.3 0 3.8 3]); set(gcf, 'PaperSize', [4 3]);
        axis equal
        xlabel('{\itx} (m)')
        ylabel('y (m)')
        % title(['DEM ']);
        hold on
        plot(cpts(:,1),cpts(:,2),'r.','Markersize',2,'Linewidm',2)
        % plot(out{1,1}.x(~out{1,2}.z(:))*1e-3,out{1,1}.y(~out{1,2}.z(:))*1e-3,'r>','Linewidm',2)
        % saveas(gcf,'DEM','fig')
        ofile=[texttar,'',textref,'cpts'];
        print('-dpng','-r400',ofile) 
        try
        saveas(gcf,ofile,'fig')
        end

        % plot control points in reference DEM hillshade in lon lat
        % hills=hillshade(demrefmesh,lonmesh,latmesh,'plotit','altitude',1,'azimuth',10);
        hills=hillshade(demrefmesh,lonmesh,latmesh,'azimuth',10,'altitude',10,'plotit');
        set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.3 0 3.8 3]); set(gcf, 'PaperSize', [4 3]);
        axis square
        hl=xlabel('Longitude ($^{\circ}$)');
        set(hl, 'Interpreter', 'latex');
        hl=ylabel('Latitude ($^{\circ}$)');
        set(hl, 'Interpreter', 'latex');
        hold on
        % plot(LONc,LATc,'r.','Markersize',20,'Linewidm',2)
        plot(LONc,LATc,'r.','Markersize',0.1,'Linewidm',2)
        ofile=[texttar,'',textref,'cptslat'];
        print('-dpng','-r400',ofile) 
        try
        saveas(gcf,ofile,'fig')
        end
            end %if 0 %plot hillshade

        %control surfaces histogram; dh is distance instead of height difference
        dh=out{1,4}.S{1,4}.height_diff;
        dhaf=out{1,4}.S{1,5}.height_diff;
        figure;
        hold all
        set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.3 0 3.8 3]); set(gcf, 'PaperSize', [4 3]);
        % edges=[-15:1:10];
        edges=[-15:0.1:10];
        % histo1=histogram(dh,edges,'Normalization','pdf');
        histo1=histogram(dh,'Normalization','pdf');
        histo=histogram(dhaf,'Normalization','pdf');
        histo.FaceColor='r';
        alpha(histo,.7)
        legend('Before Coregistration','After Coregistration','Location','NorthWest')
        % legend('After Coregistration','Before Coregistration','Location','NorthWest')
        box on
        ofile=[texttar,'',textref,'pdf'];
        print('-dpng','-r400',ofile) 
        try
        saveas(gcf,ofile,'fig')
        end

        fprintf ('\n height mean std before coregistration: %f %f\n ',out{1,4}.S{1,4}.height_mean,out{1,4}.S{1,4}.height_std)
        fprintf ('\n height mean std after coregistration: %f %f\n ',out{1,4}.S{1,5}.height_mean,out{1,4}.S{1,5}.height_std)

        % out{1,4}.S{1,4}.height_mean
        % out{1,4}.S{1,4}.height_std
        % out{1,4}.S{1,5}.height_mean
        % out{1,4}.S{1,5}.height_std
        
        end

    end %if plot

% close all
    	
    end % if j %idem2
end  %idem1

save test3.mat -v7.3

% removing list of idd2 (bad quality strips) that may be used as the target DEM;
% M=ismember(offsets.i,idd2)|ismember(offsets.j,idd2); %the logical matrix identifying bad pairs; 
%the logical matrix identifying bad pairs;coregisterdems may output NaNs
% M=ismember(offsets.i,idd2)|ismember(offsets.j,idd2)|isnan(offsets.dx)|isnan(offsets.dy)|isnan(offsets.dz);
M=ismember(offsets.i,idd2)|ismember(offsets.j,idd2)|isnan(offsets.dx)|isnan(offsets.dy)|isnan(offsets.dz)|(offsets.dxe==0)|(offsets.dye==0)|(offsets.dze==0);

offsets.i(M)=[]; %dem 1 index
offsets.j(M)=[]; %dem 2 index
offsets.dz(M)=[]; %[n×1 double] z offset (dem 2 - dem 1) -> should be dem1 -dem2; otherwise results are wrong!
offsets.dx(M)=[]; % [n×1 double] x offset (dem 2 - dem 1) -> see testadjustoffsets.m.
offsets.dy(M)=[]; % [n×1 double] y offset (dem 2 - dem 1)
offsets.dze(M)=[]; % [n×1 double] z offset 1-sigma error
offsets.dxe(M)=[]; % [n×1 double] x offset 1-sigma error
offsets.dye(M)=[]; % [n×1 double] y offset 1-sigma error
offsets.mean_dz_coreg(M)=[]; % [n×1 double] mean diff in z after corgestration
offsets.median_dz_coreg(M)=[]; % [n×1 double] median diff in z after corgestration
offsets.sigma_dz_coreg(M)=[]; % [n×1 double] std dev of diff in z after corgestration

%save test3.mat -v7.3

if exist('offsets','var')
%[dZ,dX,dY] = adjustOffsets(offsets);
prctilemax=90;  %90 percentile as the threshold;
offsetErrMax_dx=prctile([offsets.dxe(:)],prctilemax); 
offsetErrMax_dy=prctile([offsets.dye(:)],prctilemax); 
offsetErrMax_dz=prctile([offsets.dze(:)],prctilemax); 
offsetErrMax=max([offsetErrMax_dx,offsetErrMax_dy,offsetErrMax_dz]);
fprintf(['\n offsetErrMax = ',num2str(offsetErrMax),'\n']) 
offsetErrMax=min([2, offsetErrMax]);
fprintf(['\n offsetErrMax = ',num2str(offsetErrMax),'\n']) 

min_sigma_dz_coregMax=prctile([abs(offsets.sigma_dz_coreg(:))],prctilemax); %50 percentile
fprintf(['\n min_sigma_dz_coregMax= ',num2str(min_sigma_dz_coregMax),'\n'])
min_sigma_dz_coregMax=min([10, min_sigma_dz_coregMax]);
fprintf(['\n min_sigma_dz_coregMax= ',num2str(min_sigma_dz_coregMax),'\n'])

min_abs_mean_dz_coregMax=prctile([abs(offsets.mean_dz_coreg(:))],prctilemax); %50 percentile
fprintf(['\n min_abs_mean_dz_coregMax= ',num2str(min_abs_mean_dz_coregMax),'\n']) 
min_abs_mean_dz_coregMax=min([1, min_abs_mean_dz_coregMax]);
fprintf(['\n min_abs_mean_dz_coregMax= ',num2str(min_abs_mean_dz_coregMax),'\n']) 

min_abs_median_dz_coregMax=prctile([abs(offsets.median_dz_coreg(:))],prctilemax); %50 percentile
fprintf(['\n min_abs_median_dz_coregMax= ',num2str(min_abs_median_dz_coregMax),'\n']) 
min_abs_median_dz_coregMax=min([5, min_abs_median_dz_coregMax]);
fprintf(['\n min_abs_median_dz_coregMax= ',num2str(min_abs_median_dz_coregMax),'\n']) 

%[dZ,dX,dY] = adjustOffsets(offsets,'offsetErrMax',offsetErrMax); %0.5 0.1
%[dZ,dX,dY] = adjustOffsets(offsets,'offsetErrMax',offsetErrMax,'min_abs_mean_dz_coregMax',min_abs_mean_dz_coregMax,'min_sigma_dz_coregMax',min_sigma_dz_coregMax,'min_abs_median_dz_coregMax',min_abs_median_dz_coregMax); %0.5 0.1
[dZ,dX,dY] = adjustOffsets(offsets,'offsetErrMax'); 
%  dX4Sg(j,1:3)=p;
dX4Sg1=[dZ,dX,dY];
[m,~]=size(dX4Sg1);
dX4Sg(1:m,:)=dX4Sg1;
else
    warning('\n coreg3.m: no offsets for re-adjustment! \n')
end

% remove strips that are not coregistered.
if 1  %
%remove strips that are not coregistered (idd3) and bad quanlity (idd2).
idd3=find(any(isnan(dX4Sg),2));%find the id of any dZ dX dY that are nan;
idd=unique([idd2(:);idd3(:);]);
idregion(idd)=[];%XYb(idd)=[];dzxy(idd)=[];
dX4Sg(idd,:)=[];%rmsreg2(idd)=[];
end
data0r=[];%SAVE MEMORY.

fprintf ('\n Done with coregistration.\n ')
toc

return
end
