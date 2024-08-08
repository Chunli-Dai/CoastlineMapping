function [txy,flag]=coregmonostrip(refimage,tarimage,odircoregi,tag);
%Direct coregister of a mono image to a strip image.
% Part of 
%input:
% output: txy: horizontal offsets;
%	  flag: coregistration success 1; failure 0;
	constant % get flagplot
	flag=1; %output flag: 1 is sucess; 0 is failure
	flagplot2=0;%manual control of flagplot; 1 plot; 0 not plot

%       odircoregi=[deblank(odircoreg),'/s',num2str(is),'/'];
        if ~exist(odircoregi,'dir')
          mkdir(odircoregi)
        end
        [refimagep,tarimagep,dataref,datatar]=prepareMJ(refimage,tarimage,odircoregi,tag);

        str=['time setsm -Coreg 1 -image ',refimagep,' -image ', tarimagep, ' -outpath ', odircoregi];
        fprintf([str,'\n'])
        [status, cmdout]=system(str);

        %plot the animation of images and control points before and after coregistration
        %refers to /home/dai.56/arcticdemapp/river/riverwork/coregtest1/plotcontrolpts.m
        filecpt=[odircoregi,'/txt/GCPs_Image_ID_1_level_0.txt'];
        if ~exist(filecpt,'file')|| (isempty(dataref.z) || isempty(datatar.z))
           fprintf(['\n Mono image coregistration failure : ',num2str([]),' ',tarimage,' ',filecpt,'\n'])
	   flag=0;txy=0;
        else

	
	if flagplot==1 ||flagplot2==1
        cpts=load(filecpt);
        %Before Coregistration
    rangtar=[min(datatar.x) max(datatar.x) min(datatar.y) max(datatar.y)];
    rangref=[min(dataref.x) max(dataref.x) min(dataref.y) max(dataref.y)];
    rangeov=[max(rangtar(1),rangref(1)),min(rangtar(2),rangref(2)), max(rangtar(3),rangref(3)),min(rangtar(4),rangref(4))];
    range0=rangeov;
%range0=[-2231000 -2230800 549200 549400]; %%sag river %hi
        meanx=mean(rangeov(1:2));meany=mean(rangeov(3:4));range0=[meanx-100,meanx+100,meany-100,meany+100];

        [co]=testgif(dataref,datatar,range0,cpts,refimagep,tarimagep,1);
	end %plot
        
       %get coregistration parameters
        %txy=[-1.86 5.04 ];
        coregfile=[odircoregi,'coreg_result.txt'];
        c=textread(coregfile,'%s','delimiter','\n');    
        r=find(~cellfun(@isempty,strfind(c,tarimagep)));
        %Ty[meter]       Tx[meter]       avg_roh(average correlation)
        c2=c{r};
        r1=strfind(c2,'tif');c2([1:r1(1)+2])='';
        [tmp]=sscanf(c2, '%f',[1,5]);

        txy=[tmp(4), tmp(3)]; rho=tmp(5);
        %txy.txy=[tmp(4), tmp(3)]; txy.rho=tmp(5);
        
	if flagplot==1||flagplot2==1
        datatarc=datatar;datatarc.x=datatar.x-txy(1);datatarc.y=datatar.y-txy(2);
        [cptn, cptm]=size(cpts);
        if cptm>=4&cptn>=1;
        cpts(:,3)=cpts(:,3)-txy(1); cpts(:,4)=cpts(:,4)-txy(2);
        end
        [co]=testgif(dataref,datatarc,range0,cpts,refimagep,tarimagep,2);
	end %plot

        system(['rm ',refimagep, ' ',tarimagep])
        system(['rm ',odircoregi, '/tmp/*'])
        system(['rm ',odircoregi, '/*.envi*'])
        
        %ps=[0 txy]; %zxy

        end  % exist

	return
end
