function [p,perr,dzmms,dz,cpts,iter]=readcoregsetsm(in)
%read the coregistration parameters from setsm code.
% in: input structure
	p=[NaN NaN NaN];perr=p;dzmms=p;
	dz=[];cpts=[];
	iter= 1;%initialize

	odircoregi=in.odircoregi;

	if ~isfield(in, 'flag') %DEM coreg results

	idem2=in.idem2;
%odircoregi=in.odircoregi;
	tarimagep=in.tarimagep;
	datatarz=in.datatarz; tardem=in.tardem;refdem=in.refdem;

    %plot the animation of images and control points before and after coregistration
    %refers to /home/dai.56/arcticdemapp/river/riverwork/coregtest1/plotcontrolpts.m
    filecpt=[odircoregi,'/DEM_gcps_1.txt'];

    try

        if ~exist(filecpt,'file')
            fprintf(['\n setsm DEM coregistration failure id: ',num2str([idem2]),' ',tarimagep,' ',filecpt,'\n'])
            iter=49;   
            p=[NaN NaN NaN];
        else
            cpts=load(filecpt); %x y z

            %get coregistration parameters
            %txy=[-1.86 5.04 ];
            coregfile=[odircoregi,'DEM_coreg_result.txt'];
            c=textread(coregfile,'%s','delimiter','\n');   
            [~,name,~] =fileparts([strtrim(tarimagep)]);
            r=find(~cellfun(@isempty,strfind(c,name(1:end-3))));
            %Sigma0[meter]   Tx[meter]       Ty[meter]       Tz[meter]       Tx_std[meter]   Ty_std[meter]   Tz_std[meter]  
            %Dist(mean)      Dist_std(mean)  Dist(med.)      Dist_std(med.)  dh_mean         dh_med.         dh_std         
            %NumberOfCPs     processing time
            if ~isempty(r)
                c2=c{r};
                r1=strfind(c2,'dem');c2([1:r1(end)+2])='';
                [tmp]=sscanf(c2, '%f',[1,16]);

                txy=-[tmp(2), tmp(3)];   
                p=-[tmp(4) txy]; %zxy
                perr=[tmp(7) tmp(5) tmp(6) ];
                %mean median and std of dz over control points;
                dzmms=[tmp(8), tmp(10),tmp(9)]; % (8,9,10): dH_cp(mean)     dH_cp_std(mean) dH_cp(med.)
    %           dzmms=tmp(12:14); %dh_mean         dh_med.         dh_std  
            else
                iter=49;   
                p=[NaN NaN NaN];
            end %r isempty

        %Inf and NaN are not permitted.
            if any(isnan(p)|~isfinite(p))
               fprintf(['\n setsm DEM coregistration failure id: ',num2str([idem2]),' ',tarimagep,' ',filecpt,'\n'])
               iter=49;   
            else
	       %for plot
	       if isempty(datatarz)
	          dz=[];
	       else
               z2out = interp2(tardem.x'-p(2) ,tardem.y-p(3),double(datatarz)-p(1) ,refdem.x',refdem.y,'*linear');
               dz=z2out-refdem.z;
               M=isnan(z2out)|refdem.z==-9999;dz(M)=nan;
	       end
            end
        end  % exist

    catch e

           fprintf('readcoregsetsm.m There was an error! The message was:\n%s',e.message);
           fprintf(['\n setsm DEM coregistration failure id: ',num2str([idem2]),' ',tarimagep,' ',filecpt,'\n'])
           iter=49;
           p=[NaN NaN NaN];

    end %try 
	
	elseif in.flag ==1 %orthoimage 

	try
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
        p=[0 txy]; %zxy

    catch e

           fprintf('readcoregsetsm.m There was an error! The message was:\n%s',e.message);
           fprintf(['\n setsm DEM coregistration failure id: ',num2str([idem2]),' ',tarimagep,' ',filecpt,'\n'])
           iter=49;
           p=[NaN NaN NaN];

    end %try 
        
	end %if

	return
end
