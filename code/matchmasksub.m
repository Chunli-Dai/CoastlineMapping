function [z_tideh,iisv]=matchmasksub(ik,iselopc,xid,yid,xout,yout,jump,un,isv,datarsv,tideh)

	
        nsuby=length(yout);nsubx=length(xout);

        isel=iselopc(ik);
        xid1=max([1,xid(ik)-2]):min([xid(ik)+2,nsubx]);yid1=max([1,yid(ik)-2]):min([yid(ik)+2,nsuby]);
	xt=xout(xid1);yt=yout(yid1);    jumpt=jump(yid1,xid1);
    npts=length(jumpt(:));
    
    %search the image that has the same coastline as the final mask
    id=strread(un{isel},'%d');%id of overlapping files
    simi=zeros(length(id),1);
    for j=1:length(id)
        i=id(j);
        iisv=find(isv==i);
        datar=datarsv(iisv);
        if ~isempty(datar.z)
           Mj=interp2(datar.x,datar.y,datar.z,xt,yt','*nearest',-1);
           df=jumpt-Mj;
           simi(j)=sum(sum(df==0))/npts;%ratio of same mask pixels
           if 0
           figure;imagesc(xt*1e-3,yt*1e-3,jumpt-Mj);
           title([num2str(j),'; tide height ',num2str(tideh(iisv)),' m']);colorbar;colormap jet;caxis([-1 1])
           end
        end
    end
    [~,js]=max(simi);
    i=id(js);iisv=find(isv==i);
	z_tideh=tideh(iisv);
	%z(ik)=tideh(iisv);

return
end
