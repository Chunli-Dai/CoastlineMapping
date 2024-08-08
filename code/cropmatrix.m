function [wm]=cropmatrix(data,M,flag)
%crop matrix using the M matrix (1 valid data, 0 non valid data); 
%data.z and M should have the same size;
if nargin==2 %old [wm]=cropmatrix(data,M);
jump=M;
xout=data.x;yout=data.y;
[X,Y]=meshgrid(xout,yout);
xmin=min(X(jump==1));xmax=max(X(jump==1));
ymin=min(Y(jump==1));ymax=max(Y(jump==1));
idx=xout>=xmin&xout<=xmax;idy=yout>=ymin&yout<=ymax;
wm.x=xout(idx);wm.y=yout(idy);wm.z=data.z(idy,idx);
elseif nargin==3
    %e.g. flag=2 %M=rang0;
    rang2=M;
    Mx=data.x>=rang2(1)&data.x<=rang2(2);My=data.y>=rang2(3)&data.y<=rang2(4);

    %crop the area around the target cluster 
    wm.x=data.x(Mx);wm.y=data.y(My);
    wm.z=data.z(My,Mx); % %input mapped area
    
end

return
end