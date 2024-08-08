function [co]=testgif(data1r,data2r,range0,cptsin,refimagep,tarimagep,flagstep)
% in image coregistratino with setsm, data1r is reference image, data2r is target image.
% https://www.mathworks.com/matlabcentral/answers/94495-how-can-i-create-animated-gif-images-in-matlab
%flagstep: 1 before coregistration; 2 after coregistration
co=[];

cpts=[0 0 0 0];
if (nargin == 2) 
  %common range
    rangtar=[min(data2r.x) max(data2r.x) min(data2r.y) max(data2r.y)];  
    rangref=[min(data1r.x) max(data1r.x) min(data1r.y) max(data1r.y)];   
    rangeov=[max(rangtar(1),rangref(1)),min(rangtar(2),rangref(2)), max(rangtar(3),rangref(3)),min(rangtar(4),rangref(4))];
%   reduce resolution
    range0=rangeov;
elseif (nargin == 4)
    %use range0
    %use cpts
	cpts=cptsin;
elseif (nargin >4)
        cpts=cptsin;
[~,filename1,ext1]=fileparts(refimagep);text1=filename1(1:13);
[~,filename2,ext2]=fileparts(tarimagep);text2=filename2(1:13);
end

[cptn, cptm]=size(cpts);

h = figure;
axis equal manual % this ensures that getframe() returns a consistent size
set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 6 4]);
box on
%filename = 'imagecoreg.gif';
if flagstep==1
filename =strrep(tarimagep,'_prep.tif','_prepPre.gif');
elseif flagstep==2
filename =strrep(tarimagep,'_prep.tif','_prepPost.gif');
end

for n = 1:2
    if n==1
        imagesc(data1r.x*1e-3,data1r.y*1e-3,data1r.z);colorbar;title('ref');axis equal;axis(range0*1e-3);view(0,-90)
        if cptm>=2&cptn>=1;
	hold on;plot(cpts(:,1)*1e-3,cpts(:,2)*1e-3,'g>')
        else
          warning(['control points empty from setsm: ',tarimagep])
	end
	M=data1r.z(:)~=0;mean1=nanmean(data1r.z(M));	std1=nanstd(data1r.z(M));
	fprintf(['\n n=1,[mean1-3*std1, mean1+3*std1]:',num2str([mean1-3*std1, mean1+3*std1])])
	try
	caxis([mean1-3*std1, mean1+3*std1]);
	end
	if (nargin >4)
	  title(['Reference image ', text1]);
	end
  elseif n==2
        imagesc(data2r.x*1e-3,data2r.y*1e-3,data2r.z);colorbar;title('tar');axis equal;axis(range0*1e-3);view(0,-90);
        if cptm>=4&cptn>=1;
	hold on;plot(cpts(:,3)*1e-3,cpts(:,4)*1e-3,'r>')
        else
          warning(['control points empty from setsm: ',tarimagep])
	end
	M=data2r.z(:)~=0;mean1=nanmean(data2r.z(M)); std1=nanstd(data2r.z(M));
	fprintf(['\n n=2,[mean1-3*std1, mean1+3*std1]:',num2str([mean1-3*std1, mean1+3*std1])])
	try
	caxis([mean1-3*std1,mean1+3*std1]);
	end
	if (nargin >4)
	  title(['Target image ', text2]);
	end
    end
    
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end
return
end
