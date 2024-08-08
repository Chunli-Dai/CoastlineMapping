%in case of too many repeats, use only the summer and latest (in time) novmax*3 measurements for grouping.
function [idd,idkp]=selectdate(str,novmax)
%str: date and sensor string, e.g., 'WV02_20150617151345' 'WV02_20150617'
%novmax: the maximum number of strings to be kept.
%idd: output, index of str to be deleted.
% idkp: index of string to be kept and in order.

n=length(str);
idall=1:n;

%Get the date string
month=zeros(n,1);year=zeros(n,1);
for i=1:n
    month(i)=str2double(str{i}(10:11));
    year(i)=str2double(str{i}(6:9));
end

% Sort the string from hottest to coldest month 
% North canada: 07 08 06 09 05 10, 11 04 03 02 12 01 ^_^
% Greenland: 06, 07 08 09 05 11 10 04 03 02 12 01
% https://en.wikipedia.org/wiki/Climate_of_Greenland 
%Month	 Jan	Feb	 Mar	Apr	May	Jun	 Jul  Aug	Sep	Oct	 Nov  Dec	Year
%Celsius 15.3   16.0 16.0 19.1 24.8 30.1 26.3 25.2 24.9 19.3 21.6 15.9  

ids=[];

id=find(month==7);
[~,id2]=sort(year(id),'descend'); %sort year reversely
id=id(id2);
ids=[ids(:);id(:)];

id=find(month==8);
[~,id2]=sort(year(id),'descend'); %sort year reversely
id=id(id2);
ids=[ids(:);id(:)];

id=find(month==6);
[~,id2]=sort(year(id),'descend'); %sort year reversely
id=id(id2);
ids=[ids(:);id(:)];

id=find(month==9);
[~,id2]=sort(year(id),'descend'); %sort year reversely
id=id(id2);
ids=[ids(:);id(:)];

id=find(month==5);
[~,id2]=sort(year(id),'descend'); %sort year reversely
id=id(id2);
ids=[ids(:);id(:)];

id=find(month==10);
[~,id2]=sort(year(id),'descend'); %sort year reversely
id=id(id2);
ids=[ids(:);id(:)];

id=find(month==11);
[~,id2]=sort(year(id),'descend'); %sort year reversely
id=id(id2);
ids=[ids(:);id(:)];

id=find(month==4);
[~,id2]=sort(year(id),'descend'); %sort year reversely
id=id(id2);
ids=[ids(:);id(:)];

id=find(month==3);
[~,id2]=sort(year(id),'descend'); %sort year reversely
id=id(id2);
ids=[ids(:);id(:)];

id=find(month==2);
[~,id2]=sort(year(id),'descend'); %sort year reversely
id=id(id2);
ids=[ids(:);id(:)];

id=find(month==12);
[~,id2]=sort(year(id),'descend'); %sort year reversely
id=id(id2);
ids=[ids(:);id(:)];

id=find(month==1);
[~,id2]=sort(year(id),'descend'); %sort year reversely
id=id(id2);
ids=[ids(:);id(:)];

if length(ids) ~= n;
	fprintf('\n selectdate.m: ids length does not match n.')
end

idd=ids(novmax+1:end); %to be deleted strings;

strkp=str; strkp(idd)=[]; %kept strings not sorted.

if novmax <= length(ids)
idkp=ids(1:novmax);
else %keep all strings
	idkp=ids;
end
strkp2=str(idkp); % kept strings sorted by hot months and latest year.
% fprintf(['\n Kept strings sorted by hot months and latest year:']);
% strkp2{:}

return
end
