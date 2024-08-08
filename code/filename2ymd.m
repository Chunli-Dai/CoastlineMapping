function ymd=filename2ymd(filename)
%input filename: character string
%output ymd: yyyymmdd,character string

mfile=filename;
[~,filename,~]=fileparts(mfile);

if strcmp(filename(1:6),'SETSM_')
%/fs/byo/howat-data3/data1/chunliwork/landslide/site1Eureka/stripdata/SETSM_WV01_20131005_1020010026874D00_10200100251D3C00_seg1_2m_v3.0_mdf.txt
ids=6;
ymd=filename(ids+(6:13));;
else
%WV01_20131005_1020010026874D00_10200100251D3C00_seg1_2m
ymd=filename(6:13);
end

return
end
