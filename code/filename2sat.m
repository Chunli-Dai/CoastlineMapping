function sat=filename2sat(filename)
%input filename: character string
%output ymd: XXXX,character string, e.g. WV01

mfile=filename;
[~,filename,~]=fileparts(mfile);

if strcmp(filename(1:6),'SETSM_')
%/fs/byo/howat-data3/data1/chunliwork/landslide/site1Eureka/stripdata/SETSM_WV01_20131005_1020010026874D00_10200100251D3C00_seg1_2m_v3.0_mdf.txt
ids=6;
sat=filename(ids+(1:4));
else
%WV01_20131005_1020010026874D00_10200100251D3C00_seg1_2m
sat=filename(1:4);
end

return
end
