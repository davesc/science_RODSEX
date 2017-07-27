% steve data 09/27

% load('/Volumes/ThunderBay/RODSEX/data_rough_20131023/day_vort/RD_20130927.mat');
load('/Volumes/OWC HD/RODSEX/data_rough_20131023/day_vort/RD_20130927.mat');


RD1 = RD;

startday0927 = 270;  % day of the year

anan = zeros(24,1);
lrec = zeros(24,1);
nn = 0;
vort_out = zeros(24*24576,2);
for ii = 1:24;
    vort_out((ii-1)*24576 + (1:24576),1) = startday0927 + (ii-1)/24 + (0:(24576-1)).'/24/3600/8;  % time in fractional days of the year 
    vort_out((ii-1)*24576 + (1:24576),2) = RD1(ii).vort;
    nn=nn+1;
anan(nn) = length(find(isnan(RD1(ii).vort)));
lrec(nn) = length(RD1(ii).vort);
end


fid  = fopen('steve_vort_0927.txt','w+');
fprintf(fid,'%3.8f   %f\n',vort_out.')
fclose(fid)




