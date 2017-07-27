% steve data 10/14 - 10/15

% load ~/work/RODSEX/data_rough_20131023/day_vort/RD_20131014.mat
load /Volumes/ThunderBay/RODSEX/data_rough_20131023/day_vort/RD_20131014.mat
RD3 = RD; %use 7-24 hours
% load ~/work/RODSEX/data_rough_20131023/day_vort/RD_20131015.mat
load /Volumes/ThunderBay/RODSEX/data_rough_20131023/day_vort/RD_20131015.mat
RD4 = RD; %use 1-6 hours

startday1014 = 287;  % day of the year

anan = zeros(24,1);
lrec = zeros(24,1);
nn = 0;
vort_out = zeros(24*24576,2);
for ii = 8:24;
    vort_out((ii-8)*24576 + (1:24576),1) = startday1014 + (ii-1)/24 + (0:(24576-1)).'/24/3600/8;  % time in fractional days of the year 
    vort_out((ii-8)*24576 + (1:24576),2) = RD3(ii).vort;
    nn=nn+1;
anan(nn) = length(find(isnan(RD3(ii).vort)));
lrec(nn) = length(RD3(ii).vort);
end

for ii = 1:7;
    vort_out((ii+24-8)*24576 + (1:24576),1) = startday1014 + 1 + (ii-1)/24 + (0:(24576-1)).'/24/3600/8;  % time in fractional days of the year 
    vort_out((ii+24-8)*24576 + (1:24576),2) = RD4(ii).vort;
    nn = nn+1;
anan(nn) = length(find(isnan(RD4(ii).vort)));
lrec(nn) = length(RD4(ii).vort);
end
    

fid  = fopen('steve_vort_1014_1015.txt','w+');
fprintf(fid,'%3.8f   %f\n',vort_out.')
fclose(fid)




