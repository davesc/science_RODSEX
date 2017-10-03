
load /Volumes/Data2/RODSEX/data_rough_20131023/day_UV//UV_20130929.mat
rl = load('~/Dropbox/RODSEX/ROD_survey/ROD_fit_location_rot.txt');

figure(2); clf
for ii = 40
 U = sqrt(UV(1).v(ii,:).^2 +UV(1).u(ii,:).^2);
 scatter(rl(:,4),rl(:,5),100,U,'filled')
 pause(0.2)
end

% ring inst coords
x = rl(:,4);
y = rl(:,5);
% center of the ring
cx= mean(x);
cy = mean(y);

% angle of each intrument from shore normal
th = atan2d((x-cx),(y-cy));
I = th<0;
th(I) = th(I)+360;

% interpolate the vlocity magnitudes
Uex = [U(end-1:end), U, U(1:2)];
xUex = -1:16;
xi = 1:.1:15;
Uexi = interp1(xUex,Uex,xi,'spline');

% interpolate the angles 
thi = interp1(1:15,[th;th(1)],xi);

figure(10);clf
plot(thi,Uexi,'o')


