%  function [Usw, R] = hb06_lippmann_ratio(U,V,P, doffu, doffp);
%   - Falk Feddersen
%   - returns Spectrum Sm so that   2*sum(Sm)*df = var(U)
%     as required by Trowbridge appendix for estimating dissipation 
%


function [Usw, R, U_igvar, V_igvar, E_igvar, depth] = vortex_lippmann_ratio(U, V, ETA, depth, dt)




fs = 1/dt;  % Hz



if (all(isnan(U)) | all(isnan(V)) | all(isnan(ETA)) ),
   Usw = nan;
   R = nan;
   U_igvar = nan;   
   V_igvar = nan;   
   E_igvar = nan;   
   depth = nan;
   
   return;
   
   % finish this
end;


% does this work for buried paros sensors??

% load(bathy_file); 
% load mv_snap_eta.mat
% depth = h + mean(B.');  % This is the mean water depth

if (depth <= 0.0),
   disp('** Warning!  depth < 0.  Abort')
   Usw = nan;
   R = nan;
   U_igvar = nan;   
   V_igvar = nan;   
   E_igvar = nan;   
   depth = nan;   
   return;
end;   


%%%%%%%%%%%%%%%%%%%%%%

nfft = 256;

chnks = floor(size(U,1)/nfft)
lap = .5;

size(U)
[Suu,fm] = mywelch(U,dt,chnks,lap); 
[Svv,fm] = mywelch(V,dt,chnks,lap); 
[See,fm] = mywelch(ETA,dt,chnks,lap); 


% The result is now that  sum(Suu)*df = var(U) - actually < var(U) because of detrending


% DEPTH CORRECTION -------------------------------------------------------------
% find correction/conversion rcoefs at each f-band 
% to calc sea surface elevation and convert velocities to pressure units 


% SPECTRAL WEIGHTED AVERAGES & STATS -----------------------------------------
% find indices of freq bands


g = 9.81;    % gravity, cm/s^2
fcutoff = 0.3 ;  % cutoff frequency 
omega = 2*pi*fm;   % This is radian frequency

warning off MATLAB:divideByZero
k = wavenum( omega , depth);
warning on MATLAB:divideByZero



i_ig = find(fm>0.004 & fm<0.03);
df = abs(fm(2) - fm(1));



U_igvar = sum(Suu(i_ig,:))*df;
V_igvar = sum(Svv(i_ig,:))*df;
E_igvar = sum(See(i_ig,:))*df;

Utot_igvar = U_igvar + V_igvar;

size((Utot_igvar./E_igvar))
size(g./depth)
R = (Utot_igvar./E_igvar)./(g./depth);  % as in Lippmann et al. 99 eq 16

Usw = sqrt( Utot_igvar .* (1-1./R) );

%keyboard






function k = wavenum(om,h);
% function that takes the radian wave frequency and
% a vector of depths and returns the wavenumber at that
% depth by solving the dispersion relationship
% using newtons method

% the initial guess will be the shallow water wavenumber
% frequency runs down (om needs to be column) and cross-shore depth runs
% across (h needs to be a row)
%%
%h = h.';
om = om.';

som=size(om);
sh=size(h);

om=repmat(om,[1 sh(2)]);
h=repmat(h,[som(1) 1]);

g = 9.81;

%size(om)
%size(sqrt(g*h))
k = om./sqrt(g*h);


f = g*k.*tanh(k.*h) - om.^2;
%%
for a=1:10 
  dfdk = g*k.*h.*(sech(k.*h)).^2 + g*tanh(k.*h);
  k = k - f./dfdk;
  f = g*k.*tanh(k.*h) - om.^2;
end



