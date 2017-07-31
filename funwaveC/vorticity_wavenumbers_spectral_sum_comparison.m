%% cleaned up version of vorticity_wavenumbers.m

%% load bathy and dimensions
h = load('~/Dropbox/RODSEX/survey/funwaveC_bathy/bathy_RODSEX_0925_09281600_1D_D_dx1p3.depth');

% xi_frf x-coord of model bathy in frf coords
load ~/Dropbox/RODSEX/survey/funwaveC_bathy/bathy_RODSEX_0925_09281600_1D_xifrf_D_dx1p3.mat

% get rid of the sponge layer and the swash
% the sponge is 5 gridpoints long, and the swash seems to extend another 10
xi_frf(570:end) = NaN; 

dx = 1.33333;

dy = 1.33333;

xi_frf = [xi_frf, xi_frf(end)-dx] + dx/2;


numx = 585;
numy = 1200;


%% ring averaged vorticity 

numfiles = 1900;
% make averaging areas
RING = struct;
for ii=1:25;
    RING(ii).bin_length = (ii-1)*2+1; % odd numbers
    RING(ii).i = (iring-(ii-1)):(iring+(ii-1));
    RING(ii).j = 1:RING(ii).bin_length;
    RING(ii).steps = floor(numy/RING(ii).bin_length);  
    RING(ii).vort = zeros(numfiles,RING(ii).steps);
end

datadir = '/Volumes/ThunderBay/fC_RODSEX_0928_D3_dx1p3/';
% datadir = '~/Dropbox/RODSEX/funwaveC/';


n=0;
for ii = 3599:(3599+numfiles)
    n=n+1;
    fprintf('%g\n',ii)
    load(sprintf('%ssnap_vort_l2_%4.0f.mat',datadir,ii))

    for jj = 1:length(RING) 
        for kk = 0:(RING(jj).steps-1);
            RING(jj).vort(n,kk+1) = sum(sum(vort(RING(jj).i,RING(jj).j ...
                + (RING(jj).bin_length*kk))))/(RING(jj).bin_length.^2);
        end
    end
 
end


%% get vorticity std over various "ring" averaging areas 

% vorticity std inside "ring", over space (one alongshore transect)
% and time 
stdvort = zeros(length(RING),1);
% length of one size of averaging area
rsize = zeros(length(RING),1);
for ii = 1:length(RING)
    stdvort(ii) = std(RING(ii).vort(:));
    rsize(ii) = RING(ii).bin_length*dx;
end

figure(3); clf
plot(rsize,stdvort)

%% tmp


datadir = '/Volumes/ThunderBay/fC_RODSEX_0928_D3_dx1p3/';
% datadir = '~/Dropbox/RODSEX/funwaveC/';

% load a file and calucluate one sided spectrum size
load(sprintf('%ssnap_vort_l2_%4.0f.mat',datadir,3599))
sd = size(vort,2);
if mod(sd,2)==0               % only half the data is good
    stop=sd/2+1;                % we are making a onesided spectrum
else						  
    stop=(sd+1)/2;
end

% cross-shore indecies for region near ring
ixs = 510-20:510+20; 

sd2 = length(ixs);
if mod(sd2,2)==0               % only half the data is good
    stop2=sd2/2+1;                % we are making a onesided spectrum
else						  
    stop2=(sd2+1)/2;
end



numfiles = 1900;
vort_ring = zeros(numfiles,60);
wavenum_spec_vort = zeros(size(vort,1),stop);
wavenum_spec_vort_xshore_short = zeros(numy,stop);
n=0;
for ii = 3599:(3599+numfiles)
    n=n+1;
    fprintf('%g\n',ii)
    load(sprintf('%ssnap_vort_l2_%4.0f.mat',datadir,ii))    
    
    [mpsd,wavenums]=mypsd(vort(:,:).',1/dy);
    wavenum_spec_vort = wavenum_spec_vort + mpsd.'/numfiles;
    
    [mpsd2,wavenums_x_short]=mypsd(vort(ixs,:),1/dx);
    wavenum_spec_vort_xshore_short = wavenum_spec_vort_xshore_short + mpsd2.'/numfiles;
    
end




%% ??????  
% check all vars here



short_for_interp = mean(wavenum_spec_vort_xshore_short);
short_for_interp(1) = short_for_interp(2); % remove mean and extrapolate


wavenum_spec_vort_xshore_short_interp = interp1(wavenums_x_short,short_for_interp,wavenums,'linear','extrap');
wavenum_spec_vort_xshore_short_interp(1) = mean(wavenum_spec_vort_xshore_short(:,1)); % put the mean back in


wavenum_spec_vort_short_combined = wavenum_spec_vort_xshore_short_interp + mean(wavenum_spec_vort(ixs,:));



%% compare ring size vs vorticity variance with wavenumber range vs the
% spectrum integrated over that range (variance)

spec_partial_std = zeros(size(rsize));
spec_partial_var_wavenum_limit = 1./(rsize);
wavenum_limit =  zeros(size(rsize));
dk = wavenums(2) - wavenums(1);
for ii = 1:length(rsize)
    [dump,imax] = min(abs(wavenums-spec_partial_var_wavenum_limit(ii)));
    wavenum_limit(ii) = wavenums(imax);
    spec_partial_std(ii) = ...
                   sqrt(sum(wavenum_spec_vort_short_combined(2:imax)*dk));
end

figure(9); clf
plot(spec_partial_var_wavenum_limit, stdvort, ...
     wavenum_limit, spec_partial_std)
xlabel('1/ringSize, wavenum\_limit (1/m)')
ylabel('vort std')


