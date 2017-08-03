%% cleaned up version of vorticity_wavenumbers.m

%% load bathy and dimensions
h = load('~/Dropbox/RODSEX/survey/funwaveC_bathy/bathy_RODSEX_0925_09281600_1D_D_dx1p3.depth');




% xi_frf x-coord of model bathy in frf coords
load ~/Dropbox/RODSEX/survey/funwaveC_bathy/bathy_RODSEX_0925_09281600_1D_xifrf_D_dx1p3.mat
% % get rid of the sponge layer and the swash
% % the sponge is 5 gridpoints long, and the swash extends another 10
% xi_frf(570:end) = NaN; 
xi_frf = [xi_frf, xi_frf(end)-dx] + dx/2;


% model params
dx = 1.33333;
dy = 1.33333;
numx = 585;
numy = 1200;

% index for ring location
iring = 510;

%% vorticity: "ring" averaged, and wavenumber spectra

numfiles = 1900;
datadir = '/Volumes/DAVIDCLARK/fC_RODSEX_0928_D3_dx1p3/';
% datadir = '/Volumes/ThunderBay/fC_RODSEX_0928_D3_dx1p3/';
% datadir = '~/Dropbox/RODSEX/funwaveC/';


% make averaging areas
RING = struct;
for ii=1:25;
    RING(ii).bin_length = (ii-1)*2+1; % odd numbers
    RING(ii).i = (iring-(ii-1)):(iring+(ii-1));
    RING(ii).j = 1:RING(ii).bin_length;
    RING(ii).steps = floor(numy/RING(ii).bin_length);  
    RING(ii).vort = zeros(numfiles,RING(ii).steps);
end



% load a file and calucluate one sided spectrum size
load(sprintf('%ssnap_vort_l2_%4.0f.mat',datadir,3599))
% stop = alongshore spectrum length
sd = size(vort,2);
if mod(sd,2)==0               % only half the data is good
    stop=sd/2+1;                % we are making a onesided spectrum
else						  
    stop=(sd+1)/2;
end

% cross-shore indecies for region near ring, where the cross-shore
% wavenumber spectrum will be calculated. The ring location at index 510 is
% roughly the max alongshore vorticity variance (near k=0.027 1/m), and
% decreases by about 50% at +-30 cross-shore grid-points.
ixs = 510-30:510+30; 
% stop2 = cross-shore spectrum length
sd2 = length(ixs);
if mod(sd2,2)==0               % only half the data is good
    stop2=sd2/2+1;                % we are making a onesided spectrum
else						  
    stop2=(sd2+1)/2;
end

%window the crosshore vort before spectrum
w = repmat(hanning(length(ixs)),1,size(vort,2));

% initialize vars
n=0;
vort_ring = zeros(numfiles,60);
wavenum_spec_vort = zeros(size(vort,1),stop);
wavenum_spec_vort_xshore_short = zeros(numy,stop2);
% total_var = 0;

for ii = 3599:(3599+numfiles)
    n=n+1;
    fprintf('%g\n',ii)
    load(sprintf('%ssnap_vort_l2_%4.0f.mat',datadir,ii))
    
    % ring averages
    for jj = 1:length(RING) 
        for kk = 0:(RING(jj).steps-1);
            RING(jj).vort(n,kk+1) = sum(sum(vort(RING(jj).i,RING(jj).j ...
                + (RING(jj).bin_length*kk))))/(RING(jj).bin_length.^2);
        end
    end
    
    % alongshore wavenumber spectra
    [mpsd,wavenums]=mypsd(vort.',dy);
    [mpsd,wavenums]=mywelch(vort.',dy,1,0);
    wavenum_spec_vort = wavenum_spec_vort + mpsd.'/numfiles;
    
    % cross-shore wavenumber spectra, using cross-shore hanning window
    winvort0 = detrend(vort(ixs,:).*w);
    var0 = var(vort(ixs,:));
    var1 = var(winvort0);
    winvort1 = winvort0.*repmat(sqrt(var0./var1),length(ixs),1);
    [mpsd2,wavenums_x_short]=mypsd(winvort1,dx);
    wavenum_spec_vort_xshore_short = wavenum_spec_vort_xshore_short + mpsd2.'/numfiles;
    
    
    
end


%% get vorticity variance over various "ring" averaging areas 

% vorticity var inside "ring", over space (one alongshore transect)
% and time 
varvort = zeros(length(RING),1);
% length of one size of averaging area
rsize = zeros(length(RING),1);
for ii = 1:length(RING)
    varvort(ii) = var(RING(ii).vort(:));
    rsize(ii) = RING(ii).bin_length*dx;
end

figure(3); clf
plot(rsize,varvort)



%% combine cross- and alongshore wavenumber spectra for comparison with
% ring averages




short_for_interp = mean(wavenum_spec_vort_xshore_short);
short_for_interp(1) = short_for_interp(2); % remove mean and extrapolate


wavenum_spec_vort_xshore_short_interp = interp1(wavenums_x_short,short_for_interp,wavenums,'spline','extrap');
wavenum_spec_vort_xshore_short_interp(1) = mean(wavenum_spec_vort_xshore_short(:,1)); % put the mean back in



%% compare ring size vs vorticity variance with wavenumber range vs the
% spectrum integrated over that range (variance)

ring_wavenumber = 1./(rsize);
spec_partial_var = zeros(size(rsize));
wavenum_limit =  zeros(size(rsize));
spec_partial_var_alongshore_only = zeros(size(rsize));
spec_partial_var_xshore_only = zeros(size(rsize));
dk = wavenums(2) - wavenums(1);
for ii = 1:length(rsize)
    [dump,imax] = min(abs(wavenums-ring_wavenumber(ii)));
    wavenum_limit(ii) = wavenums(imax);
    spec_partial_var(ii) = (sum(wavenum_spec_vort_xshore_short_interp(2:imax) + mean(wavenum_spec_vort(RING(ii).i,2:imax)))*dk);
    spec_partial_var_alongshore_only(ii) = (sum(mean(wavenum_spec_vort(RING(ii).i,2:imax),1))*dk);
    spec_partial_var_xshore_only(ii) = (sum(wavenum_spec_vort_xshore_short_interp(2:imax))*dk);
    
end


%% save vars for later


save ~/Dropbox/RODSEX/funwaveC/vorticity_wavenumbers_spectral_sum_comparison_data.mat ...
    RING dx dy h iring ixs numx numy ring_wavenumber rsize spec_partial_var ...
    varvort vort_ring wavenum_limit wavenum_spec_vort ...
    wavenum_spec_vort_xshore_short wavenum_spec_vort_xshore_short_interp ...
    wavenums wavenums_x_short xi_frf spec_partial_var_alongshore_only ...
    spec_partial_var_xshore_only


%%
load ~/Dropbox/RODSEX/funwaveC/vorticity_wavenumbers_spectral_sum_comparison_data.mat

%% figure: ring size vs wavenumber spectrum integral 
figure(9); clf
plot(ring_wavenumber, varvort, ...
     wavenum_limit, spec_partial_var_xshore_only ,wavenum_limit, spec_partial_var_alongshore_only)
xlabel('1/ringSize, wavenum\_limit (1/m)')
ylabel('vort variance (1/s^2)')
legend('ring average','spectrum cross-shore','spectrum alognshore','location','southeast')

%% compare total variance in cross and alongshore wavenumber spectra

datadir = '/Volumes/DAVIDCLARK/fC_RODSEX_0928_D3_dx1p3/';
% datadir = '/Volumes/ThunderBay/fC_RODSEX_0928_D3_dx1p3/';
% datadir = '~/Dropbox/RODSEX/funwaveC/';
load(sprintf('%ssnap_vort_l2_%4.0f.mat',datadir,3599+numfiles))

% variance of last snapshot (all the snaps I checked had total variance
% that was within a couple percent of eachother in the ring ixs area)
v2 = vort(ixs,:);
sprintf('total variance in ixs region = %f',var(v2(:)))

% total variance in alongshore wavenumber spectra
alongshore_spectra_total_var = sum(mean(wavenum_spec_vort(ixs,:)))*(wavenums(2)-wavenums(1))

% total variance in xshore wavenumber spectra
xshore_spectra_total_var = sum(mean(wavenum_spec_vort_xshore_short))*(wavenums_x_short(2)-wavenums_x_short(1))
xshore_spectra_interp_total_var = sum(wavenum_spec_vort_xshore_short_interp)*(wavenums(2)-wavenums(1))
