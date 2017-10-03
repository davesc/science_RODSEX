%% cleaned up version of vorticity_wavenumbers.m

%% load bathy and dimensions
h = load('~/Dropbox/RODSEX/survey/funwaveC_bathy/bathy_RODSEX_0925_09281600_1D_D_dx1p3.depth');

% model params
dx = 1.33333;
dy = 1.33333;
numx = 585;
numy = 1200;

% index for ring location
iring = 510;


% xi_frf x-coord of model bathy in frf coords
load ~/Dropbox/RODSEX/survey/funwaveC_bathy/bathy_RODSEX_0925_09281600_1D_xifrf_D_dx1p3.mat
% % get rid of the sponge layer and the swash
% % the sponge is 5 gridpoints long, and the swash extends another 10
xi_frf(570:end) = NaN; 
xi_frf = [xi_frf, xi_frf(end)-dx] + dx/2;




%% vorticity: "ring" averaged, and wavenumber spectra

numfiles = 1900;
% datadir = '/Volumes/DAVIDCLARK/fC_RODSEX_0928_D3_dx1p3/';
datadir = '/Volumes/ThunderBay/fC_RODSEX_0928_D3_dx1p3/';
% datadir = '~/Dropbox/RODSEX/funwaveC/';


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



% make averaging areas
RING = struct;
for ii=1:25;
    RING(ii).bin_length = (ii-1)*2+1; % odd numbers
    RING(ii).i = (iring-(ii-1)):(iring+(ii-1));
    RING(ii).j = 1:RING(ii).bin_length;
    RING(ii).jsteps = floor(numy/RING(ii).bin_length);
    tmp = (length(ixs) - RING(ii).bin_length)/2;
    RING(ii).isteps = -tmp:1:tmp;
%     RING(ii).vort = zeros(numfiles,RING(ii).steps);
    RING(ii).meanvort2 = 0;
    RING(ii).meanvort = 0;
    RING(ii).meanvort_all = zeros(numfiles,1);
    RING(ii).varvort = 0;
    RING(ii).varvort_valid = 0;
%     RING(ii).avg_filter = ones(RING(ii).bin_length) ...
%         / ( (length(ixs) - RING(ii).bin_length + 1) * (numy - RING(ii).bin_length + 1) * numfiles); 
    RING(ii).avg_area = false(length(ixs),numy);
    RING(ii).avg_area(ceil(RING(ii).bin_length/2):(length(ixs)-floor(RING(ii).bin_length/2)), ceil(RING(ii).bin_length/2):(numy-floor(RING(ii).bin_length/2))) = true; 
    RING(ii).avg_area_num = sum(RING(ii).avg_area(:));
%     RING(ii).filter_weight = 1 / ( (length(ixs) - RING(ii).bin_length + 1) * (numy - RING(ii).bin_length + 1) * numfiles * RING(ii).bin_length.^2);
%     RING(ii).filter_weight = 1 / (length(ixs) * numy * numfiles * RING(ii).bin_length.^2);
end

% initialize vars
n=0;
vort_ring = zeros(numfiles,60);
wavenum_spec_vort = zeros(size(vort,1),stop);
wavenum_spec_vort_xshore_short = zeros(numy,stop2);
wavenum_spec_vort_mw = zeros(size(vort,1),stop);
wavenum_spec_vort_xshore_short_mw = zeros(numy,stop2);
% total_var = 0;
sz = length(ixs) * numy;

for ii = 3599:(3599+numfiles)
    n=n+1;
    fprintf('%g\n',ii)
    load(sprintf('%ssnap_vort_l2_%4.0f.mat',datadir,ii))
    vort1 = vort(ixs,:);
    vort2 = vort1.^2;
    
    % ring averages
    for jj = 1:length(RING) 
        
%         % convolve averaging filter and add to total
%         RING(jj).meanvort2 = RING(jj).meanvort2 + sum(sum(conv2(vort2(ixs,:),RING(jj).avg_filter,'valid')));
%         RING(jj).meanvort = RING(jj).meanvort + sum(sum(conv2(vort(ixs,:),RING(jj).avg_filter,'valid')));

%         % use image box filter to get means
%         tmp = imboxfilt(vort1, RING(jj).bin_length, 'NormalizationFactor',RING(jj).filter_weight);
%         RING(jj).meanvort = RING(jj).meanvort + sum(sum(tmp.*RING(jj).avg_area));
% 
%         tmp = imboxfilt(vort2, RING(jj).bin_length, 'NormalizationFactor',RING(jj).filter_weight);
%         RING(jj).meanvort2 = RING(jj).meanvort2 + sum(sum(tmp.*RING(jj).avg_area));
%         
        % use image box filter to get means 
        tmp1 = imboxfilt(vort1, RING(jj).bin_length,'padding','symmetric');
        RING(jj).varvort = RING(jj).varvort + (sum(sum((tmp1 - sum(sum(tmp1))/sz ).^2))/sz  / numfiles);
        RING(jj).meanvort_all(ii) = sum(sum(tmp1))/sz;
        
        % this uses the valid portion of the convolution in the center of
        % the sampling area (both methods yield similar results)
        RING(jj).varvort_valid = RING(jj).varvort_valid + (sum(sum((tmp1(RING(jj).avg_area) - sum(sum(tmp1(RING(jj).avg_area)))/RING(jj).avg_area_num ).^2))/RING(jj).avg_area_num  / numfiles);

        
%         
%         RING(jj).meanvort2 = RING(jj).meanvort2 + sum(sum(tmp1.^2))/(length(ixs) * numy * numfiles);

        % I think this is wrong: I should be taking the square after the
        % ring area average
%         tmp2 = imboxfilt(vort2, RING(jj).bin_length,'padding','symmetric');
%         RING(jj).meanvort2 = RING(jj).meanvort2 + sum(sum(tmp2))/(length(ixs) * numy * numfiles);

        
% %         for kk = 0:(RING(jj).jsteps-1);
% %             RING(jj).vort(n,kk+1) = sum(sum(vort(RING(jj).i,RING(jj).j ...
% %                 + (RING(jj).bin_length*kk))))/(RING(jj).bin_length.^2);
% %         end
%         for kk = 0:(RING(jj).jsteps-1);
%             for ll = 1:length(RING(jj).isteps);
%                 RING(jj).meanvort2 = RING(jj).meanvort2 ...
%                                 + sum(sum(vort(RING(jj).i + RING(jj).isteps(ll), RING(jj).j + (RING(jj).bin_length*kk)).^2))...
%                                 /(RING(jj).bin_length.^2 * numfiles * length(RING(jj).isteps) * length(RING(jj).jsteps));
%                 RING(jj).meanvort = RING(jj).meanvort ...
%                                 + sum(sum(vort(RING(jj).i + RING(jj).isteps(ll), RING(jj).j + (RING(jj).bin_length*kk))))...
%                                 /(RING(jj).bin_length.^2 * numfiles * length(RING(jj).isteps) * length(RING(jj).jsteps));
%             end
%         end
    end
    
    % TODO: averaging spectra is missing some variance. Need to include the
    % variation in the mean, 1st bin, of each spectra
    
    % alongshore wavenumber spectra
    [mpsd,wavenums]=mypsd(vort.',dy);
    % mypsd and mywelch are giving the same results 
    [mpsd_mw,wavenums_mw]=mywelch(vort.',dy,1,0);
    wavenum_spec_vort = wavenum_spec_vort + mpsd.'/numfiles;
    wavenum_spec_vort_mw = wavenum_spec_vort_mw + mpsd_mw.'/numfiles;
    
    
    % cross-shore wavenumber spectra, using cross-shore hanning window
    var0 = var(vort(ixs,:));
    winvort0 = detrend(vort(ixs,:)).*w;
    var1 = var(winvort0);
    winvort1 = winvort0.*repmat(sqrt(var0./var1),length(ixs),1);
    [mpsd2,wavenums_x_short]=mypsd(winvort1,dx);
    % mypsd and mywelch are giving the same results
    [mpsd2_mw,wavenums_x_short_mw]=mywelch(vort(ixs,:),dx,1,0);
    wavenum_spec_vort_xshore_short = wavenum_spec_vort_xshore_short + mpsd2.'/numfiles;
    wavenum_spec_vort_xshore_short_mw = wavenum_spec_vort_xshore_short_mw + mpsd2_mw.'/numfiles;
    
    
    
end


%% get vorticity variance over various "ring" averaging areas 

% % vorticity var inside "ring", over space (one alongshore transect)
% % and time 
% varvort = zeros(length(RING),1);
% % length of one size of averaging area
% rsize = zeros(length(RING),1);
% for ii = 1:length(RING)
%     varvort(ii) = var(RING(ii).vort(:));
%     rsize(ii) = RING(ii).bin_length*dx;
% end

% vorticity var inside "ring", over space (one alongshore transect)
% and time 
varvort = zeros(length(RING),1);
varvort_valid = zeros(length(RING),1);
meanvort = zeros(length(RING),1);
meanvort2 = zeros(length(RING),1);
% length of one size of averaging area
rsize = zeros(length(RING),1);
for ii = 1:length(RING)
%     RING(ii).varvort = RING(ii).meanvort2 - RING(ii).meanvort;
    varvort(ii) = RING(ii).varvort;
    varvort_valid(ii) = RING(ii).varvort_valid;
    var_meanvort_all(ii) = var(RING(ii).meanvort_all);
%     meanvort(ii) = RING(ii).meanvort;
%     meanvort2(ii) = RING(ii).meanvort2;
    rsize(ii) = RING(ii).bin_length*dx;
end

figure(3); clf
plot(rsize,varvort, rsize, varvort_valid,'--','linewidth',1.5)
legend('convolution with mirrored edges','convolution with central valid region')
xlabel('ring size')
ylabel('vorticity variance (1/s^2)')

%% combine cross- and alongshore wavenumber spectra for comparison with
% ring averages




short_for_interp = mean(wavenum_spec_vort_xshore_short);
short_for_interp(1) = short_for_interp(2); % remove mean and extrapolate


wavenum_spec_vort_xshore_short_interp = interp1(wavenums_x_short,short_for_interp,wavenums,'spline','extrap');
wavenum_spec_vort_xshore_short_interp(1) = mean(wavenum_spec_vort_xshore_short(:,1)); % put the mean back in



%% compare ring size vs vorticity variance with wavenumber range vs the
% spectrum integrated over that range (variance)

ring_wavenumber = 1./(rsize);
% spec_partial_var = zeros(size(rsize));
wavenum_limit =  zeros(size(rsize));
spec_partial_var_alongshore_only = zeros(size(rsize));
spec_partial_var_xshore_only = zeros(size(rsize));
dk = wavenums(2) - wavenums(1);
for ii = 1:length(rsize)
    [dump,imax] = min(abs(wavenums-ring_wavenumber(ii)));
    wavenum_limit(ii) = wavenums(imax);
%     spec_partial_var(ii) = (sum(wavenum_spec_vort_xshore_short_interp(2:imax) + mean(wavenum_spec_vort(RING(ii).i,2:imax)))*dk);
%     spec_partial_var_alongshore_only(ii) = (sum(mean(wavenum_spec_vort(RING(ii).i,2:imax),1))*dk);
    spec_partial_var_alongshore_only(ii) = (sum(mean(wavenum_spec_vort(ixs,2:imax),1))*dk);
    spec_partial_var_xshore_only(ii) = (sum(wavenum_spec_vort_xshore_short_interp(2:imax))*dk);
    
end


%% save vars for later


save ~/Dropbox/RODSEX/funwaveC/vorticity_wavenumbers_spectral_sum_comparison_data.mat ...
    RING dx dy h iring ixs numx numy ring_wavenumber rsize ...
    varvort vort_ring wavenum_limit wavenum_spec_vort ...
    wavenum_spec_vort_xshore_short wavenum_spec_vort_xshore_short_interp ...
    wavenums wavenums_x_short xi_frf spec_partial_var_alongshore_only ...
    spec_partial_var_xshore_only %spec_partial_var


%%
load ~/Dropbox/RODSEX/funwaveC/vorticity_wavenumbers_spectral_sum_comparison_data.mat

%% figure: ring size vs wavenumber spectrum integral 
figure(9); clf
plot(ring_wavenumber, varvort, ...
     wavenum_limit, spec_partial_var_xshore_only ,wavenum_limit, spec_partial_var_alongshore_only)
xlabel('1/ringSize, wavenum\_limit (1/m)')
ylabel('vort variance (1/s^2)')
legend('ring average','spectrum cross-shore','spectrum alognshore','location','southeast')


%% figure: plot mypsd vs mywelch wavenumber spectra

figure(10); clf
plot(wavenums,mean(wavenum_spec_vort(ixs,:)),wavenums_mw,mean(wavenum_spec_vort_mw(ixs,:)),'--',wavenums_x_short,mean(wavenum_spec_vort_xshore_short), wavenums_x_short_mw,mean(wavenum_spec_vort_xshore_short_mw),'--')
legend('along mypsd','along mywelch','cross mypsd','cross mywelch')
ylabel('vort spectra')
xlabel('wavenumber')
%% compare total variance in cross and alongshore wavenumber spectra

%datadir = '/Volumes/DAVIDCLARK/fC_RODSEX_0928_D3_dx1p3/';
 datadir = '/Volumes/ThunderBay/fC_RODSEX_0928_D3_dx1p3/';
% datadir = '~/Dropbox/RODSEX/funwaveC/';
load(sprintf('%ssnap_vort_l2_%4.0f.mat',datadir,3599+numfiles))

% variance of last snapshot (all the snaps I checked had total variance
% that was within a couple percent of eachother in the ring ixs area)
v2 = vort(ixs,:);
sprintf('total variance in ixs region = %f',var(v2(:)))

% total variance in alongshore wavenumber spectra
alongshore_spectra_total_var = sum(mean(wavenum_spec_vort(ixs,2:end)))*(wavenums(2)-wavenums(1))

% total variance in xshore wavenumber spectra
xshore_spectra_total_var = sum(mean(wavenum_spec_vort_xshore_short(:,2:end)))*(wavenums_x_short(2)-wavenums_x_short(1))
xshore_spectra_interp_total_var = sum(wavenum_spec_vort_xshore_short_interp(:,2:end))*(wavenums(2)-wavenums(1))
