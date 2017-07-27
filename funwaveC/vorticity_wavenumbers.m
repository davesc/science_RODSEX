%% get vorticity alongshore-wavenumber spectra
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

%%

% cross-shore indicies to calculate wavenums
i_index = 200:569; 
xspec = i_index*dx - 569*dx;

% area to average for "ring" vorticity
iring = 510;
avg_i = (iring-2):(iring+2); % area to average around ring x-loc
xring = iring*dx - 569*dx;
avg_j = 1:5;  % y area to average
ring_mask = true(5);
% round the corners
ring_mask([1,5,21,25])=false;  
%%

datadir = '/Volumes/ThunderBay/fC_RODSEX_0928_D3_dx1p3/';
% datadir = '~/Dropbox/RODSEX/funwaveC/';

load(sprintf('%ssnap_vort_l2_%4.0f.mat',datadir,3599))
sd = size(vort,2);
if mod(sd,2)==0               % only half the data is good
    stop=sd/2+1;                % we are making a onesided spectrum
else						  % modified on 11/11/03
    stop=(sd+1)/2;
end

numfiles = 1900;
vort_ring = zeros(numfiles,60);
wavenum_spec_vort = zeros(size(vort,1),stop);
n=0;
for ii = 3599:(3599+numfiles)
    n=n+1;
    fprintf('%g\n',ii)
    load(sprintf('%ssnap_vort_l2_%4.0f.mat',datadir,ii))
    mvort_tmp = 0;
    for jj = 0:59;
        %         [AJ,AI] = meshgrid(avg_j+(5*jj), avg_i);
%         vort_ring(n,jj+1) = sum(sum(ring_mask.*vort(avg_i,avg_j + (5*jj))))/sum(ring_mask(:));
        vort_ring(n,jj+1) = sum(sum(vort(avg_i,avg_j + (5*jj))))/25;
    end
    
    
    [mpsd,wavenums]=mypsd(vort(:,:).',1/dy);
    wavenum_spec_vort = wavenum_spec_vort + mpsd.'/numfiles;

    
end
%%

save vorticity_wavenumbers_out.mat wavenums wavenum_spec_vort vort_ring avg_i i_index xspec

%%

load ~/Dropbox/RODSEX/funwaveC/vorticity_wavenumbers_out.mat

I = find(xi_frf>100);
iindex = 1:length(xi_frf);

figure(1); clf; 
% pcolor(wavenums(2:end),xi_frf,log10(wavenum_spec_vort(:,2:end))); 
% pcolor(wavenums(2:end),iindex(I),(wavenum_spec_vort(I,2:end))); 
pcolor(wavenums(2:end),iindex(:),log10(wavenum_spec_vort(:,2:end))); 
% pcolor(wavenums(2:end),xi_frf(I),(wavenum_spec_vort(I,2:end))); 
shading flat
xlabel('wavenumber (1/m)')
ylabel('x-shore location (m), 0=shoreline')
title('log10 vorticity spectra')


figure(2); 
plot(wavenums,wavenum_spec_vort(517,:)) % index 517 is near the cross-shore peak in vorticity
xlabel('wavenumbers (1/m)')
ylabel('log10 vorticity spectra')

%% compare ring vorticty variance with summed vorticity spectra


index_ring_wavenums = find(wavenums<(1/5));

ring_vort_var = var(vort_ring(:));

dk = wavenums(2)-wavenums(1);
spec_vort_var = sum(sum(wavenum_spec_vort(avg_i,index_ring_wavenums)))*dk/length(avg_i);
spec_vort_var2 = sum(sum(wavenum_spec_vort(avg_i,:)))*dk/length(avg_i);


%% compare output from different averaging areas

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
            RING(jj).vort(n,kk+1) = sum(sum(vort(RING(jj).i,RING(jj).j + (RING(jj).bin_length*kk))))/(RING(jj).bin_length.^2);
        end
    end
 
end



%%
stdvort = zeros(length(RING),1);
rsize = zeros(length(RING),1);
for ii = 1:length(RING)
    stdvort(ii) = std(RING(ii).vort(:));
    rsize(ii) = RING(ii).bin_length*dx;
end

figure(3); clf
plot(rsize,stdvort)


%%
load vorticity_wavenumbers_out.mat

save vorticity_wavenumbers_out.mat wavenums wavenum_spec_vort vort_ring avg_i i_index xspec RING stdvort rsize
%%


%% cross-shore vorticity wavenumber spectra
% The vorticity changes a lot in the cross-shore 
% so I need to be careful about taking the spectra over a small enough
% cross-shore region that I don't include the cross-shore structure when
% comparing to the ring 


% first make x-shore wavenum spectra using the whole valid cross-shore
% region of the model

ix = find(isfinite(xi_frf));

datadir = '/Volumes/ThunderBay/fC_RODSEX_0928_D3_dx1p3/';
% datadir = '~/Dropbox/RODSEX/funwaveC/';

sd = length(ix);
if mod(sd,2)==0               % only half the data is good
    stop=sd/2+1;                % we are making a onesided spectrum
else						  % modified on 11/11/03
    stop=(sd+1)/2;
end
numfiles = 1900;
wavenum_spec_vort_xshore = zeros(numy,stop);
n=0;
for ii = 3599:(3599+numfiles)
    n=n+1;
    fprintf('%g\n',ii)
    load(sprintf('%ssnap_vort_l2_%4.0f.mat',datadir,ii))
    [mpsd,wavenums_x]=mypsd(vort(ix,:),1/dx);
    wavenum_spec_vort_xshore = wavenum_spec_vort_xshore + mpsd.'/numfiles;

    
end

%%
load vorticity_wavenumbers_out.mat

save vorticity_wavenumbers_out.mat wavenums wavenum_spec_vort vort_ring avg_i i_index xspec RING stdvort rsize wavenum_spec_vort_xshore wavenums_x ix
%%

figure(4); clf
pcolor(wavenums_x,(1:numy)*dy,wavenum_spec_vort_xshore)
shading flat




%% combine cross- and alognshore wavenumber spectra

wavenum_spec_vort_xshore_interp = interp1(wavenums_x,mean(wavenum_spec_vort_xshore),wavenums,'linear','extrap');

wavenum_spec_vort_combined = wavenum_spec_vort_xshore_interp + mean(wavenum_spec_vort);

figure(5);clf
plot(wavenums_x(2:end),(mean(wavenum_spec_vort_xshore(:,2:end))))
hold on
plot(wavenums(2:end),(mean(wavenum_spec_vort(:,2:end))))
plot(wavenums(2:end),wavenum_spec_vort_combined(2:end))
xlabel('wavenumber (1/m)')
ylabel('vorticity spectral density')
legend('x-shore','alongshore','combined')

dk = wavenums(2)-wavenums(1);
[dump,i5m] = min(abs(wavenums-(1/(5*dx))));
ring_vort_variance_from_wavenumber_spectra = sum(wavenum_spec_vort_combined(2:i5m))*dk
total_vort_variance_from_wavenumber_spectra = sum(wavenum_spec_vort_combined(2:end))*dk

%% "surfzone" cross-shore wavenumber spectra

ixs = 510-20:510+20;

datadir = '/Volumes/ThunderBay/fC_RODSEX_0928_D3_dx1p3/';
% datadir = '~/Dropbox/RODSEX/funwaveC/';

sd = length(ixs);
if mod(sd,2)==0               % only half the data is good
    stop=sd/2+1;                % we are making a onesided spectrum
else						  % modified on 11/11/03
    stop=(sd+1)/2;
end
numfiles = 1900;
wavenum_spec_vort_xshore_short = zeros(numy,stop);
n=0;
for ii = 3599:(3599+numfiles)
    fprintf('%g\n',ii)
    load(sprintf('%ssnap_vort_l2_%4.0f.mat',datadir,ii))
    [mpsd,wavenums_x_short]=mypsd(vort(ixs,:),1/dx);
    wavenum_spec_vort_xshore_short = wavenum_spec_vort_xshore_short + mpsd.'/numfiles;
end


%%

figure(6); clf
pcolor(wavenums_x_short,(1:numy)*dy,wavenum_spec_vort_xshore_short)
shading flat

sum(mean(wavenum_spec_vort_xshore_short))*(wavenums_x_short(2)-wavenums_x_short(1))



figure(7);clf
plot(wavenums_x(2:end),(mean(wavenum_spec_vort_xshore(:,2:end))))
hold on
plot(wavenums(2:end),(mean(wavenum_spec_vort(:,2:end))))
plot(wavenums(2:end),wavenum_spec_vort_combined(2:end))
plot(wavenums_x_short(2:end),mean(wavenum_spec_vort_xshore_short(:,2:end)),'mo')
xlabel('wavenumber (1/m)')
ylabel('vorticity spectral density')
legend('x-shore','alongshore','combined','xshort')

%%

short_for_interp = mean(wavenum_spec_vort_xshore_short);
short_for_interp(1) = short_for_interp(2); % remove mean


wavenum_spec_vort_xshore_short_interp = interp1(wavenums_x_short,short_for_interp,wavenums,'linear','extrap');
wavenum_spec_vort_xshore_short_interp(1) = mean(wavenum_spec_vort_xshore_short(:,1)); % put the mean back in


wavenum_spec_vort_short_combined = wavenum_spec_vort_xshore_short_interp + mean(wavenum_spec_vort(ixs,:));

%%

%%
save vorticity_wavenumbers_out.mat wavenums wavenum_spec_vort vort_ring avg_i i_index xspec RING stdvort rsize wavenum_spec_vort_xshore wavenums_x ix wavenum_spec_vort_xshore_short wavenums_x_short wavenum_spec_vort_short_combined


%%
load vorticity_wavenumbers_out.mat

ring_vort_variance_from_wavenumber_spectra_short = sum(wavenum_spec_vort_short_combined(2:i5m))*dk
total_vort_variance_from_wavenumber_spectra_short = sum(wavenum_spec_vort_short_combined(2:end))*dk


% STILL NEED TO CHECK THAT VORTICITY FROM A SINGLE "RING" LOCATION IS
% SIMILAR TO USING MANY ALONGSHORE LOCATIONS


%% animation
figure(8); clf
load(sprintf('%ssnap_vort_l2_%4.0f.mat',datadir,3599))
hp = pcolor((1:1200)*dy,xi_frf,vort);
shading flat
axis equal
caxis([-.2 .2])
colormap jet
% 
% for ii = 3599:(3599+numfiles)
%     fprintf('%g\n',ii)
%     load(sprintf('%ssnap_vort_l2_%4.0f.mat',datadir,ii))
%     set(hp,'cdata',vort)
%     drawnow
% end


for ii = 3599:10:(3599+numfiles)
    fprintf('%g\n',ii)
    vort_tmp = zeros(size(vort));
    for jj = ii:(ii+10);
    load(sprintf('%ssnap_vort_l2_%4.0f.mat',datadir,jj))
    vort_tmp = vort_tmp + vort/10;
    end
    set(hp,'cdata',vort_tmp)
    drawnow
end

%% compare ring size vs vorticity variance with wavenumber range vs the
% spectrum integrated over that range (variance)

spec_partial_std = zeros(size(rsize));
spec_partial_var_wavenum_limit = 1./(rsize);
wavenum_limit =  zeros(size(rsize));
dk = wavenums(2) - wavenums(1);
for ii = 1:length(rsize)
    [dump,imax] = min(abs(wavenums-spec_partial_var_wavenum_limit(ii)));
    wavenum_limit(ii) = wavenums(imax);
    spec_partial_std(ii) = sqrt(sum(wavenum_spec_vort_short_combined(2:imax)*dk));
end

figure(9); clf
plot(spec_partial_var_wavenum_limit, stdvort, wavenum_limit, spec_partial_std)
xlabel('1/ringSize, wavenum\_limit (1/m)')
ylabel('vort std')