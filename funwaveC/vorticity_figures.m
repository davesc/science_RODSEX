%% vorticity wavenumber plots

load ~/Dropbox/RODSEX/funwaveC/vorticity_wavenumbers_spectral_sum_comparison_data.mat

%% alongshore wavenumber spectrum

figure(1); clf
% note when spectrum is processed with mywelch instead of mypsd need
% wavenum_spec_vort(2:end-1,) to avoid nans in pcolor
pcolor(wavenums(2:end),xi_frf(2:end-1),log10(wavenum_spec_vort(2:end-1,2:end))); 
% pcolor(wavenums(2:end),xi_frf,log10(wavenum_spec_vort(:,2:end))); 
colormap jet; 
colorbar;
shading flat;
ylabel('FRF xshore (m)')
xlabel('vort wavenumber spectrum (m/s^2)')
caxis([-7 -.5])
hold on
plot([0 0.7],[xi_frf(510), xi_frf(510)],'--k')
text(0.55,xi_frf(510)+10,'ring location')

%% wavenumber spectra

figure(2); clf
plot(wavenums,mean(wavenum_spec_vort(ixs,:)), ...
    wavenums_x_short,mean(wavenum_spec_vort_xshore_short), ...
    wavenums,wavenum_spec_vort_xshore_short_interp,'--')
legend('along','cross','cross-interp')
title('xshore vort wavenum spectra, centered on ring +-30m')
ylabel('vorticity spectra (m/s^2)')
xlabel('wavenumbers (1/m)')


%% ring spatial-averaged vorticity 

figure(3); clf

subplot(121)
plot(ring_wavenumber, varvort, ...
     wavenum_limit, spec_partial_var_xshore_only ,wavenum_limit, spec_partial_var_alongshore_only)
xlabel('1/ringSize, wavenum\_limit (1/m)')
ylabel('vort variance (1/s^2)')
legend('ring average','spectrum cross-shore','spectrum alognshore','location','southeast')

subplot(122)
interp_varvort = interp1(ring_wavenumber, varvort, wavenum_limit);
plot(wavenum_limit, spec_partial_var_alongshore_only./interp_varvort, ...
     wavenum_limit, spec_partial_var_xshore_only./interp_varvort)
xlabel('wavenumber(1/m)')
ylabel('spectralSum / ringAverage ')
legend('alongshore','cross-shore')


%% snap vorticity plot, quick n dirty animation from old files on Dropbox

datadir = '~/Dropbox/RODSEX/funwaveC/';

figure(4); clf

for ii = 1:100
load(sprintf('%ssnap_vort_l2_%4.0f.mat',datadir,3599+ii))


pcolor((1:1200)*dy,xi_frf,vort);
axis([1 400 100 500])
shading flat
colormap jet
% axis equal 
colorbar
caxis([-.2, .2])
pause(.1)

end

%% quick-n-dirty sea-surface elevation animation

eta = load_fCbinary_file('/Volumes/ThunderBay/fC_RODSEX_0928_D3_dx1p3/snap_eta_snap_l2_3600.fCdat');
etasmall = eta(:,1:5:21);
eddy = load_fCbinary_file('/Volumes/ThunderBay/fC_RODSEX_0928_D3_dx1p3/snap_eddy_snap_l2_3600.fCdat');
ind_eddy = find(eddy(:,1:5:21));
xind = repmat(1:size(eddy,1),1,size(eddy,2));
h = load('/Volumes/ThunderBay/fC_RODSEX_0928_D3_dx1p3/bathy_RODSEX_0925_09281600_1D_D_dx1p3.depth');

figure(1);
clf
plot(-h); 
hold on
hp1 = plot(etasmall);
hp2 = plot(xind(ind_eddy),etasmall(ind_eddy),'ko');
title('eta, 5 cross-shore lines, 6.67 m spacing, circles indicate breaking eddy viscocity')
xlabel('xshore bins')
ylabel('eta, surface elevation (m)')

for ii = 3600:4000;
eta = load_fCbinary_file(sprintf('/Volumes/ThunderBay/fC_RODSEX_0928_D3_dx1p3/snap_eta_snap_l2_%d.fCdat',ii));
etasmall = eta(:,1:5:21);
eddy = load_fCbinary_file(sprintf('/Volumes/ThunderBay/fC_RODSEX_0928_D3_dx1p3/snap_eddy_snap_l2_%d.fCdat',ii));
ind_eddy = find(eddy(:,1:5:21));
hold off
plot(-h); 
hold on
plot(etasmall);
plot(xind(ind_eddy),etasmall(ind_eddy),'ko');
title('eta, 5 cross-shore lines, 6.67 m spacing, circles indicate breaking eddy viscosity')
xlabel('xshore bins')
ylabel('eta, surface elevation (m)')
% set(hp1,'ydata',etasmall)
% set(hp2,'xdata',xind(ind_eddy),'ydata',etasmall(ind_eddy))
pause(.15)

end





