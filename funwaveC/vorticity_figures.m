%% vorticity wavenumber plots

load ~/Dropbox/RODSEX/funwaveC/vorticity_wavenumbers_spectral_sum_comparison_data.mat

%% alongshore wavenumber spectrum

figure(1); clf
pcolor(wavenums(2:end),xi_frf,log10(wavenum_spec_vort(:,2:end))); 
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
plot(wavenum_limit, spec_partial_var./interp_varvort)
xlabel('wavenumber(1/m)')
ylabel('spectralSum / ringAverage ')



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

%%





