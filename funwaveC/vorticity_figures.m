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
caxis([-7 0])
hold on
plot([0 0.7],[xi_frf(510), xi_frf(510)],'--k')
text(0.55,xi_frf(510)+10,'ring location')

%% xshore wavenumber spectrum

figure(2); clf
plot(wavenums_x_short,mean(wavenum_spec_vort_xshore_short))
title('xshore vort wavenum spectra, centered on ring +-30m')
ylabel('vorticity spectra (m/s^2)')
xlabel('wavenumbers (1/m)')


%% ring spatial-averaged vorticity 

figure(3); clf

subplot(121)
plot(spec_partial_var_wavenum_limit, stdvort, ...
     wavenum_limit, spec_partial_std)
xlabel('1/ringSize, wavenum\_limit (1/m)')
ylabel('vort std (1/s)')
legend('ring average','spectrum','location','southeast')

subplot(122)
interp_stdvort = interp1(spec_partial_var_wavenum_limit, stdvort, wavenum_limit);
plot(wavenum_limit, spec_partial_std./interp_stdvort)
xlabel('wavenumber(1/m)')
ylabel('spectralSum / ringAverage ')







