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

%%

