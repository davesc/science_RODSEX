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
            RING(jj).vort(n,kk+1) = sum(sum(vort(RING(jj).i,RING(jj).j + (RING(jj).bin_length*kk))))/(RING(jj).bin_length.^2);
        end
    end
 
end


