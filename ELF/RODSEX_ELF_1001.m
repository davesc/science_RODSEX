%% RODSEX ELF spectra using fft on 10/01, 2013


PSD = struct;

% data_dir = '/Volumes/OWC HD/RODSEX';
data_dir = '/Volumes/DAVID_CLARK/RODSEX';

days = ['1001'];
    start_hour = 0;
    end_hour = 23;

    inst = [23
            33
            43
            63
            93];
    
    Lrec = 6120; %(51 minutes of 2 Hz data), the actual record is 51.2 minutes
    Lrec2  = 60; % 60 minutes of 60 averaged data, to be filled / interpolated
    
    p = zeros(Lrec2*24*size(days,1),1);     % presusre vector
    u = zeros(Lrec2*24*size(days,1),1);     % cross-shore velocity vector
    v = zeros(Lrec2*24*size(days,1),1);     % alongshore velocity vector
    pbar = zeros(24*size(days,1),1);
    ubar = zeros(24*size(days,1),1);
    vbar = zeros(24*size(days,1),1);
    
for kk = 1:length(inst);
    N = 0; % counter for hours since I'm crossing between days
    for ii=1:size(days,1);
        if ii==1;
            sth = start_hour;
        else
            sth = 1;
        end
        if ii==size(days,1);
            enh = end_hour;
        else
            enh=24;
        end
        for jj = sth:enh; % for each hour of the day
            N=N+1;
            % grab 51 minutes of data, make 60 second averages
            % fill the 9 minutes at the end of the hour with NaNs
            
            % load 51.2 minutes of data
            try
            ptmp = load_elgar_puv_binary(sprintf(...
                '%s/%04.0f/%04.0f%02.0f00.p%02.0f',...
                data_dir,str2num(days(ii,:)),str2num(days(ii,:)),jj-1,inst(kk)));
            catch
                ptmp = NaN(Lrec,1);
            end
            try
            utmp = load_elgar_puv_binary(sprintf(...
                '%s/%04.0f/%04.0f%02.0f00.u%02.0f',...
                data_dir,str2num(days(ii,:)),str2num(days(ii,:)),jj-1,inst(kk)));
            catch
                utmp = NaN(Lrec,1);
            end
            try
            vtmp = load_elgar_puv_binary(sprintf(...
                '%s/%04.0f/%04.0f%02.0f00.v%02.0f',...
                data_dir,str2num(days(ii,:)),str2num(days(ii,:)),jj-1,inst(kk)));
            catch
                vtmp = NaN(Lrec,1);
            end
            
            % take one minute means, and fill the last 9 minutes of the hour with NaNs
            p((1:Lrec2)+Lrec2*(N-1)) = [mean(reshape(ptmp(1:Lrec),[120 51])).'; NaN(9,1)];
            u((1:Lrec2)+Lrec2*(N-1)) = [mean(reshape(utmp(1:Lrec),[120 51])).'; NaN(9,1)];
            v((1:Lrec2)+Lrec2*(N-1)) = [mean(reshape(vtmp(1:Lrec),[120 51])).'; NaN(9,1)];
            
            pbar(N) = nanmean(ptmp);
            ubar(N) = nanmean(utmp);
            vbar(N) = nanmean(vtmp);

        end
    end
        
        
%         % Lomb-Scargle spectrum, with 127 minute chunks, 75% overlap, linear
%         % detrending
%         
%         if sum(isfinite(p)) > sum(isnan(p))*5; % only process spectra if there is enough good data
%             [psd_p,fr]=mywelch_plomb(p,60,42,0.75,1/120);
%         else
%             psd_p = NaN;
%         end
%         
%         if sum(isfinite(u)) > sum(isnan(u))*5; % only process spectra if there is enough good data
%             [psd_u,fr]=mywelch_plomb(u,60,42,0.75,1/120);
%         else
%             psd_u = NaN;
%         end
%         
%         if sum(isfinite(v)) > sum(isnan(v))*5; % only process spectra if there is enough good data
%             [psd_v,fr]=mywelch_plomb(v,60,42,0.75,1/120);
%         else
%             psd_v = NaN;
%         end
%         
%         % Lomb-Scargle spectrum, with 127 minute chunks, 75% overlap, QUADRATIC
%         % detrending        
%         if sum(isfinite(p)) > sum(isnan(p))*5; % only process spectra if there is enough good data
%             [psd_p_quad,fr]=mywelch_plomb_quad(p,60,42,0.75,1/120);
%         else
%             psd_p_quad = NaN;
%         end
%         
%         if sum(isfinite(u)) > sum(isnan(u))*5; % only process spectra if there is enough good data
%             [psd_u_quad,fr]=mywelch_plomb_quad(u,60,42,0.75,1/120);
%         else
%             psd_u_quad = NaN;
%         end
%         
%         if sum(isfinite(v)) > sum(isnan(v))*5; % only process spectra if there is enough good data
%             [psd_v_quad,fr]=mywelch_plomb_quad(v,60,42,0.75,1/120);
%         else
%             psd_v_quad = NaN;
%         end


% fft spectrum, with 127 minute chunks, 75% overlap, linear
        % detrending
        
        if sum(isfinite(p)) > sum(isnan(p)); % only process spectra if there is enough good data
            Inan = find(isnan(p));
            Irl = find(isfinite(p));
            p(Inan) = interp1(Irl,p(Irl),Inan,'linear',0);
            [psd_p,fr_p]=mywelch(p,60,42*size(days,1),0.75);
            [psd_p_quad,fr_pq]=mywelch_quad(p,60,42*size(days,1),0.75);
        else
            psd_p = NaN;
            psd_p_quad = NaN;
            disp(sprintf('bad p, %2.0f',inst(kk)))
        end
        
        if sum(isfinite(u)) > sum(isnan(u)); % only process spectra if there is enough good data
            Inan = find(isnan(u));
            Irl = find(isfinite(u));
            u(Inan) = interp1(Irl,u(Irl),Inan,'linear',0);
            [psd_u,fr_u]=mywelch(u,60,42*size(days,1),0.75);
            [psd_u_quad,fr_uq]=mywelch_quad(u,60,42*size(days,1),0.75);
        else
            psd_u = NaN;
            psd_u_quad = NaN;
            disp(sprintf('bad u, %2.0f',inst(kk)))
        end
        
        if sum(isfinite(v)) > sum(isnan(v)); % only process spectra if there is enough good data
            Inan = find(isnan(v));
            Irl = find(isfinite(v));
            v(Inan) = interp1(Irl,v(Irl),Inan,'linear',0);
            [psd_v,fr_v]=mywelch(v,60,42*size(days,1),0.75);
            [psd_v_quad,fr_vq]=mywelch_quad(v,60,42*size(days,1),0.75);
        else
            psd_v = NaN;
            psd_v_quad = NaN;
            disp(sprintf('bad v, %2.0f',inst(kk)))
        end
        


    
    PSD(kk).date = days(ii,:);
    PSD(kk).freq_p = fr_p;
    PSD(kk).freq_pq = fr_pq;
    PSD(kk).freq_u = fr_u;
    PSD(kk).freq_uq = fr_uq;
    PSD(kk).freq_v = fr_v;
    PSD(kk).freq_vq = fr_vq;
    PSD(kk).v = psd_v;
    PSD(kk).u = psd_u;
    PSD(kk).p = psd_p;
    PSD(kk).v_quadratic_detrend = psd_v_quad;
    PSD(kk).u_quadratic_detrend = psd_u_quad;
    PSD(kk).p_quadratic_detrend = psd_p_quad;
    PSD(kk).inst = inst(kk);
    PSD(kk).pbar = pbar;
    PSD(kk).ubar = ubar;
    PSD(kk).vbar = vbar;
%     
%               if kk==2;
%                 return
%             end  
 end
    
    
    
    
    
%%

figure(4); clf
clrs = hsv(length(inst));
leg = {};
vbar = 0;
hp = [];

subplot 311
for kk = 1:length(inst);
hptmp = plot(PSD(kk).freq_u,PSD(kk).u,'color',clrs(kk,:),'linewidth',1.5)
hp(kk) = hptmp(1);
hold on
leg{kk} = sprintf('q#%2.0f, V=%1.02f, U=%1.02fm/s',inst(kk),nanmean(PSD(kk).vbar), nanmean(PSD(kk).ubar));
vbar = vbar + nanmean(PSD(kk).vbar);
end
vbar = vbar/kk;
ylabel('U psd (m^2/s)')
xlabel('frequency (Hz)')



legend(hp,leg)
xlim([0 .005])
title(sprintf('RODSEX 10/01')); %, Vbar = %01.02f m/s',vbar))

subplot 312
for kk = 1:length(inst);
plot(PSD(kk).freq_v,PSD(kk).v,'color',clrs(kk,:),'linewidth',1.5)
hold on
end
ylabel('V psd (m^2/s)')
xlabel('frequency (Hz)')
    xlim([0 .005])
    
subplot 313
for kk = 1:length(inst);
plot(PSD(kk).freq_v,PSD(kk).v+PSD(kk).u,'color',clrs(kk,:),'linewidth',1.5)
hold on
end
ylabel('V+U psd (m^2/s)')
xlabel('frequency (Hz)')
    xlim([0 .005])
    
    
orient tall
print -dpdf RODSEX_ELF_1001_3sLine.pdf
