%% RODSEX ELF spectra using lomb scargle on 0927, 2013


PSD = struct;

data_dir = '/Volumes/OWC HD/RODSEX/ring_and_standalone_data';
% data_dir = '/Volumes/DAVID_CLARK/RODSEX';

days = ['0927'];
start_hour = 0;
end_hour = 23;

% length of record in hours
if size(days,1) == 1;
    rec_hours = end_hour - start_hour + 1;
elseif size(days,1) > 1;
    rec_hours = (size(days,1)-1)*24 - start_hour + 1 + end_hour
else
    rec_hours = NaN;
end


inst = [21
    22
    23
    24];


Lrec = 6120; %(51 minutes of 2 Hz data), the actual record is 51.2 minutes
Lrec2  = 60; % 60 minutes of 60 averaged data, to be filled / interpolated
Lrec8 = 24480; % 51 min of 8 Hz data

%     p = zeros(Lrec2*24*size(days,1),1);     % presusre vector
%     u = zeros(Lrec2*24*size(days,1),1);     % cross-shore velocity vector
%     v = zeros(Lrec2*24*size(days,1),1);     % alongshore velocity vector
%     pbar = zeros(24*size(days,1),1);
%     ubar = zeros(24*size(days,1),1);
%     vbar = zeros(24*size(days,1),1);



p = NaN(Lrec2*rec_hours,1);     % presusre vector
u = NaN(Lrec2*rec_hours,1);     % cross-shore velocity vector
v = NaN(Lrec2*rec_hours,1);     % alongshore velocity vector
pbar = NaN(24*rec_hours,1);
ubar = NaN(24*rec_hours,1);
vbar = NaN(24*rec_hours,1);



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
                if inst(kk)<15 % if data is from a ring adv
                    utmp = mean(reshape(utmp(1:Lrec8),[4 (Lrec8/4)])).'; % average 8Hz down to 2Hz
                end
                    
            catch
                utmp = NaN(Lrec,1);
            end
            try
                vtmp = load_elgar_puv_binary(sprintf(...
                    '%s/%04.0f/%04.0f%02.0f00.v%02.0f',...
                    data_dir,str2num(days(ii,:)),str2num(days(ii,:)),jj-1,inst(kk)));
                if inst(kk)<15 % if data is from a ring adv
                    vtmp = mean(reshape(vtmp(1:Lrec8),[4 (Lrec8/4)])).'; % average 8Hz down to 2Hz
                end
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
    
    % number of data chunks to use for spectrum estimate
    chnks = 1.75 * rec_hours;
    
    % Lomb-Scargle spectrum, with 127 minute chunks, 75% overlap, linear
    % detrending
    
    if sum(isfinite(p)) > sum(isnan(p))*2; % only process spectra if there is enough good data
        [psd_p,fr_p]=mywelch_plomb(p,60,chnks,0.75,1/120);
        [psd_p_quad,fr_pq]=mywelch_plomb_quad(p,60,chnks,0.75,1/120);
    else
        psd_p = NaN;
        psd_p_quad = NaN;
        fr_p = NaN;
        fr_pq = NaN;
        disp(sprintf('bad p, %2.0f',inst(kk)))
    end
    
    if sum(isfinite(u)) > sum(isnan(u))*2; % only process spectra if there is enough good data
        [psd_u,fr_u]=mywelch_plomb(u,60,chnks,0.75,1/120);
        [psd_u_quad,fr_uq]=mywelch_plomb_quad(u,60,chnks,0.75,1/120);
    else
        psd_u = NaN;
        psd_u_quad = NaN;
        disp(sprintf('bad u, %2.0f',inst(kk)))
        fr_u = NaN;
        fr_uq = NaN;
    end
    
    if sum(isfinite(v)) > sum(isnan(v))*2; % only process spectra if there is enough good data
        [psd_v,fr_v]=mywelch_plomb(v,60,chnks,0.75,1/120);
        [psd_v_quad,fr_vq]=mywelch_plomb_quad(v,60,chnks,0.75,1/120);
    else
        psd_v = NaN;
        psd_v_quad = NaN;
        fr_v = NaN;
        fr_vq = NaN;
        disp(sprintf('bad v, %2.0f',inst(kk)))
    end
    
    %
    % % fft spectrum, with 127 minute chunks, 75% overlap, linear
    %         % detrending
    %
    %         if sum(isfinite(p)) > sum(isnan(p)); % only process spectra if there is enough good data
    %             Inan = find(isnan(p));
    %             Irl = find(isfinite(p));
    %             p(Inan) = interp1(Irl,p(Irl),Inan,'linear',0);
    %             [psd_p,fr_p]=mywelch(p,60,chnks*size(days,1),0.75);
    %             [psd_p_quad,fr_pq]=mywelch_quad(p,60,chnks*size(days,1),0.75);
    %         else
    %             psd_p = NaN;
    %             psd_p_quad = NaN;
    %             fr_p = NaN;
    %             fr_pq = NaN
    %             disp(sprintf('bad p, %2.0f',inst(kk)))
    %         end
    %
    %         if sum(isfinite(u)) > sum(isnan(u)); % only process spectra if there is enough good data
    %             Inan = find(isnan(u));
    %             Irl = find(isfinite(u));
    %             u(Inan) = interp1(Irl,u(Irl),Inan,'linear',0);
    %             [psd_u,fr_u]=mywelch(u,60,chnks*size(days,1),0.75);
    %             [psd_u_quad,fr_uq]=mywelch_quad(u,60,chnks*size(days,1),0.75);
    %         else
    %             psd_u = NaN;
    %             psd_u_quad = NaN;
    %             disp(sprintf('bad u, %2.0f',inst(kk)))
    %             fr_u = NaN;
    %             fr_uq = NaN
    %         end
    %
    %         if sum(isfinite(v)) > sum(isnan(v)); % only process spectra if there is enough good data
    %             Inan = find(isnan(v));
    %             Irl = find(isfinite(v));
    %             v(Inan) = interp1(Irl,v(Irl),Inan,'linear',0);
    %             [psd_v,fr_v]=mywelch(v,60,chnks*size(days,1),0.75);
    %             [psd_v_quad,fr_vq]=mywelch_quad(v,60,chnks*size(days,1),0.75);
    %         else
    %             psd_v = NaN;
    %             psd_v_quad = NaN;
    %             fr_v = NaN;
    %             fr_vq = NaN;
    %             disp(sprintf('bad v, %2.0f',inst(kk)))
    %         end
    %
    
    
    
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

subplot 211
N = 0;
for kk = 1:length(inst);
    hptmp = plot(PSD(kk).freq_u,PSD(kk).u,'color',clrs(kk,:),'linewidth',1.5);
    hp(kk) = hptmp(1);
    hold on
    leg{kk} = sprintf('#%02.0f, V=%01.2f U=%01.2f',inst(kk),nanmean(PSD(kk).vbar),nanmean(PSD(kk).ubar));
    vbar_tmp = nanmean(PSD(kk).vbar);
    if isfinite(vbar_tmp);
        vbar = vbar + vbar_tmp;
        N = N+1;
    end
end
vbar = vbar/N;
ylabel('U psd (m^2/s)')
xlabel('frequency (Hz)')



legend(hp,leg)
xlim([0 .005])
title(sprintf('RODSEX %s/%s, Vbar = %01.02f m/s',days(1,1:2),days(1,3:4),vbar))

subplot 212
for kk = 1:length(inst);
    plot(PSD(kk).freq_v,PSD(kk).v,'color',clrs(kk,:),'linewidth',1.5)
    hold on
end
ylabel('V psd (m^2/s)')
xlabel('frequency (Hz)')
xlim([0 .005])


orient tall
eval(sprintf('print -dpdf RODSEX_ELF_%s_xshore_20sline.pdf',days(1,:)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% RODSEX ELF spectra using fft on 1009, 2013


PSD = struct;



inst = [51
    52
    07
    55
    54
    56];


Lrec = 6120; %(51 minutes of 2 Hz data), the actual record is 51.2 minutes
Lrec2  = 60; % 60 minutes of 60 averaged data, to be filled / interpolated

%     p = zeros(Lrec2*24*size(days,1),1);     % presusre vector
%     u = zeros(Lrec2*24*size(days,1),1);     % cross-shore velocity vector
%     v = zeros(Lrec2*24*size(days,1),1);     % alongshore velocity vector
%     pbar = zeros(24*size(days,1),1);
%     ubar = zeros(24*size(days,1),1);
%     vbar = zeros(24*size(days,1),1);



p = NaN(Lrec2*rec_hours,1);     % presusre vector
u = NaN(Lrec2*rec_hours,1);     % cross-shore velocity vector
v = NaN(Lrec2*rec_hours,1);     % alongshore velocity vector
pbar = NaN(24*rec_hours,1);
ubar = NaN(24*rec_hours,1);
vbar = NaN(24*rec_hours,1);



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
                if inst(kk)<15 % if data is from a ring adv
                    utmp = mean(reshape(utmp(1:Lrec8),[4 (Lrec8/4)])).'; % average 8Hz down to 2Hz
                end
                    
            catch
                utmp = NaN(Lrec,1);
            end
            try
                vtmp = load_elgar_puv_binary(sprintf(...
                    '%s/%04.0f/%04.0f%02.0f00.v%02.0f',...
                    data_dir,str2num(days(ii,:)),str2num(days(ii,:)),jj-1,inst(kk)));
                if inst(kk)<15 % if data is from a ring adv
                    vtmp = mean(reshape(vtmp(1:Lrec8),[4 (Lrec8/4)])).'; % average 8Hz down to 2Hz
                end
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
    
    % number of data chunks to use for spectrum estimate
    chnks = 1.75 * rec_hours;
    
    % Lomb-Scargle spectrum, with 127 minute chunks, 75% overlap, linear
    % detrending
    
    if sum(isfinite(p)) > sum(isnan(p))*2; % only process spectra if there is enough good data
        [psd_p,fr_p]=mywelch_plomb(p,60,chnks,0.75,1/120);
        [psd_p_quad,fr_pq]=mywelch_plomb_quad(p,60,chnks,0.75,1/120);
    else
        psd_p = NaN;
        psd_p_quad = NaN;
        fr_p = NaN;
        fr_pq = NaN;
        disp(sprintf('bad p, %2.0f',inst(kk)))
    end
    
    if sum(isfinite(u)) > sum(isnan(u))*2; % only process spectra if there is enough good data
        [psd_u,fr_u]=mywelch_plomb(u,60,chnks,0.75,1/120);
        [psd_u_quad,fr_uq]=mywelch_plomb_quad(u,60,chnks,0.75,1/120);
    else
        psd_u = NaN;
        psd_u_quad = NaN;
        disp(sprintf('bad u, %2.0f',inst(kk)))
        fr_u = NaN;
        fr_uq = NaN;
    end
    
    if sum(isfinite(v)) > sum(isnan(v))*2; % only process spectra if there is enough good data
        [psd_v,fr_v]=mywelch_plomb(v,60,chnks,0.75,1/120);
        [psd_v_quad,fr_vq]=mywelch_plomb_quad(v,60,chnks,0.75,1/120);
    else
        psd_v = NaN;
        psd_v_quad = NaN;
        fr_v = NaN;
        fr_vq = NaN;
        disp(sprintf('bad v, %2.0f',inst(kk)))
    end
    
    %
    % % fft spectrum, with 127 minute chunks, 75% overlap, linear
    %         % detrending
    %
    %         if sum(isfinite(p)) > sum(isnan(p)); % only process spectra if there is enough good data
    %             Inan = find(isnan(p));
    %             Irl = find(isfinite(p));
    %             p(Inan) = interp1(Irl,p(Irl),Inan,'linear',0);
    %             [psd_p,fr_p]=mywelch(p,60,chnks*size(days,1),0.75);
    %             [psd_p_quad,fr_pq]=mywelch_quad(p,60,chnks*size(days,1),0.75);
    %         else
    %             psd_p = NaN;
    %             psd_p_quad = NaN;
    %             fr_p = NaN;
    %             fr_pq = NaN
    %             disp(sprintf('bad p, %2.0f',inst(kk)))
    %         end
    %
    %         if sum(isfinite(u)) > sum(isnan(u)); % only process spectra if there is enough good data
    %             Inan = find(isnan(u));
    %             Irl = find(isfinite(u));
    %             u(Inan) = interp1(Irl,u(Irl),Inan,'linear',0);
    %             [psd_u,fr_u]=mywelch(u,60,chnks*size(days,1),0.75);
    %             [psd_u_quad,fr_uq]=mywelch_quad(u,60,chnks*size(days,1),0.75);
    %         else
    %             psd_u = NaN;
    %             psd_u_quad = NaN;
    %             disp(sprintf('bad u, %2.0f',inst(kk)))
    %             fr_u = NaN;
    %             fr_uq = NaN
    %         end
    %
    %         if sum(isfinite(v)) > sum(isnan(v)); % only process spectra if there is enough good data
    %             Inan = find(isnan(v));
    %             Irl = find(isfinite(v));
    %             v(Inan) = interp1(Irl,v(Irl),Inan,'linear',0);
    %             [psd_v,fr_v]=mywelch(v,60,chnks*size(days,1),0.75);
    %             [psd_v_quad,fr_vq]=mywelch_quad(v,60,chnks*size(days,1),0.75);
    %         else
    %             psd_v = NaN;
    %             psd_v_quad = NaN;
    %             fr_v = NaN;
    %             fr_vq = NaN;
    %             disp(sprintf('bad v, %2.0f',inst(kk)))
    %         end
    %
    
    
    
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

figure(5); clf
clrs = hsv(length(inst));
leg = {};
vbar = 0;
hp = [];

subplot 211
N = 0;
for kk = 1:length(inst);
    hptmp = plot(PSD(kk).freq_u,PSD(kk).u,'color',clrs(kk,:),'linewidth',1.5);
    hp(kk) = hptmp(1);
    hold on
    leg{kk} = sprintf('#%02.0f, V=%01.2f U=%01.2f',inst(kk),nanmean(PSD(kk).vbar),nanmean(PSD(kk).ubar));
    vbar_tmp = nanmean(PSD(kk).vbar);
    if isfinite(vbar_tmp);
        vbar = vbar + vbar_tmp;
        N = N+1;
    end
end
vbar = vbar/N;
ylabel('U psd (m^2/s)')
xlabel('frequency (Hz)')



legend(hp,leg)
xlim([0 .005])
title(sprintf('RODSEX %s/%s, Vbar = %01.02f m/s',days(1,1:2),days(1,3:4),vbar))

subplot 212
for kk = 1:length(inst);
    plot(PSD(kk).freq_v,PSD(kk).v,'color',clrs(kk,:),'linewidth',1.5)
    hold on
end
ylabel('V psd (m^2/s)')
xlabel('frequency (Hz)')
xlim([0 .005])


orient tall
eval(sprintf('print -dpdf RODSEX_ELF_%s_xshore_50sline.pdf',days(1,:)))




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%% RODSEX ELF spectra using fft on 1009, 2013


PSD = struct;


inst = [71
    72
    73
    74];


Lrec = 6120; %(51 minutes of 2 Hz data), the actual record is 51.2 minutes
Lrec2  = 60; % 60 minutes of 60 averaged data, to be filled / interpolated

%     p = zeros(Lrec2*24*size(days,1),1);     % presusre vector
%     u = zeros(Lrec2*24*size(days,1),1);     % cross-shore velocity vector
%     v = zeros(Lrec2*24*size(days,1),1);     % alongshore velocity vector
%     pbar = zeros(24*size(days,1),1);
%     ubar = zeros(24*size(days,1),1);
%     vbar = zeros(24*size(days,1),1);



p = NaN(Lrec2*rec_hours,1);     % presusre vector
u = NaN(Lrec2*rec_hours,1);     % cross-shore velocity vector
v = NaN(Lrec2*rec_hours,1);     % alongshore velocity vector
pbar = NaN(24*rec_hours,1);
ubar = NaN(24*rec_hours,1);
vbar = NaN(24*rec_hours,1);



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
                if inst(kk)<15 % if data is from a ring adv
                    utmp = mean(reshape(utmp(1:Lrec8),[4 (Lrec8/4)])).'; % average 8Hz down to 2Hz
                end
                    
            catch
                utmp = NaN(Lrec,1);
            end
            try
                vtmp = load_elgar_puv_binary(sprintf(...
                    '%s/%04.0f/%04.0f%02.0f00.v%02.0f',...
                    data_dir,str2num(days(ii,:)),str2num(days(ii,:)),jj-1,inst(kk)));
                if inst(kk)<15 % if data is from a ring adv
                    vtmp = mean(reshape(vtmp(1:Lrec8),[4 (Lrec8/4)])).'; % average 8Hz down to 2Hz
                end
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
    
    % number of data chunks to use for spectrum estimate
    chnks = 1.75 * rec_hours;
    
    % Lomb-Scargle spectrum, with 127 minute chunks, 75% overlap, linear
    % detrending
    
    if sum(isfinite(p)) > sum(isnan(p))*2; % only process spectra if there is enough good data
        [psd_p,fr_p]=mywelch_plomb(p,60,chnks,0.75,1/120);
        [psd_p_quad,fr_pq]=mywelch_plomb_quad(p,60,chnks,0.75,1/120);
    else
        psd_p = NaN;
        psd_p_quad = NaN;
        fr_p = NaN;
        fr_pq = NaN;
        disp(sprintf('bad p, %2.0f',inst(kk)))
    end
    
    if sum(isfinite(u)) > sum(isnan(u))*2; % only process spectra if there is enough good data
        [psd_u,fr_u]=mywelch_plomb(u,60,chnks,0.75,1/120);
        [psd_u_quad,fr_uq]=mywelch_plomb_quad(u,60,chnks,0.75,1/120);
    else
        psd_u = NaN;
        psd_u_quad = NaN;
        disp(sprintf('bad u, %2.0f',inst(kk)))
        fr_u = NaN;
        fr_uq = NaN;
    end
    
    if sum(isfinite(v)) > sum(isnan(v))*2; % only process spectra if there is enough good data
        [psd_v,fr_v]=mywelch_plomb(v,60,chnks,0.75,1/120);
        [psd_v_quad,fr_vq]=mywelch_plomb_quad(v,60,chnks,0.75,1/120);
    else
        psd_v = NaN;
        psd_v_quad = NaN;
        fr_v = NaN;
        fr_vq = NaN;
        disp(sprintf('bad v, %2.0f',inst(kk)))
    end
    
    %
    % % fft spectrum, with 127 minute chunks, 75% overlap, linear
    %         % detrending
    %
    %         if sum(isfinite(p)) > sum(isnan(p)); % only process spectra if there is enough good data
    %             Inan = find(isnan(p));
    %             Irl = find(isfinite(p));
    %             p(Inan) = interp1(Irl,p(Irl),Inan,'linear',0);
    %             [psd_p,fr_p]=mywelch(p,60,chnks*size(days,1),0.75);
    %             [psd_p_quad,fr_pq]=mywelch_quad(p,60,chnks*size(days,1),0.75);
    %         else
    %             psd_p = NaN;
    %             psd_p_quad = NaN;
    %             fr_p = NaN;
    %             fr_pq = NaN
    %             disp(sprintf('bad p, %2.0f',inst(kk)))
    %         end
    %
    %         if sum(isfinite(u)) > sum(isnan(u)); % only process spectra if there is enough good data
    %             Inan = find(isnan(u));
    %             Irl = find(isfinite(u));
    %             u(Inan) = interp1(Irl,u(Irl),Inan,'linear',0);
    %             [psd_u,fr_u]=mywelch(u,60,chnks*size(days,1),0.75);
    %             [psd_u_quad,fr_uq]=mywelch_quad(u,60,chnks*size(days,1),0.75);
    %         else
    %             psd_u = NaN;
    %             psd_u_quad = NaN;
    %             disp(sprintf('bad u, %2.0f',inst(kk)))
    %             fr_u = NaN;
    %             fr_uq = NaN
    %         end
    %
    %         if sum(isfinite(v)) > sum(isnan(v)); % only process spectra if there is enough good data
    %             Inan = find(isnan(v));
    %             Irl = find(isfinite(v));
    %             v(Inan) = interp1(Irl,v(Irl),Inan,'linear',0);
    %             [psd_v,fr_v]=mywelch(v,60,chnks*size(days,1),0.75);
    %             [psd_v_quad,fr_vq]=mywelch_quad(v,60,chnks*size(days,1),0.75);
    %         else
    %             psd_v = NaN;
    %             psd_v_quad = NaN;
    %             fr_v = NaN;
    %             fr_vq = NaN;
    %             disp(sprintf('bad v, %2.0f',inst(kk)))
    %         end
    %
    
    
    
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

figure(6); clf
clrs = hsv(length(inst));
leg = {};
vbar = 0;
hp = [];

subplot 211
N = 0;
for kk = 1:length(inst);
    hptmp = plot(PSD(kk).freq_u,PSD(kk).u,'color',clrs(kk,:),'linewidth',1.5);
    hp(kk) = hptmp(1);
    hold on
    leg{kk} = sprintf('#%02.0f, V=%01.2f U=%01.2f',inst(kk),nanmean(PSD(kk).vbar),nanmean(PSD(kk).ubar));
    vbar_tmp = nanmean(PSD(kk).vbar);
    if isfinite(vbar_tmp);
        vbar = vbar + vbar_tmp;
        N = N+1;
    end
end
vbar = vbar/N;
ylabel('U psd (m^2/s)')
xlabel('frequency (Hz)')



legend(hp,leg)
xlim([0 .005])
title(sprintf('RODSEX %s/%s, Vbar = %01.02f m/s',days(1,1:2),days(1,3:4),vbar))

subplot 212
for kk = 1:length(inst);
    plot(PSD(kk).freq_v,PSD(kk).v,'color',clrs(kk,:),'linewidth',1.5)
    hold on
end
ylabel('V psd (m^2/s)')
xlabel('frequency (Hz)')
xlim([0 .005])


orient tall
eval(sprintf('print -dpdf RODSEX_ELF_%s_xshore_70sline.pdf',days(1,:)))



