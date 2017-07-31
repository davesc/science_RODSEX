%% RODSEX ELF spectra sum, using lomb scargle , all days, all ADVs


data_dir = '/Volumes/OWC HD/RODSEX/ring_and_standalone_data';
% data_dir = '/Volumes/DAVID_CLARK/RODSEX';

% maximum frequency included in total ELF power
fmax = 0.002;

% structure to save data
ELF = struct;

all_days = ['0925'
    '0926'
    '0927'
    '0928'
    '0929'
    '0930'
    '1001'
    '1002'
    '1003'
    '1004'
    '1005'
    '1006'
    '1007'
    '1008'
    '1009'
    '1010'
    '1011'
    '1012'
    '1013'
    '1014'
    '1015'
    '1016'
    '1017'
    '1018'
    '1019'
    '1020'
    '1021'
    '1022'
    '1023'];


start_hour = 0;
end_hour = 23;

%     % length of record in hours
%     if size(days,1) == 1;
%         rec_hours = end_hour - start_hour + 1;
%     elseif size(days,1) > 1;
%         rec_hours = (size(days,1)-1)*24 - start_hour + 1 + end_hour
%     else
%         rec_hours = NaN;
%     end

% in this case I use 24 hours from each day
rec_hours = 24;


% all the instruments
inst = [21
    22
    23
    24
    27
    28
    29
    31
    32
    33
    41
    42
    43
    51
    52
    05
    10
    55
    54
    56
    61
    62
    63
    71
    72
    73
    74
    93];


Lrec = 6120; %(51 minutes of 2 Hz data), the actual record is 51.2 minutes
Lrec2  = 60; % 60 minutes of 60 averaged data, to be filled / interpolated
Lrec8 = 24480; % 51 min of 8 Hz data


for aa = 1:size(all_days,1);
    
    
    days = all_days(aa,1:4);
 
    
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
        
        
        
        ELF.date(aa,kk) = {days};
        ELF.tmat(aa,kk) = datenum(2013,str2num(days(1:2)),str2num(days(3:4)),12,0,0); % noon on the day
        ELF.inst(aa,kk) = inst(kk);
        
        if isfinite(fr_p(1))
        dfp = fr_p(3)-fr_p(2);
        IELFp = find(fr_p<=fmax);
        ELF.p_power(aa,kk) = sum(psd_p(IELFp)*dfp);
        else
            ELF.p_power(aa,kk) = NaN;
        end
        
        if isfinite(fr_u(1))
        dfu = fr_u(3)-fr_u(2);
        IELFu = find(fr_u<=fmax);
        ELF.u_power(aa,kk) = sum(psd_u(IELFu)*dfu);
        else
            ELF.u_power(aa,kk) = NaN;
        end
        
        if isfinite(fr_v(1))
        dfv = fr_v(3)-fr_v(2);
        IELFv = find(fr_v<=fmax);
        ELF.v_power(aa,kk) = sum(psd_v(IELFv)*dfv);
        else
        ELF.v_power(aa,kk)  = NaN
        end
        
        ELF.pbar(aa,kk) = nanmean(pbar);
        ELF.ubar(aa,kk) = nanmean(ubar);
        ELF.vbar(aa,kk) = nanmean(vbar);

%
        %               if kk==2;
        %                 return
        %             end
    end
    
    
    
    
end


save RODSEX_ELF_power_all.mat ELF