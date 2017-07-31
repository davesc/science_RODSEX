% rotate ADV data into component tangent to the ring for vorticity

% adv_ang = [Levi measurement, Bill measurement]
% adv_ang = degrees counterclockwise to rotate ADV so x-arm faces outward normal on the ring
adv_ang = [244.0, 245.0
    217.5, 217.5
    193.0, 192.5
    167.0, 167.0
    142.0, 142.0
    116.0, 116.5
    91.8, 91.5
    65.6,66.0
    41.0,41.5
    15.8, 16.5
    349.0, 349.5
    324.2,324.5
    298.0,298.0
    271.0,271.5]
% 01=244.0,244.0  extra measurement for adv 1

% note: when we KVH'd the ring all the ADVs were aligned between 359 and 6
% degrees. Levi thinks the ADVs were closer to 0 degrees and the larger
% numbers were due to the iron pipes thowing off the compas.

adv_ang2 = mean(adv_ang,2)*pi/180;

adv_ang3 = [adv_ang2(11:end); adv_ang2(1:10)];
adv_num3 = [11:14,1:10].';

%%



dth = 360/14*pi/180;

%%

% adv's are close to pointing magnetic north, assume x-arm points to 0 mag

S_rot = 84*pi/180; % steve's rotation to x-cross-shore on Oct 9

Trot = adv_ang2 + S_rot; % total counter clockwise rotation to make ADV's ring-normal

% dir1 = '/Volumes/Data2/RODSEX/data_rough_20131023';
%  dir1 = '~/work/RODSEX/data_rough_20131023';
dir1 = '/Volumes/ThunderBay/RODSEX/data_rough_20131023';


DIR2 = ['0925'
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

hr = (0:23).';

L = 24576;
dt = 1/8;
t = (0:L-1).'*dt;  % time in sec.. not sure if this should start at 0 or 1/8


th = (dth:dth:2*pi).';
xy = [2.5*cos(th), 2.5*sin(th)]; % sensor positions, assuming #14 is offshore


s = 2*2.5*tan(dth/2);  % length of a side on 14 sided polygon outside the circle

array_area = 2.5*s/2*14;

array_area2 = s^2*14/(4*tan(pi/14));

%  rot = [cos(Trot), -sin(Trot)
%         sin(Trot),  cos(Trot)];

RDH = struct;


nn=1;
for ii = 1:size(DIR2,1);
    RD = struct; %structure for each days data
    for jj = 1:length(hr);
        
        
        u = zeros(L,14);
        v = zeros(L,14);
        ct = zeros(L,14); % circle tangent velocity
        cn = zeros(L,14); % circle normal velocity
        load_check = ones(4,1);
        for kk = 1:14
            try
                fid = fopen([dir1,'/',DIR2(ii,:),'/',DIR2(ii,:),sprintf('%02.0f',hr(jj)),'00.u',sprintf('%02.0f%',kk)]);
                u(:,kk) = fread(fid,'float32') / 100;
                fclose(fid);
            catch
                u(:,kk) = NaN;
                %     load_check(1) = 0;
            end
            
            
            try
                fid = fopen([dir1,'/',DIR2(ii,:),'/',DIR2(ii,:),sprintf('%02.0f',hr(jj)),'00.v',sprintf('%02.0f%',kk)]);
                v(:,kk) = fread(fid,'float32') / 100;
                fclose(fid);
            catch
                v(:,kk) = NaN;
                %     load_check(2) = 0;
            end
            
            
            
            
            
            try
                ct(:,kk) =  v(:,kk)*cos(Trot(kk)) - u(:,kk)*sin(Trot(kk));
                cn(:,kk) = u(:,kk)*cos(Trot(kk)) + v(:,kk)*sin(Trot(kk));
            catch
                ct(:,kk) = NaN;
                cn(:,kk) = NaN;
                %     load_check(3) = 0;
            end
            
        end
        
        
        I = find(isnan(ct(1,:)));
        for mm = 1:length(I);
            if I(mm)>1 & I(mm)<14;
                ct(:,I(mm)) = (ct(:,I(mm)-1)+ct(:,I(mm)+1))/2;
            elseif I(mm)==1 ;
                ct(:,I(mm)) = (ct(:,14)+ct(:,2))/2;
            elseif I(mm)==14;
                ct(:,I(mm)) = (ct(:,13)+ct(:,1))/2;
            end
        end
        
        % ct(isnan
        
        
        
        
        p = zeros(L,2);
        for kk = 1:2
            try
                fid = fopen([dir1,'/',DIR2(ii,:),'/',DIR2(ii,:),sprintf('%02.0f',hr(jj)),'00.q',sprintf('%02.0f%',kk)]);
                p(:,kk) = fread(fid,'float32') / 100;  % pressure in meters depth
                fclose(fid);
            catch
                p(:,kk) = NaN;
                %     load_check(4) = 0;
            end
        end
        
        
        
        
        chnks = 16;
        pover = .5;
        dt = 1/8;
        % doffu = 0.75; % this is doff from the ring, but needs to be corrected for the bed height
        % doffp = 0.40; % ????? not sure what doffp is
        %
        doffu = 0.35; %height above bed, when bed is at the top of the pressure sensor
        doffp = 0; % say pressure sensor is at the bed
        
        if sum(isnan([p(:,2);u(:,10);v(:,10)]))==0;
            %  USING ADV 10 AND P2
            [fm,SEee,SEUU,Ztest,direc,spread,mdirec,mspread,Sxx,Syy,Sxy,ma1,mb1,ma2,mb2]=get_spectra_RODSEX_fixed_corrections(p(:,2),u(:,10),v(:,10),dt,chnks,pover,doffu,doffp);
            [Urot, R, U_igvar, V_igvar, E_igvar, depth] = vortex_lippmann_ratio(u(:,10),v(:,10), p(:,2), mean(p(:,2)), dt);
        else
            fm=[NaN; NaN];
            SEee=NaN;
            SEUU=NaN;
            Ztest=NaN;
            direc=NaN;
            spread=NaN;
            mdirec=NaN;
            mspread=NaN;
            Sxx=NaN;
            Syy=NaN;
            Sxy=NaN;
            ma1=NaN;
            mb1=NaN;
            ma2=NaN;
            mb2=NaN;
            Urot=NaN;
            R=NaN;
        end
        
        fm = fm.';
        df = fm(2)-fm(1);
        ss = find((fm>=(1/20)) & (fm<=(1/4)));  % sea swell band
        HsigEE = 4*sqrt(sum(SEee(ss)).*df);   % Hsig from pressure
        HsigUU = 4*sqrt(sum(SEUU(ss)).*df);   % Hsig from u and v velocities
        
        [dump, Ip] = max(SEUU);
        Tp = 1./fm(Ip);  % peak period
        
        Tm = 1/(nansum(SEUU.*fm*df)/nansum(SEUU*df));  % spectal mean period
        Sxy_mean = nansum(SEUU.*Sxy*df)/nansum(SEUU*df);  % spectal mean Sxy
        
        vort = - (1/array_area * sum(ct,2)*s); % vorticity (1/s)   ????? I put in a negative sign, but not sure it's right ????
        mdiv = 1/array_area * sum(cn,2)*s; % mean divergence (1/s)
        [Svort,fm]=mywelch(vort,dt,chnks,pover);
        [Sdiv,fm]=mywelch(sum(cn,2),dt,chnks,pover);
        
        % calc potential vorticity
        pvort = vort./(p(:,2)-doffp); % using doffp = 0 is not quite right, so there may be long term errors as the bed elevation changed
        [Spvort,fm]=mywelch(pvort,dt,chnks,pover);
        
        % simple u and v spectra, using combined data from 5, 7, and 10
        SUUtmp = mywelch(u(:,[5,7,10]),dt,chnks,pover);
        SUU = sum(SUUtmp,2)/3;
        SVVtmp = mywelch(v(:,[5,7,10]),dt,chnks,pover);
        SVV = sum(SVVtmp,2)/3;
        
        
        RDH.U(nn,1:14) = mean(u);
        RDH.V(nn,1:14) = mean(v);
        RDH.stdV(nn,1:14) = std(v);
        RDH.stdU(nn,1:14) = std(u);
        RDH.mvort(nn,1) = mean(vort);
        RDH.stdvort(nn,1) = std(vort);
        RDH.mpvort(nn,1) = mean(pvort);
        RDH.stdpvort(nn,1) = std(pvort);
        % RDH(jj).mdiv = mdiv;
        RDH.t_hour(nn,1) = mean(t + hr(jj)*3600)/24; % mean time in hours of the day
        RDH.tmat(nn,1) = mean(datenum(2013,str2num(DIR2(ii,1:2)),str2num(DIR2(ii,3:4)),hr(jj),0,t));
        RDH.P(nn,1:2) = mean(p);
        
        RDH.mdirec(nn,1) = -mdirec;
        RDH.mspread(nn,1:2) = mspread;
        RDH.HsigEE(nn,1)  = HsigEE;
        RDH.HsigUU(nn,1) = HsigUU;
        RDH.Tp(nn,1) = Tp;
        RDH.Tm(nn,1) = Tm;
        RDH.Sxy(nn,1)  = Sxy_mean;
        RDH.Urot(nn,1) = Urot; %lippmann Urot
        RDH.R(nn,1) = R; % lippmann R
        
        nn = nn+1;
        
        RD(jj).t_sec = t + hr(jj)*3600;
        RD(jj).tmat = datenum(2013,str2num(DIR2(ii,1:2)),str2num(DIR2(ii,3:4)),hr(jj),0,t);
        RD(jj).vort = vort;
        RD(jj).pvort = pvort;
        RD(jj).mdiv = mdiv;
        RD(jj).fm = fm;
        RD(jj).SEUU = SEUU;
        RD(jj).SEee = SEee;
        RD(jj).SUU = SUU;
        RD(jj).SVV = SVV;
        RD(jj).Svort = Svort;
        RD(jj).Spvort = Spvort;
        RD(jj).Sdiv = Sdiv;
        RD(jj).v = v(:,5);
        RD(jj).u = u(:,5);
        RD(jj).p = p(:,2);
        
        
        

%         save(sprintf('/Volumes/Data2/RODSEX/data_rough_20131023/day_vort/RD_2013%s.mat',DIR2(ii,:)),'RD')
%         save(sprintf('~/work/RODSEX/data_rough_20131023/day_vort/RD_2013%s.mat',DIR2(ii,:)),'RD')
        save(sprintf('/Volumes/ThunderBay/RODSEX/data_rough_20131023/day_vort/RD_2013%s.mat',DIR2(ii,:)),'RD')
        
    end
end

%  save ~/work/RODSEX/data_rough_20131023/RD_hourly.mat RDH
% save /Volumes/Data2/RODSEX/data_rough_20131023/RD_hourly.mat RDH
save /Volumes/ThunderBay/RODSEX/data_rough_20131023/RD_hourly.mat RDH

%%


%
% figure(3); clf
% plot(t,vort)
% xlabel('time (s)')
% ylabel('vorticity (1/s)')



%
%  figure(1); clf
%  t = (0:length(u)-1).'/5;
%  plot(t,u,t,v)


