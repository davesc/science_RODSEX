% find waves from ring pressure data 

% sample period
dt = 1/8;
% record length
Lrec = 24576;

% 1st derivative
d1 = [-1, 0, 1] / (2*dt);
% 2nd derivative
d2 = [1, -2, 1] / (2*dt.^2);
% third derivative 
d3 = [-1, 2, 0, -2, 1] / (2*dt.^3);


% wave "face size" cut off, where only bigger waves are condsidered
wave_size_cutoff = 0.4; %(m)
% size of region (in bins) to search for the bottom or top of a wave
sblock = 8;
% vorticity bins to average before the wave trough
vavg_before = -7:0;
% vorticity bins to average after the wave crest
vavg_after = 0:7;

% 0.5 Hz low pass filter 
cut = 0.5; % cutoff frequency in Hz
Lwin = 128; % cutoff window length
frq = [0:(Lrec/2 - 1)].'/Lrec/dt;
[dump,Icut] = min(abs(frq-cut));
flt = zeros(Lrec/2,1);
flt(frq<=cut) = 1;
win = hanning(2*Lwin);
win = win(Lwin+1:end);
flt((Icut-Lwin/2+1):(Icut+Lwin/2)) = win;
flt = [flt; flipud(flt)];

%TODO: add depth correction prior to filtering
%      maybe cut off spectrum (initially) at the high frequency minimum 

jj = 1;
t = RD(jj).tmat;
p0 = RD(jj).p;
fp = fft(p0);
p = real(ifft(flt.*fp));
pd1 = conv(p,fliplr(d1),'same');
pd2 = conv(p,fliplr(d2),'same');
pd3 = conv(p,fliplr(d3),'same');



figure(1); clf
ax(1) = subplot(211);
plot(t,p,t,pd1,t,pd2,t,pd3/10)
legend('p','pd1','pd2','pd3')

ax(2) = subplot(212);
Iface = (pd1>0 & pd3<0);
ts = t;
ts(~Iface) = NaN;
ps = p;
ps(~Iface) = NaN;
% plot( t,p0,'g',t,p,'b',ts,ps,'r','linewidth',1.5)
plot(t,p,'b',ts,ps,'r','linewidth',1.5)
% plot(t,p,ts,ps,t,p0,'linewidth',1.5)

linkaxes(ax,'x')


Istart = find(diff(Iface)>.5);
Iend = find(diff(Iface)<-.5);
Iend = Iend(Iend>(Istart(1)));
Istart = Istart(Istart<Iend(end));

% TODO: figure out how to estimate wave direction when a wave crosses the
% ring

wvs = []; % mean_time, face_height, change in vorticity, hr, Itrough, Icrest
N = 0;

for kk=1:length(Istart)
    I1 = Istart(kk);
    I2 = Iend(kk);
    [dump,Itmp] = min(p0((I1-sblock):(I1+2)));
    % bottom of the wave
    I11 = I1 - (sblock+1) + Itmp;
    % top of the wave
    [dump,Itmp] = max(p0((I2-2):(I2+(sblock-1))));
    I22 = I2 - 2 + Itmp - 1;
    h = p0(I22)-p0(I11);
    if h >= wave_size_cutoff;
        N=N+1;
        vort_diff = mean(RD(jj).vort(I2+vavg_after))...
                    - mean(RD(jj).vort(I1+vavg_before));
        tmean = mean(t(I11:I22));
        wvs(N,1:6) = [tmean, h, vort_diff,jj,I11,I22];
    end
end
%%

figure(2); clf
plot(t,p0); 
hold on; 
for ii = 1:length(wvs); 
    I = wvs(ii,5):wvs(ii,6);
%     I = Istart(ii):Iend(ii); 
    plot(t(I),p0(I),'r'); 
end


