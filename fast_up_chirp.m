clear all
close all
clc

fc=24e9;    % center frequency (10 GHz)
%T=500e-6;   % chirp period (500 us)

r_res = 1;
c = 3e8;
B = 2*c/(2*r_res);  % no 2

r_max=50;         % target range 15 km
% tau = 5.5*range2time(r_max,c);
overSamp = 4;
fs = B*overSamp;
% fs = 150e+6;
Ts=1/fs;        % sampling period (40 ns)

lambda = c/fc;

% creating baseband fmcw waveform
hwav = phased.FMCWWaveform('SweepBandwidth',B,...
    'SampleRate',fs,'SweepDirection','Up','NumSweeps',1);
% hwav = phased.FMCWWaveform('SweepTime',tau,'SweepBandwidth',B,...
%     'SampleRate',fs);
s = step(hwav);
Ls = length(s);
tau = Ts*Ls;
sweep_slope = B/(tau);

windowlength = 32;
noverlap = 16;
nfft = 32;
mult = 3;
spectrogram(s,mult*windowlength,mult*noverlap,mult*nfft,fs,'yaxis');
title('FMCW signal spectrogram');


%target
tgt_dist = [[13;0;0],[2;2;0],[5;7;0]];
tgt_vel = [[10;0;0],[0;0;0],[5;0;0]];
tgt_rcs = [2 2 2];
htgt = phased.RadarTarget('MeanRCS',tgt_rcs,'PropagationSpeed',c,...
    'OperatingFrequency',fc);
%platform (target)
htgtplatform = phased.Platform('InitialPosition',tgt_dist,'Velocity',tgt_vel);

%propogation
hchannel = phased.FreeSpace('PropagationSpeed',c,...
    'OperatingFrequency',fc,'SampleRate',fs,'TwoWayPropagation',true);

%radar system
ant_aperture = 2.5000e-05;                         % in square meter
ant_gain = aperture2gain(ant_aperture,lambda);  % in dB

tx_ppower = db2pow(11)*1e-3;                     % in watts
tx_gain = 0+ant_gain;                           % in dB

rx_gain = 26+ant_gain;                          % in dB
rx_nf = 12;                                    % in dB

htx = phased.Transmitter('PeakPower',tx_ppower,'Gain',tx_gain,...
    'InUseOutputPort',true);
hrx = phased.ReceiverPreamp('Gain',rx_gain,'NoiseFigure',rx_nf,...
    'SampleRate',fs);

% antenna_t = phased.IsotropicAntennaElement(...
%     'FrequencyRange',[5e9 40e9]);

% antenna = phased.ULA('Element',antenna_t,'NumElements',10,...
%     'ElementSpacing',lambda/2);
% antenna.Element.FrequencyRange = [5e9 40e9];

antenna = phased.ULA('NumElements',4,...
    'ElementSpacing',lambda/2);
antenna.Element.FrequencyRange = [5e9 40e9];

radiator = phased.Radiator(...
    'Sensor',antenna,...
    'OperatingFrequency',fc);

collector = phased.Collector(...
    'Sensor',antenna,...
    'OperatingFrequency',fc);

% radar system specs
radar_vel = [0;0;0];
radar_pos = [0;0;0];
hradarplatform = phased.Platform('InitialPosition',radar_pos,...
    'Velocity',radar_vel);

hspec = dsp.SpectrumAnalyzer('SampleRate',fs,...
    'PlotAsTwoSidedSpectrum',true,...
    'Title','Spectrum for received and dechirped signal',...
    'ShowLegend',true);

%simulation loop
Nsweep = 16;

for m = 1:Nsweep
    [radar_pos,radar_vel] = step(...
        hradarplatform,hwav.SweepTime);       % Radar moves during sweep
    [tgt_pos,tgt_vel] = step(htgtplatform,...
        hwav.SweepTime);                      % Target moves during sweep
    
    [tgtrng,tgtang] = rangeangle(tgt_pos,radar_pos);
    
    x = step(hwav);                           % Generate the FMCW signal
    
    [txsig,txstatus] = step(htx, x);
    txsig = step(radiator, txsig, tgtang);
    txsig = step(hchannel,txsig,radar_pos,tgt_pos,radar_vel,tgt_vel);
    
    xt = step(htgt,txsig);                       % Reflect the signal
    
    rxsig = step(collector,xt,tgtang);
    
    xt = step(hrx,rxsig);                        % Receive the signal
    xt = pulsint(xt,'coherent');
%     xt = sum(xt,2);
    nxt(:,m) = xt;
    xd = dechirp(xt,x);                       % Dechirp the signal

%     step(hspec,[xt xd]);                      % Visualize the spectrum

    xr(:,m) = xd;                             % Buffer the dechirped signal
end

hrdresp = phased.RangeDopplerResponse('PropagationSpeed',c,...
    'DopplerOutput','Speed','OperatingFrequency',fc,'SampleRate',fs,...
    'RangeMethod','FFT','SweepSlope',sweep_slope,...
    'RangeFFTLengthSource','Property','RangeFFTLength',2048,...
    'DopplerFFTLengthSource','Property','DopplerFFTLength',256);

% stretchproc = getStretchProcessor(hwav,c);

% figure;
% clf;
% plotResponse(hrdresp,xr);

% beat frequency
fb = rootmusic(pulsint(xr,'coherent'),1,fs);
% fbd = rootmusic(pulsint(xr((length(xr)/2 + 1):end,1:1:end),'coherent'),1,fs);

% range
rng_est1 = beat2range(fb,sweep_slope,c)
% rng_est2= c*(fb-fbd)/sweep_slope*4

% doppler frequency
% fd = (fbu + fbd);
% 
% % velocity
% vel_est = lambda*fd/2

%velocity
% peak_loc = val2ind(rng_est1,c/(fs*2));
% fd = rootmusic(xr(peak_loc,:),1,1/tau);
% % fd = rootmusic(xr(peak_loc,:),1,fs);
% v_est = dop2speed(fd,lambda)
% v_est = lambda*fd

%try
[xx yy] = size(xr');
% y = transpose(xr);
% powDopSave = zeros(yy,2);
% for nn = 1:yy
%    [fd pow] = rootmusic(y(:,nn),1,1/(tau));
% %     [fd pow] = rootmusic(y(:,nn),1,fs);
%     powDopSave(nn,1) = fd;
%     powDopSave(nn,2) = pow;    
% end
% 
% [val nn_est] = max(powDopSave(:,2));
% range1 = c/2*Ts*nn_est
% fd_est = powDopSave(nn_est,1);
% velocity1 = dop2speed(fd,lambda)/2

% MUSICestimator = phased.RootMUSICEstimator('SensorArray',antenna,...
%     'OperatingFrequency',fc,'NumSignalsSource','Property',...
%     'NumSignals',length(tgt_rcs),'ForwardBackwardAveraging',true);
% doa_est = step(MUSICestimator,(rxsig));
% 
% doa_est = sort(doa_est)
% 
% % arranging input/output in order of arrival for xy calculations
% [tgtang(1,:), sortaindex] = sort(tgtang(1,:));
% tgtang(1,:)
% rng_est1 = tgtrng;
% rng_est1 = rng_est1(sortaindex)
% 
% % xy graph
% tan_azangles = tand(doa_est);
% x_axis = tan_azangles.^2 + 1;
% x_axis = rng_est1./sqrt(x_axis)
% y_axis = x_axis.*tan_azangles
% 
% grid_xlen = 20;
% grid_ylen = 20;
% figure;
% plot(x_axis, y_axis,'x')
% axis([0 grid_xlen -grid_ylen grid_ylen])
% grid on
% title('Location of Targets With Respect To The Radar at (0, 0)')
% xlabel('x axis(m)')
% ylabel('y axis(m)')

% pats
%Do coherent pulse integration over pulses

y_int = sum(xr,2);


%Find peak frequency of FT of integrated output

N = 2^(nextpow2(yy));

opt = 0;

%Extract Using FFT

[f X] = fftMAG(y_int,fs,N,opt);

%Extract Using Root MUSIC

fMusic = 0:4000:fs/2;
K = 1;
[S,F] = pmusic(y_int,K,fMusic,fs); 

%Find Peaks

thresh  = 50;

[pows locs] = pkpicker(S,thresh,K,'sort');

[aa bb] = max(abs(X));

f_p = f(bb);

f_p2 = F(locs);

%Single Target Method

%Extract Range

rangeEst = c*f_p*tau/2/B

%     rangeEst2 = c*f_p2*tau/2/B

 %Extract Velocity

samp_extract = xr(bb,:);

N2 = 2^(nextpow2(8*xx));
opt2 = 0;

[f2 X2] = fftMAG(samp_extract,1/tau,N2,opt2);

[aa2 bb2] = max(abs(X2));

f_d = f2(bb2);

velocityEst = c*f_d/2/fc


%Do first fft to siphon out beat frequency

y = transpose(xr);

N = 2^(nextpow2(length(y(1,:))));

Y = zeros(xx,N);

for nn = 1:xx

    y_iter = y(nn,:);

    opt = 0;

    %Use fft

    [f X] = fftMAG(y_iter,fs,N,opt);

    Y(nn,:) = X;

end


%Extract only positive frequencies and < Rmax values of matrix

Ypos = Y(:,(N/2 +1):end);

rangeVals1 = c*f((N/2 +1):end)*tau/2/B;

%Find where max range is

rangeCheck = rangeVals1 > r_max;

[val pos] = max(rangeCheck);

Yposmax = Ypos(:,1:(pos+1));

rangeVals2 = rangeVals1(1:(pos+1));
    

%Do second fft to find dopplar frequency

N = 2^(nextpow2(pos+1));

opt = 0;

Yfinal = zeros(N,pos+1);

for nn = 1:pos+1

    cur = Yposmax(:,nn);

    [f X] = fftMAG(cur,1/tau,N,opt);
    Yfinal(:,nn) = X;
end

velocityVals = c*f/2/fc;

figure;

imagesc(rangeVals2,velocityVals,20*log10(abs(Yfinal)))

xlabel('Ranges (m)'),ylabel('Velocities')