clear all
close all
clc

fc=24e9;    % center frequency (10 GHz)
%T=500e-6;   % chirp period (500 us)

r_res = 1;
c = 3e8;
B1 = c/(2*r_res);
% B2 = c/(2.5*r_res);
B2=0.75*B1;
r_max=50;         % target range 15 km
overSamp = 4;
fs = B1*overSamp;
Ts=1/fs;        % sampling period (40 ns)

lambda = c/fc;

% creating baseband fmcw waveform
% hwav = phased.FMCWWaveform('SweepBandwidth',B1,...
%     'SampleRate',fs,'SweepDirection','Triangle','NumSweeps',2);
% hwav1 = phased.FMCWWaveform('SweepBandwidth',B1,...
%     'SampleRate',fs);
% s1 = step(hwav1);
% Ls1 = length(s1);
% tau1 = Ts*Ls1;

nB2 = [B1 B2];
hwav2 = phased.FMCWWaveform('SweepBandwidth',nB2,...
    'SampleRate',fs,'SweepDirection','Triangle','NumSweeps',4);
s2 = step(hwav2);
Ls1 = length(s2(1:end/2,:));
Ls2 = length(s2(end/2 + 1 : end,:));
tau1 = Ts*Ls1;
tau2 = Ts*Ls2;

sweep_slope1 = B1/(tau1/2);
sweep_slope2 = B2/(tau2/2);

% s = [s1;s2];
s = s2;

windowlength = 32;
noverlap = 16;
nfft = 32;
mult = 3;
% spectrogram(s2,mult*windowlength,mult*noverlap,mult*nfft,fs,'yaxis');

%target 1
tgt_dist1 = [2;0;0];
tgt_vel1 = [7;0;0];
tgt_rcs1 = [2];
htgt1 = phased.RadarTarget('MeanRCS',tgt_rcs1,'PropagationSpeed',c,...
    'OperatingFrequency',fc);
%platform (target)
htgtplatform1 = phased.Platform('InitialPosition',tgt_dist1,'Velocity',tgt_vel1);

%target 2
tgt_dist2 = [5;0;0];
tgt_vel2 = [10;0;0];
tgt_rcs2 = [2];
htgt2 = phased.RadarTarget('MeanRCS',tgt_rcs2,'PropagationSpeed',c,...
    'OperatingFrequency',fc);
%platform (target)
htgtplatform2 = phased.Platform('InitialPosition',tgt_dist2,'Velocity',tgt_vel2);

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

antenna = phased.IsotropicAntennaElement(...
    'FrequencyRange',[5e9 40e9]);

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
        hradarplatform,hwav2.SweepTime);       % Radar moves during sweep
    [tgt_pos1,tgt_vel1] = step(htgtplatform1,...
        hwav2.SweepTime);                      % Target moves during sweep
    
    [tgt_pos2,tgt_vel2] = step(htgtplatform2,...
        hwav2.SweepTime);                      % Target moves during sweep
    
    [tgtrng1,tgtang1] = rangeangle(tgt_pos1,radar_pos);
    [tgtrng2,tgtang2] = rangeangle(tgt_pos2,radar_pos);
    
%     x = step(hwav1);                           % Generate the FMCW signal
    x = step(hwav2);                           % Generate the FMCW signal
%     x = [x1;x2];
    
    [txsig,txstatus1] = step(htx, x);
    txsig = step(radiator, txsig, tgtang1);
    tx_sig1 = step(hchannel,txsig,radar_pos,tgt_pos1,radar_vel,tgt_vel1);
    tx_sig2 = step(hchannel,txsig,radar_pos,tgt_pos2,radar_vel,tgt_vel2);
    
    xt1 = step(htgt1,tx_sig1);                       % Reflect the signal
    xt2 = step(htgt1,tx_sig2);                       % Reflect the signal
    
    rx_sig1 = step(collector,xt1,tgtang1);
    rx_sig2 = step(collector,xt2,tgtang2);
%     rxsig = xt;
    
    xt1 = step(hrx,rx_sig1);                        % Receive the signal
    xt2 = step(hrx,rx_sig2);                        % Receive the signal
%     nxt(:,m) = xt;
    xd1 = dechirp(xt1,x);                       % Dechirp the signal
    xd2 = dechirp(xt2,x);                       % Dechirp the signal

%     step(hspec,[xt xd]);                      % Visualize the spectrum
%     xd = pulsint(xd);
    xr1(:,m) = xd1;                             % Buffer the dechirped signal
    xr2(:,m) = xd2;                             % Buffer the dechirped signal
end

hrdresp = phased.RangeDopplerResponse('PropagationSpeed',c,...
    'DopplerOutput','Speed','OperatingFrequency',fc,'SampleRate',fs,...
    'RangeMethod','FFT','SweepSlope',sweep_slope1,...
    'RangeFFTLengthSource','Property','RangeFFTLength',2048,...
    'DopplerFFTLengthSource','Property','DopplerFFTLength',256);

% figure;
% clf;
% plotResponse(hrdresp,xr);

% beat frequency
fbu1 = rootmusic(pulsint(xr2(1:length(xr2)/4,1:1:end),'coherent'),1,fs);
fbd1 = rootmusic(pulsint(xr2((length(xr2)/4 + 1):length(xr2)/2,1:1:end),'coherent'),1,fs);

fbu2 = rootmusic(pulsint(xr2(length(xr2)/2 + 1:3*length(xr2)/4,1:1:end),'coherent'),1,fs);
fbd2 = rootmusic(pulsint(xr2((3*length(xr2)/4 + 1):end,1:1:end),'coherent'),1,fs);

% range
% rng_est = beat2range([fbu fbd],sweep_slope,c)
%rng_est= c*(fbu-fbd)*sweep_slope/2

% doppler frequency
% fd = (fbu + fbd);
% 
% % velocity
% vel_est = lambda*fd/2

%velocity
% peak_loc = val2ind(rng_est,c/(fs*2));
% fd = rootmusic(xr(peak_loc,:),1,1/tau);
% % v_est = dop2speed(fd,lambda)
% v_est = lambda*fd

% num_sweeps = 2;
% y = [fbu1; fbd1; fbu2; fbd2];
% % for i = 1:num_sweeps
%     A(1, :) = [2*sweep_slope1/c -2/lambda];
%     A(2, :) = [2*sweep_slope2/c -2/lambda];
%     A(3, :) = [2*sweep_slope1/c -2/lambda];
%     A(4, :) = [2*sweep_slope2/c -2/lambda];
% % end
% x = A\y;


% r_hat1 = c/2*(fbu1-fbu2)*(sweep_slope1-sweep_slope2)^(-1);
% r_hat1 = c/4*(fbu1-fbd2)*(sweep_slope1)^(-1);
r_hat1 = c*(fbu1-fbd1)/(sweep_slope1*4);
% v_hat1 = lambda*(sweep_slope1*4)/(c*r_hat1) - (lambda/2)*fbu1;
% v_hat1 = lambda*(sweep_slope1*4*r_hat1)/(c) - (lambda/2)*fbu1;
v_hat1 = c*(fbu1+fbd1)/(4*(fc+B1));
    
% range
rng_est = beat2range([fbu1 fbd1],sweep_slope1,c);

%velocity
peak_loc = val2ind(rng_est,c/(fs*2));
fd = rootmusic(xr2(peak_loc,:),1,1/(tau1+tau2));
v_hat2 = lambda*fd*2;

sol = [r_hat1 rng_est;v_hat1 v_hat2]

