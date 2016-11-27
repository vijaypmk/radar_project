clear all
close all
clc

fc=24e9;    % center frequency (10 GHz)
%T=500e-6;   % chirp period (500 us)

r_res = 1;
c = 3e8;
B = c/(2*r_res);

r_max=50;         % target range 15 km
overSamp = 4;
fs = B*overSamp;
Ts=1/fs;        % sampling period (40 ns)

lambda = c/fc;

% creating baseband fmcw waveform
hwav = phased.FMCWWaveform('SweepBandwidth',B,...
    'SampleRate',fs,'SweepDirection','Triangle','NumSweeps',2);
s = step(hwav);
Ls = length(s);
tau = Ts*Ls;

sweep_slope = B/(tau/2);

windowlength = 32;
noverlap = 16;
nfft = 32;
mult = 3;
%spectrogram(s,mult*windowlength,mult*noverlap,mult*nfft,fs,'yaxis');


%target
tgt_dist = [2;0;0];
tgt_vel = [0;0;0];
tgt_rcs = [2];
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
    nxt(:,m) = xt;
    xd = dechirp(xt,x);                       % Dechirp the signal

    step(hspec,[xt xd]);                      % Visualize the spectrum

    xr(:,m) = xd;                             % Buffer the dechirped signal
end

hrdresp = phased.RangeDopplerResponse('PropagationSpeed',c,...
    'DopplerOutput','Speed','OperatingFrequency',fc,'SampleRate',fs,...
    'RangeMethod','FFT','SweepSlope',sweep_slope,...
    'RangeFFTLengthSource','Property','RangeFFTLength',2048,...
    'DopplerFFTLengthSource','Property','DopplerFFTLength',256);

figure;
clf;
plotResponse(hrdresp,xr);

% beat frequency
fbu = rootmusic(pulsint(xr(1:length(xr)/2,1:1:end),'coherent'),1,fs);
fbd = rootmusic(pulsint(xr((length(xr)/2 + 1):end,1:1:end),'coherent'),1,fs);

% range
% rng_est = beat2range([fbu fbd],sweep_slope,c)
rng_est= c*(fbu-fbd)*sweep_slope/2

% doppler frequency
% fd = (fbu + fbd);
% 
% % velocity
% vel_est = lambda*fd/2

%velocity
peak_loc = val2ind(rng_est,c/(fs*2));
fd = rootmusic(xr(peak_loc,:),1,1/tau);
% v_est = dop2speed(fd,lambda)
v_est = lambda*fd



