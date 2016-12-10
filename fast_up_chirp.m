% University of California, Santa Cruz
%
% FMCW Radar for Multi Target Range and Velocity Estimation
%---------------------------------------------------------------
% Author - Vijay Muthukumaran
% Course - EE-288 Radar, SAR and ISAR
%---------------------------------------------------------------

clear all
close all
clc
%% System Parameters
fc=24e9;            % center frequency (10 GHz)

r_res = 1;          % resolution 1 m
c = 3e8;            
B = 2*c/(2*r_res);  

r_max=50;           % maximum unambiguous range
overSamp = 4;       % oversampling amount
fs = B*overSamp;    % sampling frequecny
Ts=1/fs;            % sampling time

lambda = c/fc;      % wavelength

%% Target Specifications

graphin = input('Do you want graphical input? (1/0) ');
N = input('How many targets (max 7) would you like to input? ');

% graphical input
if graphin
% Initialize
  ranges = [];
  vels = [];

% Set up window used for input
  figure;
  clf
  grid on
  hold on
  
  axis([0 20 -20 20])
  title(['Please enter up to ' num2str(N) ' targets at min 0.5 m apart using mouse'])
  xlabel('Range (m)')
  ylabel('Velocity (m/s)')

% Loop for getting upto N target points 
  xa = 0;
  ya = 0;
  for k=1:N
    [xa,ya] = ginput(1);
    plot(xa,ya,'o')
    ranges(k,:) = [[xa;0;0]];
    vels(k,:) = [[ya;0;0]];
  end

% othersize get input from keyboard
else
  ranges = input('Enter range vector for up to targets: ');
  vels = input('Enter correponsing velocity vector: ');
  if length(ranges) ~= length(vels)
    error('You need the same number of velocities and ranges')
  end
end

%target
tgt_dist = ranges';
tgt_vel = vels';
tgt_rcs = 2*ones(N,1)';
htgt = phased.RadarTarget('MeanRCS',tgt_rcs,'PropagationSpeed',c,...
    'OperatingFrequency',fc);
%platform (target)
htgtplatform = phased.Platform('InitialPosition',tgt_dist,'Velocity',tgt_vel);

%propogation
hchannel = phased.FreeSpace('PropagationSpeed',c,...
    'OperatingFrequency',fc,'SampleRate',fs,'TwoWayPropagation',true);

%% creating baseband fmcw waveform
hwav = phased.FMCWWaveform('SweepBandwidth',B,...
    'SampleRate',fs,'SweepDirection','Up','NumSweeps',1);
s = step(hwav);
Ls = length(s);
tau = Ts*Ls;
sweep_slope = B/(tau);

figure;
windowlength = 32;
noverlap = 16;
nfft = 32;
mult = 3;
spectrogram(s,mult*windowlength,mult*noverlap,mult*nfft,fs,'yaxis');
title('FMCW signal spectrogram');

%% Antenna Specifications
%radar system
ant_aperture = 2.5000e-05;                      % in square meter
ant_gain = aperture2gain(ant_aperture,lambda);  % in dB

tx_ppower = db2pow(11)*1e-3;                    % in watts
tx_gain = 0+ant_gain;                           % in dB

rx_gain = 26+ant_gain;                          % in dB
rx_nf = 12;                                     % in dB

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
radar_vel = [0;0;0];                           % in m/s
radar_pos = [0;0;0];                           % in m
hradarplatform = phased.Platform('InitialPosition',radar_pos,...
    'Velocity',radar_vel);

%% Simulation

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
    
    xt = step(htgt,txsig);                    % Reflect the signal
    
    rxsig = step(collector,xt,tgtang);
    
    xt = step(hrx,rxsig);                     % Receive the signal
    nxt(:,m) = xt;
    xd = dechirp(xt,x);                       % Dechirp the signal

    xr(:,m) = xd;                             % Buffer the dechirped signal
end

%% Signal Processing
[xx yy] = size(xr');

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

%% Single Target Method

%Extract Range
rangeEst = c*f_p*tau/2/B

%Extract Velocity
samp_extract = xr(bb,:);

N2 = 2^(nextpow2(8*xx));
opt2 = 0;

[f2 X2] = fftMAG(samp_extract,1/tau,N2,opt2);

[aa2 bb2] = max(abs(X2));

f_d = f2(bb2);

velocityEst = c*f_d/2/fc

%% Do first fft to siphon out beat frequency

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

%% Extract only positive frequencies and < Rmax values of matrix

Ypos = Y(:,(N/2 +1):end);
rangeVals1 = c*f((N/2 +1):end)*tau/2/B;

%Find where max range is
rangeCheck = rangeVals1 > r_max;
[val pos] = max(rangeCheck);
Yposmax = Ypos(:,1:(pos+1));
rangeVals2 = rangeVals1(1:(pos+1));
    
%% Do second fft to find doppler frequency

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