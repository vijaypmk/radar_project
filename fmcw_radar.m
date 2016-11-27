clear

fc=24e9;    % center frequency (10 GHz)
%T=500e-6;   % chirp period (500 us)

r_res = 1;
c = 3e8;
B = c/(2*r_res);
%alpha=B/T;  % chirp rate (100 GHz/s)

r_max=50;         % target range 15 km
tau=2*r_max/c;      % target transit time 100 us
tau = tau*5.5;    % why?
%fb=alpha*tau;   % beat frequency 10 MHz

fs = 150e+6;
Ts=1/fs;        % sampling period (40 ns)
%N=T/Ts;         % number of samples per sweep (12,500)
%Np=(T-tau)/Ts;  % number of processed samples per sweep (10,000)
%A=alpha/fs^2;   % dimensionless chirp parameter

lambda = c/fc;
sweep_slope = B/tau;


% hwav = phased.FMCWWaveform('SweepTime',tau,'SweepBandwidth',B,...
%     'SampleRate',fs);
% creating baseband fmcw waveform
hwav = phased.FMCWWaveform('SweepTime',tau,'SweepBandwidth',B,...
    'SampleRate',fs,'SweepDirection','Triangle','NumSweeps',2);
s = step(hwav);
% subplot(211); plot(0:1/fs:tau-1/fs,real(s));
subplot(211); plot(0:1/fs:2*tau-1/fs,real(s));
xlabel('Time (s)'); ylabel('Amplitude (v)');
title('FMCW signal'); axis tight;
subplot(212); spectrogram(s,32,16,32,fs,'yaxis');
title('FMCW signal spectrogram');

% %object
car_dist = [[8;0;0]];  %%here 0.1
car_speed = [[99;0;0]];  %0
% car_rcs = db2pow(min(10*log10(car_dist)+5,20));
car_rcs = [2];

hcar = phased.RadarTarget('MeanRCS',car_rcs,'PropagationSpeed',c,...
    'OperatingFrequency',fc);
hcarplatform = phased.Platform('InitialPosition',car_dist,'Velocity',car_speed);
hchannel = phased.FreeSpace('PropagationSpeed',c,...
    'OperatingFrequency',fc,'SampleRate',fs,'TwoWayPropagation',true);


%radar system
ant_aperture = 2.5000e-05; %6.06e-4; old                         % in square meter
ant_gain = aperture2gain(ant_aperture,lambda);  % in dB

tx_ppower = db2pow(5)*1e-3;                     % in watts
tx_gain = 9+ant_gain;                           % in dB

rx_gain = 15+ant_gain;                          % in dB
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

radar_speed = [0;0;0]; %%here 0
radar_pos = [0;0;0];
hradarplatform = phased.Platform('InitialPosition',radar_pos,...
    'Velocity',radar_speed);

hspec = dsp.SpectrumAnalyzer('SampleRate',fs,...
    'PlotAsTwoSidedSpectrum',true,...
    'Title','Spectrum for received and dechirped signal',...
    'ShowLegend',true);

%simulation loop
% Nsweep = 64;
Nsweep = 16;
% xr = complex(zeros(hwav.SampleRate*hwav.SweepTime,Nsweep));

for m = 1:Nsweep
    [radar_pos,radar_vel] = step(...
        hradarplatform,hwav.SweepTime);       % Radar moves during sweep
    [tgt_pos,tgt_vel] = step(hcarplatform,...
        hwav.SweepTime);                      % Car moves during sweep
    
    [tgtrng,tgtang] = rangeangle(tgt_pos,radar_pos);
    
    x = step(hwav);                           % Generate the FMCW signal
    
    [txsig,txstatus] = step(htx, x);
    txsig = step(radiator, txsig, tgtang);
    txsig = step(hchannel,txsig,radar_pos,tgt_pos,radar_vel,tgt_vel);
    
    xt = step(hcar,txsig);                       % Reflect the signal
    
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
% 
% fm = 1/tau;
% % n=0:N-1;                    % time index
% % t=(-N/2+n)*Ts;              % time grid
% % phib=@(t) phiTX(t)-phiTX(t-tau);
% % R = c/(2*B*fm)
% 
% % [Pyy,F] = periodogram(xt,[],1024,fs,'centered');
% 
% % fb1 = rootmusic(pulsint(xr,'coherent'),1,fs);
% % fb2 = rootmusic(pulsint(xr2,'coherent'),1,fs);
% % or should I use the given range
% % rng_est = beat2range(fb_rng,sweep_slope,c)
% % rng_est1 = c*fb1/(2*sweep_slope)
% % rng_est2 = c*fb2/(2*sweep_slope)
% 
% hwavtr = clone(hwav);
% release(hwavtr);
% % hwavtr.SweepDirection = 'Triangle';
% 
% % xr = helperFMCWSimulate(Nsweep,hwavtr,hradarplatform,hcarplatform,...
% %     htx,hchannel,hcar,hrx);
% 
fbu_rng = rootmusic(pulsint(xr(1:275,1:1:end),'coherent'),1,fs);
fbd_rng = rootmusic(pulsint(xr(276:end,2:1:end),'coherent'),1,fs);

rng_est = beat2range([fbu_rng fbd_rng],sweep_slope,c)