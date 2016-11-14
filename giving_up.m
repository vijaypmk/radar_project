pd = 0.9;            % Probability of detection
pfa = 1e-6;          % Probability of false alarm
max_range = 150;    % Maximum unambiguous range
range_res = 0.5;      % Required range resolution old - 0.05
tgt_rcs = 0.05;         % Required target radar cross section old - 0.1

prop_speed = physconst('LightSpeed');   % Propagation speed
pulse_bw = prop_speed/(2*range_res);    % Pulse bandwidth
pulse_width = 1/pulse_bw;               % Pulse width
prf = prop_speed/(2*max_range);         % Pulse repetition frequency
fs = 2*pulse_bw;                        % Sampling rate
% waveform = phased.RectangularWaveform(...
%     'PulseWidth',1/pulse_bw,...
%     'PRF',prf,...
%     'SampleRate',fs);

% adc
adc_fs = fs/65536;

waveform = phased.LinearFMWaveform(...
    'PulseWidth',pulse_width,...
    'PRF',prf,...
    'SampleRate',fs,...
    'SweepBandwidth',1e9,...
    'SweepDirection','Up',...
    'Envelope','Rectangular',...
    'OutputFormat','Pulses','NumPulses',1);
% waveform = phased.LinearFMWaveform('PulseWidth',pulse_width,...
% 'SweepBandwidth',1e9,'PRF',prf,'SampleRate',fs);
figure;
disp(waveform.PulseWidth*waveform.SweepBandwidth)
plot(waveform)

noise_bw = pulse_bw;

receiver = phased.ReceiverPreamp(...
    'Gain',20,...
    'NoiseFigure',0,...
    'SampleRate',fs,...
    'EnableInputPort',true);

figure;
snr_db = [-inf, 0, 3, 10, 13];
rocsnr(snr_db,'SignalType','NonfluctuatingNoncoherent');

num_pulse_int = 10;
rocsnr([0 3 5],'SignalType','NonfluctuatingNoncoherent',...
    'NumPulses',num_pulse_int);

snr_min = albersheim(pd, pfa, num_pulse_int)

tx_gain = 20;

fc = 24e9;
lambda = prop_speed/fc;

peak_power = radareqpow(lambda,max_range,snr_min,pulse_width,...
    'RCS',tgt_rcs,'Gain',tx_gain)

transmitter = phased.Transmitter(...
    'Gain',tx_gain,...
    'PeakPower',peak_power,...
    'InUseOutputPort',true);

antenna_t = phased.IsotropicAntennaElement(...
    'FrequencyRange',[5e9 40e9]);

antenna = phased.ULA('Element',antenna_t,'NumElements',10,...
    'ElementSpacing',lambda/2);
antenna.Element.FrequencyRange = [5e9 40e9];


sensorpos = [0; 0; 0];
sensorvel = [0; 0; 0];
sensormotion = phased.Platform(...
    'InitialPosition',sensorpos,...
    'Velocity',sensorvel);

radiator = phased.Radiator(...
    'Sensor',antenna,...
    'OperatingFrequency',fc);

collector = phased.Collector(...
    'Sensor',antenna,...
    'OperatingFrequency',fc);

% tgtpos = [[2024.66;0;0],[3518.63;0;0],[3845.04;0;0]];
tgtpos = [[2.024;0;2],[3.5;0;1],[4.3;0;1.5],[10.75;0;0]];
tgtvel = [[0;0;0],[0;0;0],[0;0;0],[0;0;0]];
tgtmotion = phased.Platform('InitialPosition',tgtpos,'Velocity',tgtvel);

% tgtrcs = [1.6 2.2 1.05];
tgtrcs = [0.6 0.2 0.1 0.8];
target = phased.RadarTarget('MeanRCS',tgtrcs,'OperatingFrequency',fc);

channel = phased.FreeSpace(...
    'SampleRate',fs,...
    'TwoWayPropagation',true,...
    'OperatingFrequency',fc);

fast_time_grid = unigrid(0,1/fs,1/prf,'[)');
slow_time_grid = (0:num_pulse_int-1)/prf;

receiver.SeedSource = 'Property';
receiver.Seed = 2007;

% Pre-allocate array for improved processing speed
rxpulses = zeros(numel(fast_time_grid),num_pulse_int);

% DOA
Nsnapshots = 1000;
rng default
npower = 0.01;

% figure;
for m = 1:num_pulse_int

    % Update sensor and target positions
    [sensorpos,sensorvel] = step(sensormotion, 1/prf);
    [tgtpos,tgtvel] = step(tgtmotion, 1/prf);

    % Calculate the target angles as seen by the sensor
    [tgtrng,tgtang] = rangeangle(tgtpos,sensorpos);

    % Simulate propagation of pulse in direction of targets
    pulse = step(waveform);
    [txsig,txstatus] = step(transmitter, pulse);
    txsig = step(radiator, txsig, tgtang);
    txsig = step(channel,txsig,sensorpos,tgtpos,sensorvel,tgtvel);
%     plot(txsig)

    % Reflect pulse off of targets
    tgtsig = step(target, txsig);

    % Receive target returns at sensor
    rxsig = step(collector,tgtsig,tgtang);
%     rxasig = sensorsig(getElementPosition(antenna)/lambda,...
%     Nsnapshots,tgtang,npower);
%     rxasig = collectPlaneWave(antenna, tgtsig, tgtang);
%     rxisig = pulsint(rxsig,'noncoherent');
    rxisig = sum(rxsig')';
    rxpulses(:,m) = step(receiver,rxisig,~(txstatus>0));
end

rxasig = sensorsig(getElementPosition(antenna)/lambda,...
    Nsnapshots,tgtang,npower);
rxasig = collectPlaneWave(antenna, rxasig, tgtang);
    
oldrxpulses = rxpulses;

npower = noisepow(noise_bw,receiver.NoiseFigure,...
    receiver.ReferenceTemperature);
threshold = npower * db2pow(npwgnthresh(pfa,num_pulse_int,'noncoherent'))

threshold = threshold*500;

figure;
num_pulse_plot = 2;
helperRadarPulsePlot(rxpulses,threshold,...
    fast_time_grid,slow_time_grid,num_pulse_plot);

%check
% rxpulses = pulsint(rxpulses,'noncoherent');

matchingcoeff = getMatchedFilter(waveform);
matchedfilter = phased.MatchedFilter(...
    'Coefficients',matchingcoeff,...
    'GainOutputPort',true);
[rxpulses, mfgain] = step(matchedfilter, rxpulses);

%check
% rxpulses = pulsint(rxpulses,'noncoherent');

matchingdelay = size(matchingcoeff,1)-1;
rxpulses = buffer(rxpulses(matchingdelay+1:end),size(rxpulses,1));

threshold = threshold * db2pow(mfgain);

figure;
helperRadarPulsePlot(rxpulses,threshold,...
    fast_time_grid,slow_time_grid,num_pulse_plot);

range_gates = prop_speed*fast_time_grid/2;

tvg = phased.TimeVaryingGain(...
    'RangeLoss',2*fspl(range_gates,lambda),...
    'ReferenceLoss',2*fspl(max_range,lambda));

rxpulses = step(tvg, rxpulses);

figure;
helperRadarPulsePlot(rxpulses,threshold,...
    fast_time_grid,slow_time_grid,num_pulse_plot);

% copying rxpulse
vxpulses = rxpulses;

rxpulses = pulsint(rxpulses,'noncoherent');

helperRadarPulsePlot(rxpulses,threshold,...
    fast_time_grid,slow_time_grid,1);

[~,range_detect] = findpeaks(abs(rxpulses),'MinPeakHeight',sqrt(threshold));

true_range = tgtrng
range_estimates = range_gates(range_detect)

% [rpeaks,rlocs]=pkpicker( abs(rxpulses(:,1)), threshold, 30);

vd = vxpulses(range_detect,:);         % use only rows where the peaks are located

fd = [];
V = [];
NFFT = 1024;

figure;
for k=1:length(range_detect)
    % using dtft to find doppler frequency
    [H, W] = dtft(vd(k,:),NFFT);  
    
    tempo = 1:length(H);
    plot(tempo, abs(H))
    title('Doppler Frequecny Shift For Different Pulses')
    xlabel('Frequency in MHz')
    ylabel('Amplitude')
    grid on;
    hold on;
    
    % Find location of peak in fft
    [peakH,loc]=max(abs(H));

    % calculate doppler shift and velocity
    fd(k) = W(loc)/((1/prf)*2*pi*fc);
    V(k)=(3e8/2)*(fd(k));
end

V

% Calculate final target range and angle
[FinalRng,FinalAng] = rangeangle(tgtpos,...
    sensormotion.InitialPosition);
DeltaRng = FinalRng-tgtrng;

tgtang

estimator = phased.BeamscanEstimator('SensorArray',antenna,...
    'OperatingFrequency',fc,'ScanAngles',-90:90,...
    'DOAOutputPort',true,'NumSignals',4);
[y,sigang] = step(estimator,(rxasig));
sigang
plotSpectrum(estimator)

MUSICestimator = phased.RootMUSICEstimator('SensorArray',antenna,...
    'OperatingFrequency',fc,'NumSignalsSource','Property',...
    'NumSignals',4,'ForwardBackwardAveraging',true);
doa_est = step(MUSICestimator,(rxasig))

% xy graph
tan_azangles = tan(doa_est);
calcu = tan_azangles.^2 + 1;
calcu = real(sqrt(calcu - range_estimates))


% radar ambiguity
% [afmag_lfm,delay_lfm,doppler_lfm] = ambgfun(pulse,...
% waveform.SampleRate,waveform.PRF);
% surf(delay_lfm*1e6,doppler_lfm/1e3,afmag_lfm,...
% 'LineStyle','none');
% axis tight; grid on; view([140,35]); colorbar;
% xlabel('Delay \tau (\mus)');
% ylabel('Doppler f_d (kHz)');
% title('Linear FM Pulse Waveform Ambiguity Function');

