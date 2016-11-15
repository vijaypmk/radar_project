clear

% specifications
pd = 0.9;            % Probability of detection
pfa = 1e-6;          % Probability of false alarm
max_range = 150;    % Maximum unambiguous range
% range_res = 0.5;      % Required range resolution old - 0.05
range_res = 0.025;       % 0.1 seems to work well
tgt_rcs = 0.01;         % Required target radar cross section old - 0.1

prop_speed = physconst('LightSpeed');   % Propagation speed
pulse_bw = prop_speed/(2*range_res);    % Pulse bandwidth
pulse_width = 1/pulse_bw;               % Pulse width
prf = prop_speed/(2*max_range);         % Pulse repetition frequency
fs = 2*pulse_bw;                        % Sampling rate

tx_gain = 20;

fc = 24e9;
lambda = prop_speed/fc;


% Cartesian Space
sensorpos = [0; 0; 0];
sensorvel = [0; 0; 0];

% tgtpos = [[2024.66;0;0],[3518.63;0;0],[3845.04;0;0]];
% tgtpos = [[2.024;-1.5;0],[3.5;0;0],[4.5;-2;0],[7.5;2;0],[1;1;0]];
% tgtvel = [[30;0;0],[20;0;0],[-10;0;0],[0;0;0],[0;0;0]];
% tgtpos = [[2;7;0],[1;1;0],[0.5;-1;0]];
% tgtvel = [[0;0;0],[0;0;0],[0;0;0]];
% temp
% tgtpos = [[2.1455;1.8925;0],[3.8029;4.3224;0],[6.1418;-1.7523;0]];
% tgtrcs = [0.6 0.2 0.4 0.2 0.2];
% tgtrcs = [0.06 0.04 0.05];

graphin = input('Do you want graphical input? ');
N = input('How many targets would you like to input? ');

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
%   mr = maxrange*Rmul;
%   mv = maxvelocity*Vmul;
%   plot([minrange maxrange maxrange minrange minrange ], ...
%        maxvelocity*[-1 -1 1 1 -1],'--')
  
  axis([0 10 -10 10])
  title(['Please enter up to ' num2str(N) ' targets at min 0.5 m apart using mouse'])
%   str = sprintf('Please enter up to %N targets at min 0.5 m apart using mouse', N);
%   title(str)
  xlabel('x-axis (m)')
  ylabel('y-axis (m)')

% Loop for getting upto N target points 
  xa = 0;
  ya = 0;
  for k=1:N
    [xa,ya] = ginput(1);
%     if Ra < 0 | Ra > mr | V < -mv | V > mv 
%       break
%     end
    plot(xa,ya,'o')
%     ranges = [ranges Ra];
%     vels = [vels Va];
    ranges(k,:) = [[xa ya 0]];
  end

% othersize get input from keyboard
else
  ranges = input('Enter range vector for up to twenty targets: ');
%   vels = input('Enter correponsing velocity vector: ');
%   if length(ranges) ~= length(vels)
%     error('You need the same number of velocities and ranges')
%   end
end

tgtpos = ranges';
tgtvel = zeros(3,N);
tgtrcs = 0.01*ones(N,1)';
% tgtrcs = [0.6 0.3];
% tgtpos


% waveform things

% waveform = phased.RectangularWaveform(...
%     'PulseWidth',1/pulse_bw,...
%     'PRF',prf,...
%     'SampleRate',fs);

% adc
% adc_fs = fs/65536;

waveform = phased.LinearFMWaveform(...
    'PulseWidth',pulse_width,...
    'PRF',prf,...
    'SampleRate',fs,...
    'SweepBandwidth',1e9,...
    'SweepDirection','Up',...
    'Envelope','Rectangular',...
    'OutputFormat','Pulses','NumPulses',1);

figure;
disp(waveform.PulseWidth*waveform.SweepBandwidth)
plot(waveform)
title('Linear FM Wave')
xlabel('Time')
ylabel('Amplitude')


% Phased toolbox code to model hardware
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

snr_min = albersheim(pd, pfa, num_pulse_int)  % some phased magic?

peak_power = radareqpow(lambda,max_range,snr_min,pulse_width,...
    'RCS',tgt_rcs,'Gain',tx_gain)

transmitter = phased.Transmitter(...
    'Gain',tx_gain,...
    'PeakPower',peak_power,...
    'InUseOutputPort',true);

% old antenna
antenna_t = phased.IsotropicAntennaElement(...
    'FrequencyRange',[5e9 40e9]);

% antenna = phased.ULA('Element',antenna_t,'NumElements',10,...
%     'ElementSpacing',lambda/2);
% antenna.Element.FrequencyRange = [5e9 40e9];

antenna = phased.ULA('NumElements',10,...
    'ElementSpacing',lambda/2);
antenna.Element.FrequencyRange = [5e9 40e9];

sensormotion = phased.Platform(...
    'InitialPosition',sensorpos,...
    'Velocity',sensorvel);

radiator = phased.Radiator(...
    'Sensor',antenna,...
    'OperatingFrequency',fc);

collector = phased.Collector(...
    'Sensor',antenna,...
    'OperatingFrequency',fc);


tgtmotion = phased.Platform('InitialPosition',tgtpos,'Velocity',tgtvel);


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

% DOA stuff (not being used)
Nsnapshots = 1000;
rng default
npower = 0.01;

% figure;
% tgtpos = sort(tgtpos);
for m = 1:num_pulse_int

    % Update sensor and target positions
    [sensorpos,sensorvel] = step(sensormotion, 1/prf);
    [tgtpos,tgtvel] = step(tgtmotion, 1/prf);

    % Calculate the target angles as seen by the sensor
    [tgtrng,tgtang] = rangeangle(tgtpos,sensorpos);
    
    % sorting (to match output with input)
    [tgtrng, sortindex] = sort(tgtrng);
    tgtang(1,:) = tgtang(1,sortindex);

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
%     rxasig = collectPlaneWave(antenna, rxsig(:,1), tgtang, fc);
%     rxisig = pulsint(rxsig,'noncoherent');

    % sum individual signals with equal weights(1)
    rxisig = sum(rxsig')';
%     rxasig = collectPlaneWave(antenna, rxisig, tgtang, fc);
    rxpulses(:,m) = step(receiver,rxisig,~(txstatus>0));
end

% rxasig = sensorsig(getElementPosition(antenna)/lambda,...
%     Nsnapshots,tgtang,npower);
% rxasig = collectPlaneWave(antenna, rxasig, tgtang);
    
oldrxpulses = rxpulses;

% automatic threshold detector, more phased magic, need to understand
% what's going on
npower = noisepow(noise_bw,receiver.NoiseFigure,...
    receiver.ReferenceTemperature);
threshold = npower * db2pow(npwgnthresh(pfa,num_pulse_int,'noncoherent'))

% threshold needs to be multipled by some constant! (500)
threshold = threshold*100;

figure;
num_pulse_plot = 2;
helperRadarPulsePlot(rxpulses,threshold,...
    fast_time_grid,slow_time_grid,num_pulse_plot);

% matched filter
matchingcoeff = getMatchedFilter(waveform);
matchedfilter = phased.MatchedFilter(...
    'Coefficients',matchingcoeff,...
    'GainOutputPort',true);
[rxpulses, mfgain] = step(matchedfilter, rxpulses);

matchingdelay = size(matchingcoeff,1)-1;

% more magic
rxpulses = buffer(rxpulses(matchingdelay+1:end),size(rxpulses,1));

threshold = threshold * db2pow(mfgain);

figure;
helperRadarPulsePlot(rxpulses,threshold,...
    fast_time_grid,slow_time_grid,num_pulse_plot);

range_gates = prop_speed*fast_time_grid/2;

% tvg is poor
% tvg = phased.TimeVaryingGain(...
%     'RangeLoss',2*fspl(range_gates,lambda),...
%     'ReferenceLoss',2*fspl(max_range,lambda));
% 
% rxpulses = step(tvg, rxpulses);

figure;
helperRadarPulsePlot(rxpulses,threshold,...
    fast_time_grid,slow_time_grid,num_pulse_plot);

% copying rxpulse for velocity detection
vxpulses = rxpulses;

% integrating all pulses
rxpulses = pulsint(rxpulses,'noncoherent');

helperRadarPulsePlot(rxpulses,threshold,...
    fast_time_grid,slow_time_grid,1);

% range detection
[~,range_detect] = findpeaks(abs(rxpulses),'MinPeakHeight',sqrt(threshold));

true_range = tgtrng
range_estimates = range_gates(range_detect)


% delay = (range_detect-length(pulse))/(fs);
% [rpeaks,rlocs]=pkpicker( abs(rxpulses(:,1)), threshold, 30);

% velocity detection
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

% two methods for doa estimation
% beamscan doa estimator (works poor)
estimator = phased.BeamscanEstimator('SensorArray',antenna,...
    'OperatingFrequency',fc,'ScanAngles',-90:90,...
    'DOAOutputPort',true,'NumSignals',length(tgtrcs));
[y,sigang] = step(estimator,(rxsig));
sigang
figure;
plotSpectrum(estimator)

% super-resolution method (works good)
MUSICestimator = phased.RootMUSICEstimator('SensorArray',antenna,...
    'OperatingFrequency',fc,'NumSignalsSource','Property',...
    'NumSignals',length(tgtrcs),'ForwardBackwardAveraging',true);
doa_est = step(MUSICestimator,(rxsig));

doa_est = sort(doa_est)

% arranging input/output in order of arrival for xy calculations
[tgtang(1,:), sortaindex] = sort(tgtang(1,:));
tgtang(1,:)
range_estimates = range_estimates(sortaindex);


% xy graph
tan_azangles = tand(doa_est);
x_axis = tan_azangles.^2 + 1;
x_axis = range_estimates./sqrt(x_axis)
y_axis = x_axis.*tan_azangles

figure;
plot(x_axis, y_axis,'o')
axis([0 10 -10 10])
grid on
title('Location of Targets With Respect To The Radar at (0, 0)')
xlabel('x axis(m)')
ylabel('y axis(m)')

% tgtdoppler = 0;
% tgtLocation = global2localcoord(tgtpos,'rs',sensorpos);
% tgtazang = tgtLocation(1);
% tgtelang = tgtLocation(2);
% tgtrng = tgtLocation(3);
% tgtcell = val2ind(tgtrng,...
% physconst('LightSpeed')/(2*waveform.SampleRate));
% snapshot = shiftdim(rxpulses(tgtcell,:,:)); % Remove singleton dim
% hadresp = phased.AngleDopplerResponse('SensorArray',antenna,...
% 'OperatingFrequency',fc, ...
% 'PropagationSpeed',physconst('LightSpeed'),...
% 'PRF',prf, 'ElevationAngle',tgtelang);
% plotResponse(hadresp,snapshot);
% text(tgtazang,tgtdoppler,'+Target');


% radar ambiguity
% [afmag_lfm,delay_lfm,doppler_lfm] = ambgfun(pulse,...
% waveform.SampleRate,waveform.PRF);
% surf(delay_lfm*1e6,doppler_lfm/1e3,afmag_lfm,...
% 'LineStyle','none');
% axis tight; grid on; view([140,35]); colorbar;
% xlabel('Delay \tau (\mus)');
% ylabel('Doppler f_d (kHz)');
% title('Linear FM Pulse Waveform Ambiguity Function');

