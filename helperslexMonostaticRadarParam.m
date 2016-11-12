function helperslexMonostaticRadarParam
% This function is only in support of slexMonostaticRadarExample. 
% It may be removed in a future release.

%   Copyright 2014 The MathWorks, Inc.
    clear

    [propSpeed, fc, pulseBw, prf, fs, txGain, peakPower, ...
     matchingCoeff, metersPerSample, rangeOffset, rangeLoss, ...
     referenceLoss, target1Rcs, target1Pos, target1Vel, lambda] = calcParams();
    
    % Environment
    paramRadar.propSpeed = propSpeed;
    paramRadar.fc = fc;
    paramRadar.lambda = lambda;
    % Waveform parameters
    paramRadar.pulseBw = pulseBw;
    paramRadar.prf = prf;
    paramRadar.fs = fs;
    % Transmitter parameters
    paramRadar.txGain = txGain;
    paramRadar.peakPower =  peakPower;
    % Matched filter parameters
    paramRadar.matchingCoeff = matchingCoeff;
    % Time varying gain parameters 
    paramRadar.metersPerSample = metersPerSample;
    paramRadar.rangeOffset = rangeOffset;
    paramRadar.rangeLoss = rangeLoss;
    paramRadar.referenceLoss = referenceLoss;
    % Radar parameters
    paramRadar.target1Rcs = target1Rcs;
    paramRadar.target1Pos = target1Pos;
    paramRadar.target1Vel = target1Vel;

    assignin('base','paramRadar',paramRadar);

end

function [propSpeed, fc, pulseBw, prf, fs, txGain, peakPower, ...
          matchingCoeff, metersPerSample, rangeOffset, rangeLoss, ...
          referenceLoss, target1Rcs, target1Pos, target1Vel, lambda]  = calcParams()  
    % Environment
    propSpeed = physconst('LightSpeed');   % Propagation speed
    fc = 24e9;           % Operating frequency - old 10e9
    lambda = propSpeed/fc;


    % Constraints
    maxRange = 50;    % Maximum unambiguous range old - 5000
    rangeRes = 0.5;      % Required range resolution old - 50
    pd = 0.9;            % Probability of detection
    pfa = 1e-6;          % Probability of false alarm
    tgtRcs = 0.25;         % Required target radar cross section old - 1
    numPulseInt = 11;  % Integrate 10 pulses at a time


    % Waveform parameters
    pulseBw = propSpeed/(2*rangeRes);    % Pulse bandwidth
    pulseWidth = 1/pulseBw;               % Pulse width
    prf = propSpeed/(2*maxRange);         % Pulse repetition frequency
    fs = 2*pulseBw;    

    % Transmitter parameters
    snrMin = albersheim(pd, pfa, numPulseInt);
    txGain = 20;  % old 20
    peakPower =  ...
        radareqpow(lambda,maxRange,snrMin,pulseWidth,...
                   'RCS',tgtRcs,'Gain',txGain);

    % Matched filter parameters
    hwav = phased.RectangularWaveform(...
        'PulseWidth',1/pulseBw,...
        'PRF',prf,...
        'SampleRate',fs);
    matchingCoeff = getMatchedFilter(hwav);

    % Delay introduced due to filter
    matchingDelay = size(matchingCoeff,1)-1;

    % Time varying gain parameters 
    fastTimeGrid = unigrid(0,1/fs,1/prf,'[)');
    rangeGates = propSpeed*fastTimeGrid/2; 
    metersPerSample = rangeGates(2);
    rangeOffset = -rangeGates(2)*matchingDelay;
    rangeLoss = 2*fspl(rangeGates,lambda);
    referenceLoss = 2*fspl(maxRange,lambda);

    %Radar parameters
    target1Rcs = [0.3 0.4 0.5 0.66 0.5];
%      target1Rcs = [0.5];
%     target1Pos = [[10;0;0],[0;0;0],[0;0;0]];  % old 1988.66
%     target1Vel = [[0;0;0], [0;0;0],[0;0;0]];
%     target1Pos = [2.5;0;0];
%     target1Vel = [0;0;0];
    target1Pos = [1.98 3.53 4.5 10.45 6.1;...
                                    0 0 0 0 0; ...
                                    0 0 0 0 0];
    target1Vel = zeros(3,5);


%     target1Rcs = [0.6 2.2 1.05 .5];
%     target1Pos = [1988.66 3532.630 3845.04 1045.04;...
%                                     0 0 0 0 ; ...
%                                     0 0 0 0 ];
%     target1Vel = zeros(3,4);

end
