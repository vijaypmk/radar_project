% University of California, Santa Cruz
%
% LFM Pulse Doppler Radar for Range and Velocity
%---------------------------------------------------------------
% Author - Vijay Muthukumaran
% Course - EE-288 Radar, SAR and ISAR
%---------------------------------------------------------------

clear,  format compact
T = 7;                    % pulse length usec
W = 7;                    % bandwitdh MHz
fs = 8;                   % sampling frequency MHz
c = 0.3;                  % speed of light in Km/usec
RF = 7e9;                 % not being used

% creating chirp
[s, nn] = chirp(T, W, fs/W);

Np = 11;                  % 11 pulses
jkl = 0:(Np-1);

pri = 60;                 % PRI usec
T_0 = pri*jkl;            % time of start of each pulse usec
g = ones(1,Np);           % gains  (makes frequency linear)
T_out = [25 50];          % start and receive window in usec
T_ref = 0;                
fc = 7000;                % center frequency in MHz

% my own range and velocity
ranges = [4.5 5.7];       % Km
Ntargets = length(ranges);
amps = 10*[1 1];          % amplitude
vels = [100 50];          % velocity in m/sec

% use for user defined range and velcoity
% y = radar(s,fs,T_0,g,T_out,T_ref,fc,ranges,amps,vels);

% use for given .mat file
temp = load('r100.mat');
y = temp.y;

%% -----------------------------------------------------------------------

% Ask user input
noise = 0;
noise = input('Do you want to add noise to system for estimation analysis? (1/0)');

if noise
    % standard deviation
    sigma = input('Standard Deviation of noise? (1 - 20)');
else
    % default to 0
    sigma = 0;
end

% Creating noise for doing estimation analysis
[M,N]=size(y);
% % Add white noise (Gaussian) with std dev of sigma
yn = y + sigma*randn(M,N);
yn = yn + sigma*j*randn(M,N);

%% -------------------------------------------------------------------------

figure;
plot(nn, real(s))
title('Graph of the LFM Chirp')
xlabel('Time in usec')
ylabel('Amplitude ')
grid on;

ny = 1:length(y);
figure;
plot(ny, abs(y))
title('Graph of Returned Radar Signal and Matched Filter Output')
xlabel('Time in usec')
ylabel('Amplitude')
grid on;
hold on;

%% -----------------------------------------------------------------------

% filter section
h = conj(fliplr(s));         % use conjugate of flipped LFM signal
e = [];
% 3 pulse delay canceller MTI filter
hf = [1 -2 1];               
ymti = filter(hf,1,yn,[],1);

for i=1:Np
    % matched filter
    e(:, i) = filter(h, 1, ymti(:, i)); 
end

ex = [0:1:length(e) - 1];
plot(ex, abs(e))

%% -----------------------------------------------------------------------

% Threshold Detection (currently not being used)
% Tried different methods, explained in report
% YN = []
% for q = 1:Np
%     YN(q)= sum(abs(y(:,q)));
% end
% YNN = sum(YN);
% YNNN = (1/Np)*abs(YNN);

% ee1 = sum(e(1:length(e)/2),1);
% ee1 = sum(ee1);
% ee2 = sum(e(length(e)/2 + 1 : end),1);
% ee2 = sum(ee2);
% ee1 = (1/11)*abs(ee1);
% ee2 = (1/11)*abs(ee2);
% eet = min(ee1, ee2);

% ee = sum(e, 1);
% ee = sum(ee);
% eet = (4/11)*abs(ee);

%% -----------------------------------------------------------------------

% Calculating range for every pulse
range_thresh = 300;              %Eyeballed threshold;

R = [];
flag = 1;
for i = 1:(Np)
    [rpeaks,rlocs]=pkpicker( abs(e(:,i)), range_thresh, 30);        
    
    % Calculating time delay for range
    delay = (rlocs-length(s))/(fs) + T_out(1); 

    range = (c/2)*delay;
    
    if flag
        r0_locs = rlocs;
        flag = 0;
        R(:,1) = [range];
    end
end

% Print range obtained from first pulse
disp('Range')
disp(R(:,1)')
%% -----------------------------------------------------------------------

% Calculating for velocity
vd = e(r0_locs,:);         % use only rows where the peaks are located

fd = [];
V = [];
NFFT = Np*64;

figure;
for k=1:length(r0_locs)
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
    fd(k) = W(loc)/(pri*2*pi*fc);
    V(k)=(c*1000/2)*(fd(k)*1e6);
end

% Print velocity
disp('Velocity')
disp(V)

%% -----------------------------------------------------------------------

% for contour plot
cont = [];
contH = [];
contW = [];
nf_d = [];

% calculate dtft for range data across all pulses
for j = 1:length(e)
    cont = transpose(e);
    curr = cont(:,j);
    [contH(j,:), contW(j,:)] = dtft(curr,NFFT);
    nf_d(j,:) = contW(j,:)/(2*pi*fc*pri);
end

% Calculating axis for contour plot
Rp = linspace(2.4,7,201);
[X Y] = meshgrid(Rp, nf_d(1,:));

figure;
contour(X, Y, abs(20*log10(contH')))
title('Contour plot of Range-Doppler Frequency')
xlabel('Range in Km')
ylabel('Doppler Frequency in MHz')

figure;
mesh(X, Y, abs((contH')))
title('Mesh plot of Range-Doppler Frequency')
xlabel('Range in Km')
ylabel('Doppler Frequency')
zlabel('Magnitude')
% ------------------------------------------------------------------------