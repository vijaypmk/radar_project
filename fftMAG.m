function [f X] = fftMAG(x,fs,N,opt)

L = length(x);
%N = 2^(nextpow2(L));
X = fftshift(fft(x,N));
f = (fs)*linspace(-0.5, 0.5, N);

if opt == 1 %non log scale
    figure;
    plot(f,abs(X)/L,'m','linewidth',4);
    title('Magnitude of FFT');
    xlabel('Frequency (Hz)')
    ylabel('Magnitude |X(f)|');
    grid

elseif opt == 2
    figure;
    plot(f,20*log10(abs(X)/(L)),'m','linewidth',4);
    title('Magnitude of FFT');
    xlabel('Frequency (Hz)')
    ylabel('Magnitude |X(f)|  (dB)');
    grid    

else

end

