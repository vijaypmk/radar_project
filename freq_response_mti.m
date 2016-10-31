b = [1 -2 1];
a = 1;
nnft = Np*64;
% nnft = 1024;
% [f1, yDATA] = dtft(yData, nnft);
[H f2] = freqz(b, a, nnft);
figure;
% subplot(211)
%     plot(w*Fs/2/pi,20*log10(abs(H)))
%     hold on
%     plot(f,20*log10(abs(yDATA(:,1))),'m-')
%     grid, xlabel('Frequency (Hz)'),ylabel('Magnitude Squared (dB)')
%     title('Magnitude Frequency Response')
% subplot(212)
%     plot(w*Fs/2/pi,angle(H))
%     grid, xlabel('Frequency (Hz)'),ylabel('Phase (rad)')
%     title('Phase Frequency Response')

% plot(ex(1:1:10), 20*log10(abs(e(1:1:10,1))))
% hold on
plot(f2, 20*log10(abs(H)))
title('Frequecny Response of MTI Filter')
xlabel('Frequency in MHz')
ylabel('Log Magnitude')
grid on;