clear
close
clc

%% (a)
fs = 1024; Td = 2;
chirp_rate = 37.5;
fi1 = 110; fi2 = 30;

x1 = fmlin(Td*fs, fi1/fs, (fi1 + 2*chirp_rate*Td)/fs);
x2 = fmlin(Td*fs, fi2/fs, (fi2 + 2*chirp_rate*Td)/fs);

x = x1 + x2;
[instf, tt] = instfreq(x);

figure(1)
t = Td/fs:1/fs:Td-1/fs;
plot(t, instf*fs, 'b')
xlabel('Time (sec)'),ylabel('Frequency (Hz)')
title('Instant. Frequency of the sum of the 2 LFM chirp signals')
axis([0 2 0 500])

%% (b)
N = length(x);
h = gausswin(21);
tfr1 = tfrstft(x, 1:N, N, h);
h = gausswin(61);
tfr2 = tfrstft(x, 1:N, N, h);
h = gausswin(101);
tfr3 = tfrstft(x, 1:N, N, h);

figure(2)
tt = linspace(0, Td, N); 
ff = linspace(0, fs, N);
subplot(311)
imagesc(tt, ff, abs(tfr1));
axis([0 Td 0 fs/2]);
axis xy
xlabel('Time (sec)'), ylabel('Frequency (Hz)')
title('STFT of the summed signal, low pass window size = 21')
subplot(312)
imagesc(tt, ff, abs(tfr2));
axis([0 Td 0 fs/2]);
axis xy
xlabel('Time (sec)'), ylabel('Frequency (Hz)')
title('STFT of the summed signal, low pass window size = 61')
subplot(313)
imagesc(tt, ff, abs(tfr3));
axis([0 Td 0 fs/2]);
axis xy
xlabel('Time (sec)'), ylabel('Frequency (Hz)')
title('STFT of the summed signal, low pass window size = 101')