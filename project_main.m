clear
close
clc

% alpha = 1;
% fs = 1024;
% t = (0:1/fs:1-1/fs)';
% 
% gau = exp(-1*alpha*(t.^2)/2);
% gauh = hilbert(gau);
% 
% tfrstft(gauh);

fs = 1024; Td = 2; 
chirp_rate = 37.5; 
fi1 = 110; fi2 = 30; 

x1 = fmlin(Td*fs, fi1/fs, (fi1 + 2*chirp_rate*Td)/fs); 
x2 = fmlin(Td*fs, fi2/fs, (fi2 + 2*chirp_rate*Td)/fs); 
x = x1 + x2;

N = length(x); 
h = gausswin(101); 
tfr3 = tfrstft(x, 1:N, N, h);
tfr3 = flipud(tfr3);

figure(1) 
tt = linspace(0, Td, N);  
ff = linspace(0, fs, N); 
imagesc(tt, ff, abs(tfr3(1025:end,:)));
% axis([0 Td 0 fs/2]); 
% axis xy
xlabel('Time (sec)'), ylabel('Frequency (Hz)') 
title('STFT of the summed signal, low pass window size = 101')

% =========================================================================
img = abs(tfr3);
[row, col] = size(img);
img_max = max(max(img));
img_min = min(min(img));

k = mat2gray(img);
grey = round(k*256);
figure(2)
imshow(grey)