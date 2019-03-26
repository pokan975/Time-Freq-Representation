clear
close
clc

% T = 10e-3;
fs = 51.2e3;
% Ts = 1/fs;
% t = Ts:Ts:T;
ff = 1:1:512;
t = 1:1:512;
Nt = length(t);

%% Signal 1 (Frequency-Overlapping Signal)
t1 = t(75:130);
t2 = t;
t3 = t(245:445);

x1 = zeros(1,Nt);
x1(75:130) = 2.7*exp(-1*((9 * (t1 - 102).^2) / 6272) + 1i*((2 * pi * 1000 * (t1 - 70).^2) / (511^2)));
x2 = 1.4*exp(-1*(((t2 - 185).^2) / 338) + 1i*(2*pi*125*t2/511));
x3 = zeros(1,Nt);
x3(245:445) = exp(-1*(((t3 - 344).^2) / 80000) + 1i*(2*pi*125*t3/511));

x = (x1 + x2 + x3).';

%% Signal 2 (Time-Overlapping Signal)
t1 = t(100:400);
t2 = t(125:375);
t3 = t(200:300);

s1 = zeros(1,Nt);
s1(100:400) = 0.5*exp(-1*(((t1 - 250).^2) / 125000) + 1i*(2*pi*30*t1/511));
s2 = zeros(1,Nt);
s2(125:375) = 0.8*exp(-1*((9*(t2 - 250).^2) / 781250) + 1i*((pi*4000*(t2 - 250).^3) / (3*511^3)) + 1i*((2*pi*70*(t2 - 250)) / 511));
s3 = zeros(1,Nt);
s3(200:300) = exp(-1*(((t3 - 250).^2) / 450) + 1i*(2*pi*130*t3/511));

s = (s1 + s2 + s3).';

%% Wigner Distribution
figure(1)
plot(t, real(x))
xlabel('time samples'), ylabel('real(x)')
axis([0 512 -3 3])

figure(2)
plot(t, real(s))
xlabel('time samples'), ylabel('real(x)')
axis([0 512 -3 3])

[wd_x,~,f] = tfrwv(x);
figure(3)
f = f*fs;
imagesc(t, f, wd_x);
xlabel('Time Samples')
ylabel('Frequency (Hz)')
title('WD of signal x(t)')
axis xy

[wd_s,~,f] = tfrwv(s);
figure(4)
imagesc(t, f, wd_s);
xlabel('Time Samples')
ylabel('Frequency (Hz)')
title('WD of signal s(t)')
axis xy

%% Unit Energy Gaussian Window
h = gausswin(64+1,1);
% figure(5)
% tfrwv(h);

%% STFT
% stft_x = tfrstft(x, t, Nt, h);
% figure(6)
% imagesc(t, ff, abs(stft_x));
% axis xy
% xlabel('Time Samples'), ylabel('Frequency (Hz)')
% axis([1 Nt 0 Nt/2]);

figure(6)
tfrstft(x, t, Nt, h);