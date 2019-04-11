clear
close
clc

fs = 51.2e3;
ff = 1:1:256;
t = 1:1:256;
Nt = length(t);

%% Signal 1 (frequency-overlapping Signal)
t1 = t(75:130);
t2 = t;

x1 = zeros(1,Nt);
x1(75:130) = 2.7*exp(-1*((9 * (t1 - 102).^2) / 6272) + 1i*((2 * pi * 1000 * (t1 - 70).^2) / (511^2)));
x2 = 1.4*exp(-1*(((t2 - 185).^2) / 338) + 1i*(2*pi*125*t2/511));

x = (x1 + x2).';

%% Unit-energy circular Gaussian window
h = gausswin(64+1,1);

%% STFT
figure(6)
tfrstft(x, t, Nt, h);

%% Spectrogram
T = zeros(Nt, 1);
pad = zeros(Nt*3/2, 1);
x = [pad; x];

n = 0;
V = zeros(Nt, Nt);
B = zeros(Nt, Nt);
E = zeros(Nt, 1);

%% Recursion
C = zeros(1, Nt);  % candidate window width
S = E(Nt/2+1, 1);  % E0
C(0) = E(Nt/2+1, 1);

for i = 2:Nt
    ii = i - 1;
    S = E;
    C(i) = S/(2*(ii^2) + 4*ii + 1);
end

T_opt = max(C);

if n > 384
    
end