clear
close
clc

% fs = 51.2e3;
% ff = 1:1:256;
t = 1:1:256;
Ntt = length(t);

%% Signal 1 (frequency-overlapping Signal)
t1 = t(75:130);
t2 = t;

x1 = zeros(1,Ntt);
x1(75:130) = 2.7*exp(-1*((9 * (t1 - 102).^2) / 6272) + 1i*((2 * pi * 1000 * (t1 - 70).^2) / (511^2)));
x2 = 1.4*exp(-1*(((t2 - 185).^2) / 338) + 1i*(2*pi*125*t2/511));

x = x1 + x2;

%% Unit-energy circular Gaussian window
% h = gausswin(64+1,1);

%% STFT
% figure(6)
% tfrstft(x.', t, Nt, h);

%% Spectrogram
Nt = 128;
pad = zeros(1, Nt*3/2);
x = [pad, x];
NN = length(x);

% time-domain index, total index = 384+256; first 384 indeices are 0s (negative time index)
V = zeros(Nt, Nt);  % -N/2 ~ N/2-1
B = zeros(Nt, Nt);  % -N/2 ~ N/2-1
E = zeros(Nt, 1);   % -N/2 ~ N/2-1
C = zeros(1, Nt);   % candidate window width
T = zeros(1, NN);   % optimum window width
%% Recursion
for n = 1:NN
    S = E(Nt/2+1, 1);  % E0
    C(1) = E(Nt/2+1, 1);

    for i = 2:Nt/2
        ii = i - 1;  % 1 ~ N/2-1
        S = S + E(Nt/2+1+ii, 1) + E(Nt/2+1-ii, 1);
        C(i) = S/(2*(ii^2) + 4*ii + 1);
    end

    T_opt = max(C);
    T(n) = 2*T_opt;

    V(:, 1:end-1) = V(:, 2:end);
    for j = -Nt/2:Nt/2-1
        aa = 0;
        bb = 0;
        cc = 0;
        dd = 0;
        if n + 1 + j + abs(j) > 0 && n + 1 + j + abs(j) <= NN
            aa = x(n + 1 + j + abs(j));
        end
        if n + 1 - j + abs(j) > 0 && n + 1 - j + abs(j) <= NN
            bb = x(n + 1 - j + abs(j));
        end
        if n + j - abs(j) > 0 && n + j - abs(j) <= NN
            cc = x(n + j - abs(j));
        end
        if n - j - abs(j) > 0 && n - j - abs(j) <= NN
            dd = x(n - j - abs(j));
        end
        jj = j + Nt/2 + 1;  % real j index in MATLAB
        V(jj,end) = aa*conj(bb) - cc*conj(dd);
    end

    b = B(:, 1);
    B(:, 1:end-1) = B(:, 2:end);
    B(: , end) = abs(V(:, end)).^2;

    E = E + B(:, end) - b;
end

T = round(T);
T = T(end-Ntt+1:end);

figure(1)
plot(T)
