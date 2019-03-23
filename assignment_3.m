clear; close; clc

alpha = 2e6;
f0 = 20e3;
T = 10e-3;
fs = 102.4e3;
Ts = 1/fs;

%% 4.(a)
t = (Ts:Ts:T)';
x = exp(1i * 2 * pi * (f0*t + (alpha/2)*(t.*t)));
N = length(x);

figure(1)
plot(t, real(x))

[tfr,t1,f] = tfrwv(x);
figure(2)
f = f*fs;
imagesc(t, f, tfr);
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
title('WD of LFM signal x(t) = e^{j2\pi(f_0t+\alpha t^2/2)}')
axis xy

%% 4.(b)
h = gausswin(101);
[tfr, ~, ~] = tfrsp(x, t1, N, h);
ff = fs/N:fs/N:fs;
figure(3)
imagesc(t, ff, tfr)
axis([0 T 0 fs/2]);
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
title('QTFR Spectrogram of signal x(t) = e^{j2\pi(f_0t+\alpha t^2/2)}')
axis xy

%% 4.(c)
g = real(x);
[tfr, ~, ~] = tfrwv(g);
figure(4)
imagesc(t, f, tfr);
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
title('WD of real signal g(t) = cos2\pi(f_0t+\alpha t^2/2), T_d = 10ms')
axis xy

%% 4.(d)
T = 20e-3;
t = (Ts:Ts:T)';
r = cos(2*pi*(f0*t + (alpha/2)*(t.*t)));
[tfr, ~, ~] = tfrwv(r);
figure(5)
imagesc(t, f, tfr);
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
title('WD of real signal r(t) = cos2\pi(f_0t+\alpha t^2/2), T_d = 20ms')
axis xy

%% 4.(e)
z = hilbert(r);
[tfr, ~, ~] = tfrwv(z);
figure(6)
imagesc(t, f, tfr);
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
title('WD of analytic signal of z(t) = cos2\pi(f_0t+\alpha t^2/2)')
axis xy