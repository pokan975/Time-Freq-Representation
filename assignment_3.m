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

[tfr,t1,f] = tfrwv(x);
tfr = 20*log10(abs(tfr));
figure(1)
f = f*fs;
imagesc(t, f, tfr);
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
title('WD of LFM signal x(t) = e^{j2\pi(f_0t+\alpha t^2/2)}')
colormap winter
axis xy

%% 4.(b)
h = gausswin(101);
[tfr, ~, ~] = tfrsp(x, t1, N, h);
ff = fs/N:fs/N:fs;
figure(2)
imagesc(t, ff, tfr)
axis([0 T 0 fs/2]);
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
title('QTFR Spectrogram of signal x(t) = e^{j2\pi(f_0t+\alpha t^2/2)}')
axis xy

%% 4.(c)
g = real(x);
[tfr, ~, ~] = tfrwv(g);
tfr = 20*log10(abs(tfr));
tfr([1:end/2 end/2+1:end], :) = tfr([end/2+1:end 1:end/2], :);
figure(3)
f = f - fs/4;
imagesc(t, f, tfr);
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
title('WD of real signal g(t) = cos2\pi(f_0t+\alpha t^2/2), T_d = 10ms')
colormap winter
axis xy

%% 4.(d)
T = 20e-3;
t = (Ts:Ts:T)';
r = cos(2*pi*(f0*t + (alpha/2)*(t.*t)));
[tfr, ~, ~] = tfrwv(r);
tfr = 20*log10(abs(tfr));
tfr([1:end/2 end/2+1:end], :) = tfr([end/2+1:end 1:end/2], :);
figure(4)
imagesc(t, f, tfr);
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
title('WD of real signal r(t) = cos2\pi(f_0t+\alpha t^2/2), T_d = 20ms')
colormap winter
axis xy

%% 4.(e)
z = hilbert(r);
[tfr, ~, ~] = tfrwv(z);
tfr = 20*log10(abs(tfr));
tfr([1:end/2 end/2+1:end], :) = tfr([end/2+1:end 1:end/2], :);
figure(5)
imagesc(t, f, tfr);
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
title('WD of analytic signal of z(t) = cos2\pi(f_0t+\alpha t^2/2)')
colormap winter
axis xy
