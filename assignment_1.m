clear
close
clc

%% (a)
fs = 800;
t = (1/fs:1/fs:2)';
tt = 2/fs:1/fs:2-1/fs;

x1 = exp(1i*2*pi*50*t);
[instf1,t1] = instfreq(x1);
figure(1)
plot(tt,instf1*fs)
xlabel('Time (sec)'),ylabel('Frequency (Hz)')
title('Time-Frequency Representation of x_1(t) = e^{j2\pi50t}')
axis([0 2 0 200])

% x2 = exp(1i*2*pi*75*t);
x2 = fmconst(2*fs, 75/fs);
[instf2,t2] = instfreq(x2);
figure(2)
plot(tt,instf2*fs)
xlabel('Time (sec)'),ylabel('Frequency (Hz)')
title('Time-Frequency Representation of x_2(t) = e^{j2\pi75t}')
axis([0 2 0 200])

x3 = exp(1i*2*pi*(20*t + 40*(t.^2)));
% x3 = fmlin(2*fs, 20, 180);
[instf3,t3] = instfreq(x3);
figure(3)
plot(tt,instf3*fs)
xlabel('Time (sec)'),ylabel('Frequency (Hz)')
title('Time-Frequency Representation of x_3(t) = e^{j2\pi(20t+40t^2)}')
axis([0 2 0 200])

x4 = x1 + x2 + x3;
[instf4,t4] = instfreq(x4);
figure(4)
plot(tt,instf4*fs)
xlabel('Time (sec)'),ylabel('Frequency (Hz)')
title('Time-Frequency Representation of x(t) = \Sigma_{i=1}^{3}x_i(t)')
axis([0 2 0 400])

%% (b)
figure(5)
subplot(211)
plot(tt,instf4*fs)
xlabel('Time (sec)'),ylabel('Frequency (Hz)')
title('x(t) = \Sigma_{i=1}^{3}x_i(t)')
axis([0 2 0 400])
subplot(212)
plot(tt,instf1*fs,'b')
hold on
plot(tt,instf2*fs,'r')
hold on
plot(tt,instf3*fs,'g')
hold off
xlabel('Time (sec)'),ylabel('Frequency (Hz)')
legend('x_1(t)','x_2(t)','x_3(t)')
title('superimposition of x_1(t), x_2(t), x_3(t)')
axis([0 2 0 200])

% Discussion:
% the synthesized signal x(t) causes a distortion in instantaneous frequency estimation,
% the trend of synthesized frequency can be roughly observed, but
% it cannot represent the 3 different components of x1, x2, and x3 (in subplot(212))