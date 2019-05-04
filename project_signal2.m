clear; close; clc

fs = 1024;
Td = 1;
t = 1:Td*fs;
Nt = length(t);
ff = linspace(0,fs,Nt);
tt = linspace(0,Td,Nt);

window = [16,32,48,64,80,100];
w_window = max(window) + 9;

%% Signal 1 (frequency-overlapping Signal)
t1 = t(75:130);
t2 = t;
t3 = t(645:845);

x1 = zeros(1,Nt);
x1(75:130) = 2.7*exp(-1*((9 * (t1 - 102).^2) / 6272) + 1i*((2 * pi * 1000 * (t1 - 70).^2) / (511^2)));
x2 = 1.4*exp(-1*(((t2 - 385).^2) / 338) + 1i*(2*pi*125*t2/511));
x3 = zeros(1,Nt);
x3(645:845) = exp(-1*(((t3 - 744).^2) / 80000) + 1i*(2*pi*125*t3/511));

x = x1 + x2 + x3;

h = hann(16+1);
stft = tfrstft(x.', t, Nt, h);
figure(1)
imagesc(tt, ff, abs(stft));
axis xy
xlabel('Time (sec)'), ylabel('Frequency (Hz)')
axis([0 Td 0 fs/2]);


%% Candidate window size
len = length(window);
cand = cell(1,len);
for i = 1:len
    h = hann(window(i)+1);
    stft = tfrstft(x.', t, Nt, h);
    cand{1,i} = stft;
end

%% Calc L-2, L-4 norm
l2norm = cell(1,len);
l4norm = cell(1,len);
for i = 1:len
    tfr = cand{1,i};
    norm2 = zeros(1,Nt);
    norm4 = zeros(1,Nt);
    for j = 1:Nt
        norm2(j) = sum(abs(tfr(:,j)).^2);
        norm4(j) = sum(abs(tfr(:,j)).^4);
    end
    l2norm{1,i} = norm2;
    l4norm{1,i} = norm4;
end

%% Calc energy concentration
L2 = zeros(len,Nt);
L4 = zeros(len,Nt);
pad = (w_window-1)/2;
for i = 1:len
    norm2_pad = [zeros(1,pad),l2norm{i},zeros(1,pad)];
    norm4_pad = [zeros(1,pad),l4norm{i},zeros(1,pad)];
    for j = 1:Nt
        L2(i,j) = sum(norm2_pad(1,pad+j:pad*2+j));
        L4(i,j) = sum(norm4_pad(1,pad+j:pad*2+j));
    end
end

concen = L4./(L2.^2);
opt_len = zeros(1,Nt);
astft = zeros(size(stft));
for i = 1:Nt
    temp = concen(:,i);
    [m,index] = max(temp);
    opt_len(i) = window(index);
    astft(:,i) = cand{index}(:,i);
end

figure(2)
plot(tt, opt_len)
xlabel('Time (sec)'),ylabel('window length')
axis([0 Td 0 120]);

figure(3)
imagesc(tt, ff, abs(astft))
axis xy
xlabel('Time (sec)'), ylabel('Frequency (Hz)')
axis([0 Td 0 fs/2]);