clear; close; clc

fs = 1024;
Ts = 1/fs;
Td = 1;
t = 1:Td/Ts;
Nt = length(t);

tt = linspace(0,Td,Nt);
ff = linspace(0,fs,Nt);

%% Signal 1 (frequency-overlapping Signal)
sin_start = 100;
chirp_start = 600;
chirp_rate = 700;
chirp_Td = 200;
fi1 = 200; 
fi2 = 50;

% linear chirp signals
x1 = fmlin(chirp_Td, fi1/fs, (fi1 + 2*chirp_rate*chirp_Td*Ts)/fs);
x2 = fmlin(chirp_Td, fi2/fs, (fi2 + 2*chirp_rate*chirp_Td*Ts)/fs);
% sinusoid signals
x3 = fmconst(chirp_Td,100/fs);
x4 = fmconst(chirp_Td,200/fs);

x = zeros(Nt,1);
% add all signals
x(chirp_start+1:chirp_start + chirp_Td) = x1 + x2;
x(sin_start+1:sin_start + chirp_Td) = x3 + x4;


%% Candidate window size
window = [16,32,48,64,80,100];
w_window = max(window) + 9;

len = length(window);
cand = cell(1,len);
for i = 1:len
    h = hann(window(i)+1);
    stft = tfrstft(x, t, Nt, h);
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

figure(1)
plot(tt, opt_len)
xlabel('Time (sec)'),ylabel('window length')
axis([0 Td 0 120]);

%% for test only
% sinusoid signal
h = hann(opt_len(sin_start + chirp_Td/2) + 1);
stft = tfrstft(x, t, Nt, h);
tfr1 = stft(:,1:end/2);
% linear chirp signal
h = hann(opt_len(chirp_start + chirp_Td/2) + 1);
stft = tfrstft(x, t, Nt, h);
tfr2 = stft(:,end/2+1:end);

tfr = [tfr1,tfr2];

figure(2)
imagesc(tt, ff, abs(tfr))
axis xy
xlabel('Time (sec)'), ylabel('Frequency (Hz)')
axis([0 Td 0 fs/2]);

figure(3)
h = hann(64+1);
stft = tfrstft(x, t, Nt, h);
imagesc(tt, ff, abs(stft));
axis xy
xlabel('Time (sec)'), ylabel('Frequency (Hz)')
axis([0 Td 0 fs/2]);