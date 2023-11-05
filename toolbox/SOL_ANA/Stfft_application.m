
% ultilize stfft method to indicate the time dependent frequency
clc;
global fs
% close all;
% fs = 10000;              % 采样频率
dt = 1/fs;              % 采样时间间隔
t=dt.*length(yout(:,3));
f=yout(:,3);

% serial=find(t<20);
% t=t(serial);
% f=f(serial);
% Nt=2^ceil(log2(length(f)));
% Nt=2^floor(log2(length(f)));

Nt=length(t);
f1=hilbert(f);
% 
SNR=20;
window=4096;
noverlap=window-1024;

[S,F,T,P]=spectrogram(f1,window,noverlap,10*fs,fs,'MinThreshold',-SNR,'yaxis');
figure
imagesc(T,F,abs(S))
ylim([0 200])
axis xy;
box on;
colorbar;
xlabel('Time/s');ylabel('Frequency');title('Stfft')
set(gca,'fontname','times new roman','fontsize',12)

% [q,nd] = max(10*log10(P));
% hold on
% plot3(t,f(nd),q,'r','linewidth',4)
% hold off