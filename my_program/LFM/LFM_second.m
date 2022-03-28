clear,clc;
[T,B,Rmin,Rmax,R,RCS]=LFM_radar();%set the default parameter

%% parameter caclulate;
C=3e8;
K=B/T;
Rwid = Rmax-Rmin;
Twid = (2*Rwid)/C;
Fs=5*B; %过采样因子
Ts=1/Fs;
Nwid=ceil(Twid/Ts); %sample numbers

%% generate the echo
t=linspace(2*Rmin/C,2*Rmax/C,Nwid);%receive window
M = length(R);
td = ones(M,1)*t-2*R'/C*ones(1,Nwid);%delay srt
Srt=RCS*(exp(1j*pi*K*td.^2).*(abs(td)<T/2));%strain in the range (-T/2,T/2),radar echo from point targets

%% Digtal processing of pulse compression radar using FFT and IFFT
Nchirp=ceil(T/Ts);
Nfft=2^nextpow2(Nwid+Nwid-1);% why nwid nwid------------------------
Srw=fft(Srt,Nfft);
t0=linspace(-T/2,T/2,Nchirp);
St=exp(1j*pi*K*t0.^2);       
Sw=fft(St,Nfft);%%fft maxlength   
Sot=fftshift(ifft(Srw.*conj(Sw)));

%% output the wave
N0=Nfft/2-Nchirp/2;%%-------------------
Z=abs(Sot(N0:N0+Nwid-1));
Z=Z/max(Z);
Z=20*log10(Z+1e-6);
subplot(211);
plot(t*1e6,real(Srt));axis tight;
xlabel('Time in u sec');ylabel('Amplitude')
title('Radar echo without compression');
subplot(212);
plot(t*C/2,Z);
axis([10000,15000,-60,0]);
xlabel('Range in meters');ylabel('Amplitude in dB')
title('Radar echo after compression');