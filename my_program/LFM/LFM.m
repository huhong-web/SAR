clear,clc;
N = 100000;% 为采样点数
Tp = 1e-6;%脉冲信号持续时间
i = 0:Tp/N:Tp-Tp/N;
A = 1;%amplitude
B = 30e6;%bandwidth khz
k = B/Tp;
f0 = 10000;
s = A*exp(1i*2*pi.*(f0.*i+0.5*k.*i.^2));
%plot(i,real(s));
y = fft(real(s));%spectrum tool function 1. fft(data dft) spectrum real amplitude
%need to divide N
f = (0:N-1)*(1/Tp);
plot(f,abs(y)/N);        % 2 freqz h(n) = b1*e^jw+b2*e^j2w equal to ft. 3.10030 11450

