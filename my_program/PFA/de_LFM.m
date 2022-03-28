close all;clear clc;
%% setting parameter
c=3e8;
delta_r=50;    
Rref = 200;
H=100;
rmin = Rref-delta_r/2;
rmax = Rref+delta_r/2;  
B = 200e6;
Tp = 1.5e-6;
Tref=2e-6;
K = B/Tp;
range = (rmax.^2-rmin.^2).^0.5;
v=100;
fc=5e9;
lambda=c/fc;
Kd=2*v^2/(rmin*lambda);
Bmax=(range/v)*Kd;
PRF=4*Bmax;     
N = 2*delta_r*B/c;
Ts=Tp/N;
T=1e-3;
tm=-range/v:1/PRF:range/v-1/PRF;
tq=-Tp/2+2*rmin/c:Ts:Tp/2-Ts+2*rmin/c;
A0=1;

%% send ref and return signal SAR
R=(rmin.^2+(v*tm).^2).^0.5;
s_ret = zeros(length(tm),length(tq));
s_ref = zeros(length(tm),length(tq));
s = zeros(length(tm),length(tq));
for i=1:length(tm)
    %tref=-Tref/2+2*Rref/c:Ts:-Tref/2+2*Rref/c-Ts;
    %s_ref(i,:)=A0*rectpuls(tq-2*Rref/c,Tref).*(1j*2*pi.*(fc.*(tq-2*Rref/c)+0.5*K*(tq-2*Rref/c).^2));
    s_ret(i,:)=rectpuls(tq-2*R(i)/c,Tp).*exp(1j*2*pi.*(fc.*(tq-2*R(i)/c)+0.5*K*(tq-2*R(i)/c).^2));
    %s(i,:)=s_ref(i,:).*conj(s_ret(i,:));
end
Y=zeros(length(tm),length(tq));
f=zeros(1,length(tm));
for i=1:length(tm)
    tt=0:Ts:Tp-Ts;
    hk=exp(-1j*pi*tt.^2.*K);
    %N0 = 2^nextpow2(length(tq));
    Y(i,:)=ifft(fft(s_ret(i,:)).*conj(fft(hk,length(hk))));
end
x_r=tq*c/2;
y_a=v*tm;
[R,q]=meshgrid(x_r,y_a);
mesh(R,q,abs(Y));view(0,90);xlim([0 300]);
xlabel('距离向');
% ylabel('合成孔径 (慢时间), meters');
title('距离压缩后');
%% tq fft
y = fft(s);



