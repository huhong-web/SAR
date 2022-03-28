close all;
clear;
clc;
%% setting parameter
c=3e8;   
H=100; 
B = 100e6;%发送脉宽
Tp = 10e-6;%脉冲持续时间
K = B/Tp;%调频率
La=150;
range=1000;
target=[range-50 100 0;range 200 0;range+50 200 0];
v=100;
fc=5e9;
rmin=sqrt(sum(target(1,:).^2));
r_d=sqrt(sum(target(3,:).^2));
rmax=sqrt(r_d.^2+(La/2).^2);
y_min=min(target(:,2))-La/2;
y_max=max(target(:,2))+La/2;
limit=y_max-y_min;
lambda=0.24;%sending wavelength
Kd=2*v^2/(rmax*lambda);
Bmax=(La/v)*Kd;
PRF=1.2*Bmax;     
Ts=1/(1.25*B);
tm=y_min/v:1/PRF:y_max/v-1/PRF;
tq=-Tp/2+2*rmax/c:Ts:Tp/2-Ts+2*rmax/c;
A0=1;
%% send ref and return signal SAR
y_a=v*tm;
s_ret = zeros(length(tm),length(tq),length(target));
s_echo = zeros(length(tm),length(tq));
for j=1:length(target)
    for i=1:length(tm)
        y_coordi = [0,y_a(i),H];
        R=sqrt(sum((target(j,:)-y_coordi).^2));
        s_ret(i,:,j)=(abs(target(j,2)-y_a(i))<=La/2)*rectpuls(tq-2*R/c,Tp).*...
          exp(1j*2*pi.*(fc.*(tq-2*R/c)+0.5*K*(tq-2*R/c).^2));
    end
    s_echo=s_echo+s_ret(:,:,j);
end
for i=1:size(s_ret,1)
    s_echo(i,:)=s_echo(i,:).*exp(-1j*2*pi*fc.*tq);
end
tt=0:Ts:Tp-Ts;%%   how to choose the range
Y=zeros(length(tm),length(tq));
for i=1:length(tm)
    hk=exp(1j*pi*tt.^2.*K);
    Y(i,:)=ifft(fft(s_echo(i,:)).*conj(fft(hk,length(tq))));
end
x_r=sqrt((tq*c/2).^2-H.^2);
[R_delta,q]=meshgrid(x_r,y_a);
mesh(R_delta,q,abs(Y));
view(0,90);
xlabel('距离向');
% ylabel('合成孔径 (慢时间), meters');
title('距离压缩后');



