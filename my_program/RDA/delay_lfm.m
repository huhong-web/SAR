close all;
clear;
clc;
%% setting parameter
c=3e8;   
H=100; 
B = 75e6;%发送脉宽
Tp = 10e-6;%脉冲持续时间
K = B/Tp;%调频率
La=200;
range=10000;
target=[range 200 0;range-50 100 0;range+50 200 0];
v=100;
fc=5e9;
rc=sqrt(range.^2+H.^2);
r_min=sqrt(sum(target(2,:).^2));
r_d=sqrt(sum(target(3,:).^2));
rmax=sqrt(r_d.^2+(La/2).^2);
y_min=min(target(:,2))-La/2;
y_max=max(target(:,2))+La/2;
limit=y_max-y_min;
lambda=c/fc;%sending wavelength
Kd=2*v^2/(rc*lambda);
Bmax=(La/v)*Kd;
PRF=4*Bmax;     
Ts=1/(2*B); %%%%%采样率太低为什么造成双图像
fs=1/Ts;
tm=0:1/PRF:y_max/v-1/PRF;
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
        s_ret(i,:,j)=A0*(abs(target(j,2)-y_a(i))<=La/2)*rectpuls(tq-2*R/c,Tp).*...
          exp(1j*2*pi.*(fc.*(tq-2*R/c)+0.5*K*(tq-2*R/c).^2));
    end
    s_echo=s_echo+s_ret(:,:,j);
end
for i=1:size(s_ret,1)
    s_echo(i,:)=s_echo(i,:).*exp(-1j*2*pi*fc.*tq);
end
tt=0:Ts:Tp-Ts;%%   how to choose the range
Y=zeros(size(s_echo));
for i=1:length(tm)
    hk=exp(1j*pi*tt.^2*K);
    Y(i,:)=ifft(fft(s_echo(i,:)).*conj(fft(hk,length(tq))));
end
figure;
x_r=sqrt((tq*c/2).^2-H.^2);
[R_delta,q]=meshgrid(x_r,y_a);
mesh(R_delta,q,abs(Y));
view(0,90);
xlabel('距离向');
% ylabel('合成孔径 (慢时间), meters');
title('距离压缩后');

%%方位向fft平移
S_shift=Y;
fftshift(fft(S_shift),1);
hk=exp(-1j*pi*Kd*(0:1/PRF:La/v-1/PRF).^2);
Hk=fft(hk,length(tm))';
for i=1:length(tq)
    Y(:,i)=ifft(fft(Y(:,i)).*Hk);
end
figure;
[R_delta,q]=meshgrid(x_r,y_a);
mesh(x_r,y_a,abs(Y));
view(0,90);
xlabel('距离向');
% ylabel('合成孔径 (慢时间), meters');
title('方位压缩后');

%%距离徙动  由于是正侧视 所以距离走动为0 只有距离弯曲
N=6; %插值点数
f=(-0.5:1/length(tm):0.5)*PRF;%frequency range
r=(-0.5:1/length(tq):0.5)*(c/(2*fs));
r=r+rc;
delta_r=(lambda^2).*(f.^2)'*r./(8*v^2);%%range-doppler curving
for m=1:length(tm)
    for n=N/2:length(tq)
        for k=-N/2+1:N/2
            cu=floor(delta_r(m,n)*fs/c);
            du=delta_r(m,n)-cu;
            if n+cu+k>length(tq)
                S_shift(m,n)=S_shift(m,n)+S_shift(m,length(tq))*sinc(cu+du-k);
            else 
                S_shift(m,n)=S_shift(m,n)+S_shift(m,n+cu+k)*sinc(du-k);
            end       
        end
    end
end
fftshift(S_shift);
HK=repmat(Hk,1,length(tq));
Z=ifft(S_shift.*HK);
figure;
mesh(x_r,y_a,abs(Z));
view(0,90);



