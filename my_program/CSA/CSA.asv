close all;
clear;
clc;
R_nc = 20e3;
Vr = 150;
Tr = 2.5e-6;
Kr = 20e12;
f0 = 5.3e9;
Bw_dop = 80;
Fr = 60e6;
Fa = 200;
Naz = 1024;
Nrg = 320;
sita_r_c = (3.5*pi)/180;
c = 3e8;

R0 = R_nc*cos(sita_r_c);
Nr = Tr*Fr;
Bw_range = Kr*Tr;
lambda = c/f0;
fnc = 2*sin(sita_r_c)*Vr/lambda;
La_real = 0.886*2*Vr*cos(sita_r_c)/Bw_dop;
beta_w = 0.886*lambda/La_real; %实际阵列单程1 实际阵列双程0.88 合成孔径阵列双程0.44
La = beta_w*R0;
a_sr = Fr / Bw_range;
a_sa = Fa / Bw_dop;
Mamb = round(fnc/Fa);%多普勒模糊？？
NFFT_r = Nrg;
NFFT_a = Naz;

R_ref = R0;
fn_ref = fnc;

delta_R0 = 0;
delta_R1 = 120;%目标1和目标2方位向相差
delta_R2 = 80;%目标2和目标3距离向相差
%目标1
x1 = R0;
y1 = delta_R0+x1*tan(sita_r_c);
%目标2
x2 = x1;
y2 = y1+delta_R1;
%目标3
x3 = x2 + delta_R2;
y3 = y2 + delta_R2*tan(sita_r_c);%??
x_range = [x1, x2, x3];
y_azimuth = [y1,y2,y3];
n_tr = zeros(size(x_range));
%计算各自三个目标的中心穿越时刻
for i=1:length(x_range)
    n_tr(i) = (y_azimuth(i)-x_range(i)*tan(sita_r_c))/Vr;
end

% 距离(方位)向时间

% 距离
tr = 2*R0/c + (-Nrg/2:(Nrg/2-1))/Fr;%时间行向量
fr = (-NFFT_r/2:NFFT_r/2-1)*(Fr/NFFT_r);

% 方位
ta = (-Naz/2:Naz/2-1)/Fa;
fa = fnc + fftshift(-NFFT_a/2:NFFT_a/2-1)*(Fa/NFFT_a);%将零频置在起始位置

% 生成距离时间矩阵 key------------------
tr_mtx = ones(Naz,1)*tr;
ta_mtx = ta'*ones(1,Nrg);
fr_mtx = ones(Naz,1)*fr;
fa_mtx = fa'*ones(1,Nrg);

% 生成点目标原始数据
s_echo = zeros(Naz,Nrg); 

for i=1:length(x_range)
    R_n = ((x_range(i)^2*ones(Naz,Nrg))+(Vr*ta_mtx-y_azimuth(i)).^2).^0.5;
    sita = atan((Vr*(ta_mtx-n_tr(i)))./x_range(i));
    azimuth_w1 = (sinc(0.886*sita*beta_w)).^2;
    azimuth_wind = (abs(ta_mtx-n_tr(i))<=1.135*La/(2*Vr));
    azimuth_wind = azimuth_wind.*azimuth_w1;
    range_wind = (abs(tr_mtx-2*R_n./c)<=Tr/2);
    s_echo = s_echo+azimuth_wind.*range_wind.*exp(-4j*pi*R_n./lambda).*exp(1j*pi*Kr*(tr_mtx-2*R_n./c).^2);
end

figure;
subplot(2,2,1);
imagesc(real(s_echo));
xlabel('距离向');
ylabel('方位向');
title('实部信号');
text(300,-70,'图1,原始数据');

subplot(2,2,2);
imagesc(imag(s_echo));
xlabel('距离向');
ylabel('方位向');
title('虚部信号');

subplot(2,2,3);
imagesc(abs(s_echo));
xlabel('距离向');
ylabel('方位向');
title('幅度信号');

subplot(2,2,4);
imagesc(angle(s_echo));
xlabel('距离向');
ylabel('方位向');
title('相位');
figure;
[x_r,y_a]=meshgrid(tr,ta);
mesh(x_r,y_a,abs(s_echo));
view(0,90);
xlabel('距离向');
% ylabel('合成孔径 (慢时间), meters');
title('方位压缩后、、、、、、、');

% 频域谱
figure;
subplot(2,1,1);
echo_rd = fft(s_echo);
imagesc(abs(echo_rd));
xlabel('距离频域');
ylabel('方位频域');
title('信号   多普勒距离');

subplot(2,1,2);
echo_2f = fft(echo_rd,[],2);
imagesc(abs(echo_2f));
xlabel('距离频域');
ylabel('方位频域');
title('信号二维频域');

% 变换到距离多普勒域
s_rd = s_echo.*exp(-2j*pi*fnc*(ta_mtx));%????????  频谱搬移
S_rd = fft(s_rd,NFFT_a,1);
figure;
imagesc(abs(S_rd));
D_vt = sqrt(1-lambda.^2*(fa_mtx).^2./(4*Vr^2));
D_vt_r = sqrt(1-lambda.^2*(fa').^2./(4*Vr^2));
%D_fn_ref = 1;
D_fn_ref = sqrt(1-lambda.^2*(fnc).^2./(4*Vr^2));
K_src = 2*Vr^2*f0^3.*D_vt.^3./(c*R_ref*fa_mtx.^2);
K_m = Kr./(1-Kr./K_src);
s_cs = exp(1j*pi*K_m.*(D_fn_ref./D_vt-1).*(tr_mtx-2*R_ref./(c*D_vt)).^2);
S_rd1 = S_rd.*s_cs;

figure;
subplot(2,1,1);
imagesc(abs(S_rd));
title('未进行补余RCMC距离多普勒域');

subplot(2,1,2);
imagesc(abs(S_rd1));
title('进行了补余RCMC距离多普勒域');


%进行距离压缩和距离徙动校正
%s_src = exp(1j*pi*fr_mtx.^2./(K_m.*(1+1./D_vt)));%距离压缩
s_src = exp(1j*pi*fr_mtx.^2.*D_vt./(K_m*D_fn_ref));%距离压缩
%s_rm = exp(1j*4*pi*fr_mtx*R_ref./D_vt);%距离徙动
s_rm = exp(1j*4*pi*fr_mtx*R_ref.*(1./D_vt-1/D_fn_ref)./c);
s_cs2 = fftshift(s_src.*s_rm,2); %距离零频在两端 需要平移 
S_2f = fft(S_rd1,NFFT_r,2);
S2_2f = S_2f.*s_cs2; % 校正之后的二维频域

%画出二维频域频谱和相位图
figure;
subplot(2,2,1);
imagesc(abs(S2_2f));
title('二维频域幅度谱');

subplot(2,2,2);
imagesc(angle(S2_2f));
title('二维频域相位图');


%进行距离反变换
S2_rd = ifft(S2_2f,NFFT_r,2);
subplot(2,2,3);
imagesc(abs(S2_rd));
title('距离压缩和距离徙动校正之后:距离多普勒幅度');

subplot(2,2,4);
imagesc(angle(S2_rd));
title('距离压缩和距离徙动校正之后:距离多普勒相位');

% 相位补余和方位压缩

R0_rcmc = c*tr./2;
S_amf = exp(1j*4*pi.*(D_vt_r*R0_rcmc).*f0./c);%乘的先后顺序会影响结果精度
S_phase = exp(-1j*4*pi*K_m.*(1-D_vt./D_fn_ref).*(R0_rcmc./D_vt-R_ref./D_vt).^2/(c^2));
H2 = S_amf.*S_phase;
S3_rd = S2_rd.*H2; %相位补偿和方位压缩距离-多普勒
S3_2t = ifft(S3_rd,NFFT_a,1); % 二维时域

figure;
imagesc(abs(S3_rd));
xlabel('距离向');
ylabel('方位向');
title('距离多普勒域域幅度信号');

figure;
imagesc(abs(S3_2t));
xlabel('距离向');
ylabel('方位向');
title('二维时域幅度信号');

%% 计算目标1的方位中心位置
target1_x = rem(R0*tan(sita_r_c)*Fa/Vr,Naz);
if target1_x > Naz/2
    target1_x = target1_x - Naz/2;
else 
    target1_x = target1_x + Naz/2;
end
target1_x = round(target1_x);
% 计算目标2的方位
target2_x = target1_x + delta_R1*Fa/Vr;

% 计算目标3的方位
target3_x = target2_x + delta_R2*tan(sita_r_c)*Fa/Vr;

%% 计算距离向位置

target1_y = round(Nrg/2+1)+(R_ref/D_fn_ref*2-2*R0)*Fr/c;

target2_y = target1_y;

target3_y = target2_y + delta_R2*2*Fr/c;













