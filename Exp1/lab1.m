%%%%%%%%%%%%%%SFCW回波生成和一维距离像序列观测%%%%%%%%%%%%%%

clear;close all;clc
%% 雷达参数
pulses = 128; %一个脉冲串内的脉冲数129
burst =128;   %脉冲串数128
c = 3e8;
f0 = 10e9;    %脉冲初始频率
bw = 128e6;   %频带带宽
T1 = (pulses-1)/bw; % 脉冲宽度
PRF = 20e3;   %脉冲重频（两个脉冲之间，不是脉冲串之间！）
T2 = 1/PRF;   %脉冲重复周期
dr = c/(2*bw);%距离分辨率

%% 目标参数
W = 0;  % 角速度
Vr = 70;   % 径向速度
ar = 0.1;  % 径向加速度
R0 = 16e3; % 目标初始距离

%% figure 1 散射点分布
load Fighter 
plot(-Xc,Yc,'o', 'MarkerSize',8,'MarkerFaceColor', [0,0,1]);grid;
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight','Bold'); 
axis([-35 35 -30 30])
xlabel('X [m]'); ylabel('Y [m]');
%极坐标表示的散射中心
[theta,r]=cart2pol(Xc,Yc);

i = 1:pulses*burst;
T = T1/2+2*R0/c+(i-1)*T2;  %时间刻度,以每个脉冲的中间时刻+回波的延迟时间
Rvr = -Vr*T-(0.5*ar)*(T.^2);%纵向位移 
Tetw = W*T;                %旋转位移(rad)

i = 1:pulses;
df = (i-1)*1/T1;      %频率增量(每一个脉冲相对于f0的脉冲增量）
k = (2*pi*(f0+df))/c;
k_fac=ones(burst,1)*k;%波束二维向量

%% 计算后向散射电场 
Es = zeros(burst,pulses);
for scat = 1:1:length(Xc)
    arg = (Tetw - theta(scat) );
    rngterm = R0 + Rvr - r(scat)*sin(arg);
    range = reshape(rngterm,pulses,burst);
    range = range.';
    phase = k_fac.* range * 2;
    Ess = exp(-1i*phase);% 单个散射点的后向散射电场
    Es = Es+Ess;% 39个散射点加一起的后向散射电场
end

%% figure 2 一维距离像
% 对每个脉冲进行FFT得到距离像
range_profile = fftshift(fft(Es, [], 2),2);  % 对每一行(每个脉冲)做FFT

% 距离轴
range_axis = (0:pulses-1) * dr;

%% 绘制一维距离像
figure;
subplot(2,1,1);
imagesc(range_axis, 1:burst, abs(range_profile));
xlabel('距离 (m)');
ylabel('脉冲序号');
title('一维距离像序列');
colorbar;

subplot(2,1,2);
% 选择某个特定脉冲的距离像进行显示
pulse_idx = 64;  % 选择第64个脉冲
plot(range_axis, abs(range_profile(pulse_idx, :)));
xlabel('距离 (m)');
ylabel('幅度');
title(['第', num2str(pulse_idx), '个脉冲的一维距离像']);


%% figure 3 相邻距离像对比
% 选择某个特定脉冲的距离像进行显示
figure;
subplot(2,2,1);
pulse_idx = 20;  % 选择第51个脉冲
plot(range_axis, abs(range_profile(pulse_idx, :)));
xlabel('距离 (m)');
ylabel('幅度');
title(['第', num2str(pulse_idx), '个脉冲的一维距离像']);
grid on;
subplot(2,2,2);
pulse_idx = 21;  % 选择第52个脉冲
plot(range_axis, abs(range_profile(pulse_idx, :)));
xlabel('距离 (m)');
ylabel('幅度');
title(['第', num2str(pulse_idx), '个脉冲的一维距离像']);
grid on;
subplot(2,2,3);
pulse_idx = 22;  % 选择第53个脉冲
plot(range_axis, abs(range_profile(pulse_idx, :)));
xlabel('距离 (m)');
ylabel('幅度');
title(['第', num2str(pulse_idx), '个脉冲的一维距离像']);
grid on;
subplot(2,2,4);
pulse_idx = 23;  % 选择第54个脉冲
plot(range_axis, abs(range_profile(pulse_idx, :)));
xlabel('距离 (m)');
ylabel('幅度');
title(['第', num2str(pulse_idx), '个脉冲的一维距离像']);
grid on;


%% figure 4 相距较远的距离像对比 
% 选择某个特定脉冲的距离像进行显示
figure;
subplot(2,2,1);
pulse_idx = 20;  % 选择第20个脉冲
plot(range_axis, abs(range_profile(pulse_idx, :)));
xlabel('距离 (m)');
ylabel('幅度');
title(['第', num2str(pulse_idx), '个脉冲的一维距离像']);
grid on;
subplot(2,2,2);
pulse_idx = 40;  % 选择第40个脉冲
plot(range_axis, abs(range_profile(pulse_idx, :)));
xlabel('距离 (m)');
ylabel('幅度');
title(['第', num2str(pulse_idx), '个脉冲的一维距离像']);
grid on;
subplot(2,2,3);
pulse_idx = 60;  % 选择第60个脉冲
plot(range_axis, abs(range_profile(pulse_idx, :)));
xlabel('距离 (m)');
ylabel('幅度');
title(['第', num2str(pulse_idx), '个脉冲的一维距离像']);
grid on;
subplot(2,2,4);
pulse_idx = 80;  % 选择第80个脉冲
plot(range_axis, abs(range_profile(pulse_idx, :)));
xlabel('距离 (m)');
ylabel('幅度');
title(['第', num2str(pulse_idx), '个脉冲的一维距离像']);
grid on;

