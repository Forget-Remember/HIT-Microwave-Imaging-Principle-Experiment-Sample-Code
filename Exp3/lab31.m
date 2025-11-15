%%%%%%%%%%%%%%%%%%参数化运动补偿（估计目标运动轨迹）%%%%%%%%%%%%%%

clear;close all;clc
%% 雷达参数
burst = 128; %一个脉冲串内的脉冲数129
pulses =128;   %脉冲串数128
c = 3e8;
f0 = 10e9;    %脉冲初始频率
bw = 128e6;   %频带带宽
T1 = (burst-1)/bw; % 脉冲宽度
PRF = 20e3;   %脉冲重频（两个脉冲之间，不是脉冲串之间！）
T2 = 1/PRF;   %脉冲重复周期
dr = c/(2*bw);%距离分辨率

%% 目标参数
W = 0.03;  % 角速度
Vr = 70; % 径向速度
ar = 0.1;  % 径向加速度
R0 = 16e3; % 目标初始距离

%% figure 1 散射点分布
figure(1)
load Fighter 
plot(-Xc,Yc,'o', 'MarkerSize',8,'MarkerFaceColor', [0,0,1]);grid;
set(gca,'FontName', 'Arial', 'FontSize',12,'FontWeight','Bold'); 
axis([-35 35 -30 30])
xlabel('X [m]'); ylabel('Y [m]');
title('散射点分布');
%极坐标表示的散射中心
[theta,r]=cart2pol(Xc,Yc);

i = 1:pulses*burst;
T = T1/2+2*R0/c+(i-1)*T2;  %时间刻度,以每个脉冲的中间时刻+回波的延迟时间
Rvr = Vr*T+(0.5*ar)*(T.^2);%纵向位移 
Tetw = W*T;                %旋转位移(rad)

i = 1:burst;
df = (i-1)*1/T1;      %频率增量(每一个脉冲相对于f0的脉冲增量）
k = (2*pi*(f0+df))/c; %组脉冲波束向量（可以在这写2*pi，然后下面再乘个双程的2，但不应该叫波束了）
k_fac=ones(pulses,1)*k;%波束二维向量

%% 计算后向散射电场 
Es = zeros(pulses,burst);
for scat = 1:1:length(Xc)
    arg = (Tetw - theta(scat) );
    rngterm = R0 + Rvr - r(scat)*sin(arg);
    range = reshape(rngterm,burst,pulses);
    range = range.';
    phase = k_fac.* range * 2;
    Ess = exp(-1i*phase);% 单个散射点的后向散射电场
    Es = Es+Ess;% 39个散射点加一起的后向散射电场
end

%% figure 2 补偿前一维距离像
figure(2)
Es_range = fftshift(fft(Es,[],2),2); % 距离维FFT
imagesc(abs(Es_range));
title('补偿前一维距离像');
xlabel('距离单元'); ylabel('脉冲数');
colorbar;

%% figure 3 补偿前二维成像
figure(3)
Es_2d =fftshift( fft(Es_range,[],1),1); % 方位维FFT
imagesc(abs(Es_2d));
title('补偿前二维ISAR图像');
xlabel('距离单元'); ylabel('方位单元');
colorbar;

%% 一次性补偿
RP = fft(Es,[],2); % 距离像

for l=1:burst; % Cross-correlation between RPn & RPref
    %cr(l,:) = abs(xcorr(abs(RP(l,:)),(abs(RP(1,:)))));
    cr(l,:) = abs(ifft(fft(abs(RP(1,:))).* fft(abs(conj(RP(l,:))))));
    pk(l) = find((max(cr(l,:))== cr(l,:))); % Find max. ind. (range shift) range)
end

Spk = smooth((0:pulses-1),pk,0.8,'rlowess'); % smoothing the delays
RangeShifts = dr*pk; % range shifts
SmRangeShifts = dr*Spk; % range shifts
RangeDif = SmRangeShifts(2:pulses)-SmRangeShifts(1:pulses-1); % range differences

%% 低阶多项式拟合
pulse_time = (0:pulses-1) * T2; % 时间向量
% 一阶多项式拟合（线性拟合）
p1 = polyfit(pulse_time, SmRangeShifts, 1);
fitted_shifts_1 = polyval(p1, pulse_time);

% 二阶多项式拟合
p2 = polyfit(pulse_time, SmRangeShifts, 2);
fitted_shifts_2 = polyval(p2, pulse_time);

%% figure 4 低阶多项式拟合结果
figure(4)
plot(pulse_time, RangeShifts, 'bo', 'MarkerSize', 4, 'DisplayName', '原始距离偏移');
hold on;
plot(pulse_time, SmRangeShifts, 'r-', 'LineWidth', 2, 'DisplayName', '平滑后距离偏移');
plot(pulse_time, fitted_shifts_1, 'g--', 'LineWidth', 2, 'DisplayName', '一阶多项式拟合');
plot(pulse_time, fitted_shifts_2, 'm-.', 'LineWidth', 2, 'DisplayName', '二阶多项式拟合');
hold off;
grid on;
xlabel('时间 (s)');
ylabel('距离偏移 (m)');
title('低阶多项式拟合结果');
legend('show');
set(gca, 'FontSize', 12);



RangeDifAv = mean(RangeDif); % average range differences
T_burst = T(pulses+1)-T(1); % time between the bursts
Vr_Dif = (-RangeDif/T_burst); % estimated radial velocity from each RP
Vr_av = (-RangeDifAv/T_burst); % estimated radial velocity (average)

f = (f0+df); % frequency vector
T_reshape = reshape(T,pulses,burst); % prepare time matrix
F = f.'*ones(1,burst); % prepare frequency matrix
Es_comp = Es.*exp((1i*4*pi*F/c).*(Vr_av*T_reshape)); % Phase of E-field is compensated

%% 补偿后的处理
Es_comp_range = fft(Es_comp,[],2); % 补偿后的距离维FFT
%% figure 5 补偿后一维距离像
figure(5)
imagesc(abs(Es_comp_range));
title('补偿后一维距离像');
xlabel('距离单元'); ylabel('脉冲数');
colorbar;

%% figure 6 补偿后二维ISAR图像
figure(6)
compensated_2d = fftshift(fft(Es_comp_range,[],1),1); % 方位维FFT
Es_shifted = circshift(compensated_2d, -40);
imagesc(abs(Es_shifted));
title('补偿后二维ISAR图像');
xlabel('距离单元'); ylabel('方位单元');
colorbar;

%% 显示坐标轴
X = -dr*((burst)/2-1):dr:dr*burst/2;
t_slow = T(1:burst:burst*pulses);
PRF1 = 1/(T2*(pulses+1));
fs2 = -PRF1/2:PRF1/pulses:PRF1/2-PRF1/pulses;
Y = fs2*c/f0/(2*W);

%% figure 7 补偿前后对比（正确坐标）
figure(7)
subplot(1,2,1)
imagesc(X, Y, abs(Es_2d));
title('补偿前二维ISAR图像');
xlabel('距离(m)'); ylabel('方位(m)');
axis xy; colorbar;

subplot(1,2,2)
imagesc(X, Y, abs(Es_shifted));
title('补偿后二维ISAR图像');
xlabel('距离(m)'); ylabel('方位(m)');
axis xy; colorbar;


