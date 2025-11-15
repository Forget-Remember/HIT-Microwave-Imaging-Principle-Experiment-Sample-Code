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

for l=1:burst % Cross-correlation between RPn & RPref
    %cr(l,:) = abs(xcorr(abs(RP(l,:)),(abs(RP(1,:)))));
    cr(l,:) = abs(ifft(fft(abs(RP(1,:))).* fft(abs(conj(RP(l,:))))));
    pk(l) = find((max(cr(l,:))== cr(l,:))); % Find max. ind. (range shift) range)
end
pk_new=pk-pk(1);
Spk = smooth((0:pulses-1),pk_new,0.8,'rlowess'); % smoothing the delays
RangeShifts = dr*pk; % range shifts
SmRangeShifts = dr*Spk; % range shifts
RangeDif = SmRangeShifts(2:pulses)-SmRangeShifts(1:pulses-1); % range differences

%% figure 4 包络对齐距离偏移估计
figure(4)
pulse_time = (0:pulses-1) * T2;
plot(pulse_time, RangeShifts, 'bo-', 'MarkerSize', 4, 'DisplayName', '原始距离偏移');
hold on;
plot(pulse_time, SmRangeShifts, 'r-', 'LineWidth', 2, 'DisplayName', '平滑后距离偏移');
hold off;
grid on;
xlabel('时间 (s)');
ylabel('距离偏移 (m)');
title('包络对齐距离偏移估计');
legend('show');
set(gca, 'FontSize', 12);

Ne = (0:pulses-1);
Es_align = Es.*exp(-1i*2*pi*Spk.*Ne/pulses);
RP1=abs(fftshift(ifft(Es_align,[],2),2));

%% figure 5 包络对齐后一维距离像
figure(5)
Es_align_range = fftshift(fft(Es_align,[],2),2);
imagesc(abs(Es_align_range));
title('包络对齐后一维距离像');
xlabel('距离单元'); ylabel('脉冲数');
colorbar;

%% figure 6 包络对齐后二维ISAR图像
figure(6)
Es_align_2d = fftshift(fft(Es_align_range,[],1),1);
imagesc(abs(Es_align_2d));
title('包络对齐后二维ISAR图像');
xlabel('距离单元'); ylabel('方位单元');
colorbar;

%% 特显点选择
for i=1:pulses
    Un(i)=mean(abs(RP1(:,i)));
    Un2(i)=mean(abs(RP1(:,i).^2));
end
xigma=1-Un.^2./Un2;
point=find(min(xigma)==xigma);

%% figure 7 特显点选择图
figure(7)
subplot(2,1,1)
plot(Un, 'b-', 'LineWidth', 2);
hold on;
plot(point, Un(point), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold off;
grid on;
xlabel('距离单元');
ylabel('均值');
title('距离像均值');
legend('均值', '特显点位置');

subplot(2,1,2)
plot(xigma, 'g-', 'LineWidth', 2);
hold on;
plot(point, xigma(point), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold off;
grid on;
xlabel('距离单元');
ylabel('方差');
title('距离像方差（最小值对应特显点）');
legend('方差', '特显点位置');

fprintf('选择的特显点位置: 距离单元 %d\n', point);

%% 相位校正
phase=RP1(:,point)./abs(RP1(:,point));
% 重新构建相位校正后的信号
Es_phase_comp = Es_align;
for n=1:burst
    Es_phase_comp(n,:) = Es_phase_comp(n,:) ./ phase(n);
end

%% figure 8 相位校正后一维距离像
figure(8)
Es_phase_range = fftshift(fft(Es_phase_comp,[],2),2);
imagesc(abs(Es_phase_range));
title('相位校正后一维距离像');
xlabel('距离单元'); ylabel('脉冲数');
colorbar;

%% figure 9 相位校正后二维ISAR图像
figure(9)
Es_phase_2d_raw = fftshift(fft(Es_phase_range,[],1),1);
Es_phase_2d=circshift(Es_phase_2d_raw, -35,2);
imagesc(abs(Es_phase_2d));
title('相位校正后二维ISAR图像');
xlabel('距离单元'); ylabel('方位单元');
colorbar;

%% 显示坐标轴
X = -dr*((burst)/2-1):dr:dr*burst/2;
t_slow = T(1:burst:burst*pulses);
PRF1 = 1/(T2*(pulses+1));
fs2 = -PRF1/2:PRF1/pulses:PRF1/2-PRF1/pulses;
Y = fs2*c/f0/(2*W);

%% figure 10 补偿前后对比（正确坐标）
figure(10)
subplot(2,2,1)
imagesc(X, Y, abs(Es_2d));
title('补偿前二维ISAR图像');
xlabel('距离(m)'); ylabel('方位(m)');
axis xy; colorbar;

subplot(2,2,2)
imagesc(X, Y, abs(Es_align_2d));
title('包络对齐后二维ISAR图像');
xlabel('距离(m)'); ylabel('方位(m)');
axis xy; colorbar;

subplot(2,2,3)
imagesc(X, Y, abs(Es_phase_2d));
title('相位校正后二维ISAR图像');
xlabel('距离(m)'); ylabel('方位(m)');
axis xy; colorbar;

subplot(2,2,4)
% 最终补偿结果（包络对齐+相位校正）
Es_final_2d = Es_phase_2d;
imagesc(X, Y, abs(Es_final_2d));
title('最终补偿后二维ISAR图像');
xlabel('距离(m)'); ylabel('方位(m)');
axis xy; colorbar;

%% figure 11 一维距离像对比
figure(11)
subplot(2,2,1)
plot(abs(Es_range(64,:)));
title('补偿前一维距离像（第64个脉冲）');
xlabel('距离单元'); ylabel('幅度');
grid on;

subplot(2,2,2)
plot(abs(Es_align_range(64,:)));
title('包络对齐后一维距离像（第64个脉冲）');
xlabel('距离单元'); ylabel('幅度');
grid on;

subplot(2,2,3)
plot(abs(Es_phase_range(64,:)));
title('相位校正后一维距离像（第64个脉冲）');
xlabel('距离单元'); ylabel('幅度');
grid on;

subplot(2,2,4)
plot(angle(phase));
title('相位校正量');
xlabel('脉冲数'); ylabel('相位 (rad)');
grid on;

