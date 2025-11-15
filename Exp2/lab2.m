%%%%%%%%%%%%%%%LFM脉冲去斜处理和匀转台目标成像%%%%%%%%%%%
clear;clc;close all

%% 雷达参数
c = 3e8;
fc = 10e9;   %载频
BWf = 400e6; %chirp信号带宽
T1 = .8e-6;  %脉冲宽度
PRF = 100;   %脉冲重频
T = 2;       %合成孔径时间
Kchirp = BWf/T1; %调频率

%% 目标参数
v=0;a1=0;a2=0;         %径向速度、加速度、加加速度
w=0.01;b1=0;b2=0;      %角速度、角加速度、角加加速度
Xo = 0;Yo = 5000 ;     %位置坐标
Xsize = 40;Ysize = 60;  %目标尺寸
Ro = sqrt(Xo^2+Yo^2);   %参考距离

%% 距离维参数
fs = 2*BWf;                  %采样频率
M = round((T1+2*Xsize/c)*fs);  %距离维采样点数

%% 方位维参数
T2 = 1/PRF;                %脉冲重复周期
N = floor(T/T2);           %方位维采样点数

%% 采样及坐标参数
dt = 1/fs;                           %采样时间间隔
t = (-M/2*dt:dt:dt*(M/2-1) )+ 2*Ro/c;%快时间 
slow_t = 2*Ro/c:T2:2*Ro/c+(N-1)*T2;  %慢时间
fs1=-fs/2:fs/M:fs/2-fs/M;            %距离维fft后的频域坐标
fs2=-PRF/2:PRF/N:PRF/2-PRF/N;        %方位维fft后的频域坐标
X=-c/2*fs1/Kchirp;                   %距离向距离坐标
Y=fs2*c/fc/(2*w);                     %方位向距离坐标
Tetw = w*slow_t+(0.5*b1)*(slow_t.^2)+(1/6*b2)*(slow_t.^3); %旋转角度(rad)

%% 加载散射点模型
load Fighter
figure;
plot(-Xc,Yc,'o','MarkerSize',8,'MarkerFaceColor',[0,0,1]);grid;
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');
axis([-35 35 -30 30]);
xlabel('X(m)');ylabel('Y(m)');

%% 回波模型
Es = zeros(N,M);        %原始回波
for i=1:N
   s = zeros(1,M);
   for ii = 1:length(Xc)
       dy=Xc(ii)*sin(Tetw(i))+Yc(ii)*cos(Tetw(i));
       dx=Xc(ii)*cos(Tetw(i))-Yc(ii)*sin(Tetw(i));
       y = Yo+dy;            %纵向距离
       x = Xo+dx;            %横向距离
       R = sqrt(y.^2+x.^2);
       ap= exp(1i*2*pi*(fc*(t-2*R/c)+0.5*Kchirp*(t-2*R/c).^2));
       s =s + ap;
   end
   Es(i,:) = s;
end

%% dechirp处理
% 构造参考信号（以场景中心为参考）
t_ref = t - 2*Ro/c;  % 以参考距离为中心
ref_signal = exp(1i*2*pi*(fc*t_ref + 0.5*Kchirp*t_ref.^2));
Es_if=Es.*conj(ref_signal);
Ess=fftshift(fft(Es_if,[],2),2);


% RVP补偿
% 计算RVP项
Sc = exp(-1i*3*pi*fs.^2/Kchirp);
S_if=Ess.*conj(Sc);
S_final=fftshift(fft(S_if,[],1),1);



%% 补偿前距离像
Ess_before = fftshift(fft(Es_if,[],2),2);

%% 补偿后距离像
Ess_after = S_if;

%% 补偿前二维成像
image_2d_before = fftshift(fft(Ess_before,[],1),1);
%% 补偿后二维成像
image_2d_after = fftshift(fft(Ess_after,[],1),1);

%% 图像对比
Ess_before_normalized = 20*log10(abs(Ess_before) /max(max(abs(Ess_before(:)))));
Ess_after_normalized = 20*log(abs(Ess_after) /max( max(abs(Ess_after(:)))));
image_2d_before_normalized = 20*log(abs(image_2d_before) / max(abs(image_2d_before(:))));
image_2d_after_normalized = 20*log(abs(image_2d_after) / max(abs(image_2d_after(:))));
figure;
subplot(2,2,1);
imagesc(X, Y/100+0.5, (Ess_before_normalized));
xlabel('距离 (m)');
ylabel('时间(s)');
title('补偿前距离像');
clim([-30, 0]);
colorbar;

subplot(2,2,2);
imagesc(X, Y/100+0.5, (Ess_after_normalized));
xlabel('距离 (m)');
ylabel('时间(s)');
title('补偿后距离像');
clim([-50, 0]);
colorbar;

subplot(2,2,3);
imagesc(X, Y, (image_2d_before_normalized));
xlabel('距离 (m)');
ylabel('纵向距离(m)');
title('补偿前二维成像');
clim([-30, 0]);
colorbar;

subplot(2,2,4);
imagesc(X, Y, (image_2d_after_normalized));
xlabel('距离 (m)');
ylabel('纵向距离(m)');
title('补偿后二维成像');
clim([-30, 0]);
colorbar;
%% 错误示范
error_1=fft(Ess_after,[],1);
error_2=fftshift(fft(Ess_after,[],1),2);
error_3=fftshift(fftshift(fft(Ess_after,[],1),1),2);
error_1_normalized = 20*log10(abs(error_1) / max(abs(error_1(:))));
error_2_normalized = 20*log10(abs(error_2) / max(abs(error_2(:))));
error_3_normalized = 20*log10(abs(error_3) / max(abs(error_3(:))));
figure;
subplot(2,2,1);
imagesc(X, Y, error_1_normalized);
xlabel('距离 (m)');
ylabel('纵向距离(m)');
title('错误1');
clim([-30, 0]);
colorbar;

subplot(2,2,2);
imagesc(X, Y, (error_2_normalized));
xlabel('距离 (m)');
ylabel('纵向距离(m)');
title('错误2');
clim([-30, 0]);
colorbar;

subplot(2,2,3);
imagesc(X, Y, error_3_normalized);
xlabel('距离 (m)');
ylabel('纵向距离(m)');
title('错误3');
clim([-30, 0]);
colorbar;
