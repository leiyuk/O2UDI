%以下为Lc=0.02m,0.05m,0.15m,0.3m时y=delta_v的分布（准确值）
clear;
x=logspace(-2,2.6,1000);
%0.02m
Lc3=0.02;
namudac3=log10(Lc3);
miu4=0.2*(-0.3-1.4*(namudac3+1.75))+1.85;
sigma4=sqrt((0.2+0.1333*(namudac3+3.5))^2*0.04+0.16);
y3=1./x/log(10)/sqrt(2*pi)/sigma4.*exp(-(log10(x)-miu4).^2/2/sigma4^2);
max3=10^(miu4-sigma4^2*log(10));%极大值点的横坐标
Ex3=10^(miu4+sigma4^2*log(10)/2);%均值
middle3=fsolve(@(t) integral(@(x) 1./x/log(10)/sqrt(2*pi)/sigma4.*exp(-(log10(x)-miu4).^2/2/sigma4^2),0,t)-0.5,100);%中位数
%0.05m
Lc1=0.05;
namudac1=log10(Lc1);
miu3=-0.28*namudac1+1.3;
sigma3=sqrt(0.16+0.04*(0.1333*namudac1+0.66655)^2);
y1=1./x/log(10)/sqrt(2*pi)/sigma3.*exp(-(log10(x)-miu3).^2/2/sigma3^2);
max1=10^(miu3-sigma3^2*log(10));%极大值点的横坐标
Ex1=10^(miu3+sigma3^2*log(10)/2);%均值
middle1=fsolve(@(t) integral(@(x) 1./x/log(10)/sqrt(2*pi)/sigma3.*exp(-(log10(x)-miu3).^2/2/sigma3^2),0,t)-0.5,100);%中位数
%0.15m
Lc2=0.15;
namudac2=log10(Lc2);
alfa=0.3+0.4*(namudac2+1.2);
miu1=1.85+0.2*(-0.6-0.318*(namudac2+1.1));
sigma1=sqrt(0.16+0.04*(0.1+0.2*(namudac2+1.3))^2);
miu2=1.85+0.2*(-1.2);
sigma2=sqrt(0.16+0.04*0.5^2);
y2=alfa*1./x/log(10)/sqrt(2*pi)/sigma1.*exp(-(log10(x)-miu1).^2/2/sigma1^2)...
    +(1-alfa)*1./x/log(10)/sqrt(2*pi)/sigma2.*exp(-(log10(x)-miu2).^2/2/sigma2^2);
%0.3m
Lc4=0.3;
namudac4=log10(Lc4);
alfa=0.3+0.4*(namudac4+1.2);
miu5=(-0.6-0.318*(namudac4+1.1))*0.2+1.85;
sigma5=sqrt((0.1+0.2*(namudac4+1.3))^2+0.16);
miu6=(-1.2-1.333*(namudac4+0.7))*0.2+1.85;
sigma6=sqrt(0.5^2*0.04+0.16);
y4=alfa*1./x/log(10)/sqrt(2*pi)/sigma5.*exp(-(log10(x)-miu5).^2/2/sigma5^2)...
    +(1-alfa)*1./x/log(10)/sqrt(2*pi)/sigma6.*exp(-(log10(x)-miu6).^2/2/sigma6^2);
Ex4=alfa*10^(miu5+sigma5^2*log(10)/2)+(1-alfa)*10^(miu6+sigma6^2*log(10)/2);%均值

set(groot, 'defaultAxesFontName', 'Times New Roman');
figure('Color','w','Position',[50 50 1200 1000],'WindowStyle','normal');
f1=figure(1);

plot(x,y3,x,y1,x,y2,x,y4,'LineWidth',2);
xlabel('Magnitude of ejection-velocity increment (m/s)');
ylabel('Probability density');
legend('$L_c = 0.02\;\mathrm{m}$', ...
        '$L_c = 0.05\;\mathrm{m}$', ...
        '$L_c = 0.15\;\mathrm{m}$', ...
        '$L_c = 0.30\;\mathrm{m}$', ...
        'Interpreter', 'latex','Location','northeast');
xlim([0 400]);
ylim([0 0.016]);
ax = gca(f1); % 获取 f1 中当前的 Axes 句柄
set(ax, 'FontSize', 20,'LineWidth',1); % 对 Axes 对象设置 FontSize

drawnow;  % 刷新图像

% ==== 精确导出为 A4 尺寸、400 DPI ====
% A4: 8.27 × 11.69 inch


outfile = fullfile(pwd,'delta_v.png');  % 路径可自改
exportgraphics(f1, outfile, 'Resolution', 600, 'BackgroundColor','white');