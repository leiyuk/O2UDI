%以下为Lc=0.02m,0.05m,0.15m,0.3m时y=S/m的分布
clear;
format long
x=logspace(-5,0,1000);
%0.02m
Lc3=0.02;
namudac3=log10(Lc3);
miu4=-0.3-1.4*(namudac3+1.75);
sigma4=0.2+0.1333*(namudac3+3.5);
y3=1./x/log(10)/sqrt(2*pi)/sigma4.*exp(-(log10(x)-miu4).^2/2/sigma4^2);
max3=10^(miu4-sigma4^2*log(10));%极大值点的横坐标
maxf3=1/sqrt(2*pi)/log(10)/sigma4*exp(1/2*sigma4^2*log(10)^2-miu4*log(10));%概率密度函数的最大值
Ex3=10^(miu4+sigma4^2*log(10)/2);%均值
middle3=fsolve(@(t) integral(@(x) 1./x/log(10)/sqrt(2*pi)/sigma4.*exp(-(log10(x)-miu4).^2/2/sigma4^2),0,t)-0.5,0.1);%中位数
%0.05m
Lc1=0.05;
namudac1=log10(Lc1);
miu3=-0.3-1.4*(namudac1+1.75);
sigma3=0.2+0.1333*(namudac1+3.5);
y1=1./x/log(10)/sqrt(2*pi)/sigma3.*exp(-(log10(x)-miu3).^2/2/sigma3^2);
max1=10^(miu3-sigma3^2*log(10));%极大值点的横坐标
maxf1=1/sqrt(2*pi)/log(10)/sigma3*exp(1/2*sigma3^2*log(10)^2-miu3*log(10));%概率密度函数的最大值
Ex1=10^(miu3+sigma3^2*log(10)/2);%均值
middle1=fsolve(@(t) integral(@(x) 1./x/log(10)/sqrt(2*pi)/sigma3.*exp(-(log10(x)-miu3).^2/2/sigma3^2),0,t)-0.5,0.1);%中位数
%0.15m
Lc2=0.15;
namudac2=log10(Lc2);
alfa=0.3+0.4*(namudac2+1.2);
miu1=-0.6-0.318*(namudac2+1.1);
sigma1=0.1+0.2*(namudac2+1.3);
miu2=-1.2;
sigma2=0.5;
y2=alfa*1./x/log(10)/sqrt(2*pi)/sigma1.*exp(-(log10(x)-miu1).^2/2/sigma1^2)...
    +(1-alfa)*1./x/log(10)/sqrt(2*pi)/sigma2.*exp(-(log10(x)-miu2).^2/2/sigma2^2);
Ex2=alfa*10^(miu1+sigma1^2*log(10)/2)+(1-alfa)*10^(miu2+sigma2^2*log(10)/2);%均值
middle2=fsolve(@(t) integral(@(x) alfa*1./x/log(10)/sqrt(2*pi)/sigma1.*exp(-(log10(x)-miu1).^2/2/sigma1^2)...
    +(1-alfa)*1./x/log(10)/sqrt(2*pi)/sigma2.*exp(-(log10(x)-miu2).^2/2/sigma2^2),0,t)-0.5,0.2);%中位数
%0.3m
Lc4=0.3;
namudac4=log10(Lc4);
alfa=0.3+0.4*(namudac4+1.2);
miu5=-0.6-0.318*(namudac4+1.1);
sigma5=0.1+0.2*(namudac4+1.3);
miu6=-1.2-1.333*(namudac4+0.7);
sigma6=0.5;
y4=alfa*1./x/log(10)/sqrt(2*pi)/sigma5.*exp(-(log10(x)-miu5).^2/2/sigma5^2)...
    +(1-alfa)*1./x/log(10)/sqrt(2*pi)/sigma6.*exp(-(log10(x)-miu6).^2/2/sigma6^2);
Ex4=alfa*10^(miu5+sigma5^2*log(10)/2)+(1-alfa)*10^(miu6+sigma6^2*log(10)/2);%均值
middle4=fsolve(@(t) integral(@(x) alfa*1./x/log(10)/sqrt(2*pi)/sigma5.*exp(-(log10(x)-miu5).^2/2/sigma5^2)...
    +(1-alfa)*1./x/log(10)/sqrt(2*pi)/sigma6.*exp(-(log10(x)-miu6).^2/2/sigma6^2),0,t)-0.5,0.1);%中位数

set(groot, 'defaultAxesFontName', 'Times New Roman');
figure('Color','w','Position',[50 50 1200 1000],'WindowStyle','normal');
f1=figure(1);

plot(x,y3,x,y1,x,y2,x,y4,'LineWidth',2);
xlabel('Area-to-mass ratio ($\mathrm{m^2/kg}$)','Interpreter', 'latex');
ylabel('Probability density');
legend('$L_c = 0.02\;\mathrm{m}$', ...
        '$L_c = 0.05\;\mathrm{m}$', ...
        '$L_c = 0.15\;\mathrm{m}$', ...
        '$L_c = 0.30\;\mathrm{m}$', ...
        'Interpreter', 'latex','Location','northeast');
xlim([0 1]);
ylim([0 8]);
ax = gca(f1); % 获取 f1 中当前的 Axes 句柄
set(ax, 'FontSize', 20,'LineWidth',1); % 对 Axes 对象设置 FontSize

drawnow;  % 刷新图像

% ==== 精确导出为 A4 尺寸、400 DPI ====
% A4: 8.27 × 11.69 inch


outfile = fullfile(pwd,'S_m.png');  % 路径可自改
exportgraphics(f1, outfile, 'Resolution', 600, 'BackgroundColor','white');






