%绘制Gabbard Plot

clc;
clear;


miu=398600; %单位km
R=6378;%地球半径，单位Km

% 指定读取的COSMOS碎片的两行根数的文件夹路径
folderPath = '../multi_debris/cosmos_2251';
% 获取文件夹中所有 .txt 文件的信息
filePattern = fullfile(folderPath, '*.txt');
txtFiles = dir(filePattern);
% 初始化一个存储数据的单元数组（或矩阵，视情况而定）
allData = {};
NUM=length(txtFiles);  %绘制的碎片数量
X=zeros(3,NUM);

% 循环读取每个文件
for k = 1:length(txtFiles)
    % 获取文件的完整路径
    baseFileName = txtFiles(k).name;
    fullFileName = fullfile(folderPath, baseFileName);
    fidin = fopen(fullFileName, 'r');
    data='';
    N=0;
    
    uuu=0;
    while ~feof(fidin)        %判断是否文件读完
        tline=fgetl(fidin);  %逐行读入数据
        if isempty(tline)
            continue
        end
        data=[data tline];
        uuu=uuu+1;
        if uuu>10
            break;
        end
    end
    fclose(fidin);

    %提取数据
    e=str2double(data(-42+155*1:-36+155*1));
    e=e*10^(-7);%偏心率
    %M=str2double(data(-25+155*1:-18+155*1));%平近点角
    n=str2double(data(-16+155*1:-6+155*1)); %平均角速度，圈/天
    a=(86400^2*miu/(4*pi^2*n^2))^(1/3); %半长轴，单位km
    %E=fsolve(@(E)E-e*sin(E)-M,M);
    %f=2*atan(sqrt((1+e)/(1-e))*tan(E/2));
    p=a*(1-e^2);
    ra=p/(1-e);
    rp=p/(1+e);
    X(1,k)=24*60/n;         %将周期保存到矩阵第1行中，单位min
    X(2,k)=ra-R;          %将远地点高度存入矩阵的第2行
    X(3,k)=rp-R;          %将近地点高度存入矩阵的第3行
end

set(groot, 'defaultAxesFontName', 'Times New Roman');
figure('Color','w','Position',[50 50 1200 1000],'WindowStyle','normal');

f1=figure(1);
scatter(X(1,:),X(2,:),10,"filled",'Color', [0, 82, 155]/255);
hold on;
scatter(X(1,:),X(3,:),10,"filled",'Color', [230, 90, 13]/255);
hold off;
legend('Apogee','Perigee','Location','northwest');
xlabel('Obital period (minutes)');
ylabel('Altitude (km)');
xlim([97 111]);
ylim([400 1800]);

ax = gca(f1); % 获取 f1 中当前的 Axes 句柄
set(ax, 'FontSize', 20,'LineWidth',1); % 对 Axes 对象设置 FontSize

drawnow;  % 刷新图像

% ==== 精确导出为 A4 尺寸、400 DPI ====
% A4: 8.27 × 11.69 inch


outfile = fullfile(pwd,'gabbard.png');  % 路径可自改
exportgraphics(f1, outfile, 'Resolution', 600, 'BackgroundColor','white');