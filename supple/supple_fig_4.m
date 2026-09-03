clc;
clear;
% 1. 设置包含CSV文件的文件夹路径 (请替换为你的实际路径)
folderPath = 'ESA'; 

% 2. 获取文件夹下所有以 .csv 结尾的文件信息
filePattern = fullfile(folderPath, '*.csv');
csvFiles = dir(filePattern);
numFiles = length(csvFiles);

if numFiles == 0
    error('未在指定文件夹中找到CSV文件，请检查路径。');
end

% 3. 循环读取每个文件并动态生成同名变量
for i = 1:numFiles
    currentFileName = csvFiles(i).name;
    fullFilePath = fullfile(folderPath, currentFileName);
    
    % 获取去掉 '.csv' 后缀的文件名
    [~, nameWithoutExt, ~] = fileparts(currentFileName);
    
    % 将文件名转换为合法的 MATLAB 变量名
    % (例如将 "1-data test" 转换为 "x1_data_test")
    validVarName = matlab.lang.makeValidName(nameWithoutExt);
    
    try
        % 读取 F 列的数据
        columnData = readmatrix(fullFilePath, 'Range', 'F:F');
        
        % 清理因表头产生的 NaN 值（可选）
        columnData(isnan(columnData)) = [];
        
        % 提取前47个数据
        if length(columnData) >= 47
            finalData = columnData(1:47);
        else
            warning('文件 %s 的有效数据不足 47 个，已用 NaN 填充。', currentFileName);
            finalData = NaN(47, 1);
            finalData(1:length(columnData)) = columnData;
        end
        
        % ================= 核心步骤 =================
        % assignin 函数会在当前工作区直接创建一个新变量
        % 第一个参数 'caller' 表示在当前调用空间（即你的工作区）创建
        % 第二个参数是变量名，第三个参数是对应的值
        assignin('caller', validVarName, finalData);
        
        fprintf('已成功读取并创建变量: %s\n', validVarName);
        
    catch ME
        warning('读取文件 %s 时发生错误: %s', currentFileName, ME.message);
    end
end

disp('所有文件读取并转换为独立变量完成！');

set(groot, 'defaultAxesFontName', 'Times New Roman');
f1 = figure('Color','w','Position',[50 50 1200 1000], ...
    'WindowStyle','normal');

X = (1980:2026).';
fragmentationDebris = round(RF_up - RF_down + PF_up - ...
    PF_down_and_PL_up);
allDebris = round(RM_up - PF_down_and_PL_up);
barData = [fragmentationDebris(:), allDebris(:)];

b = bar(X, barData, 'grouped', 'BarWidth', 0.9);
b(1).FaceColor = [0, 82, 155]/255;
b(1).EdgeColor = 'none';
b(2).FaceColor = [230, 90, 13]/255;
b(2).EdgeColor = 'none';

legend(b, {'Fragmentation debris','All cataloged debris'}, ...
    'Location','northwest');
xlabel('Reference year');
ylabel('Cataloged debris count');
xlim([1979.5 2026.5]);
maxCount = max(barData, [], 'all', 'omitnan');
yLimitMax = max(5000, ceil(1.05 * maxCount / 5000) * 5000);
ylim([0 yLimitMax]);

% 设置坐标轴与年度刻度
ax = gca;
% 强制将 Y 轴的指数设为 0，从而取消科学计数法
ax.YAxis.Exponent = 0;
ax.XTick = [1980:5:2025];
ax.YGrid = 'on';
ax.XGrid = 'off';
ax.GridAlpha = 0.15;
set(ax, 'FontSize', 20, 'LineWidth', 1, 'Box', 'on');

drawnow;  % 刷新图像

% ==== 精确导出为 A4 尺寸、400 DPI ====
% A4: 8.27 × 11.69 inch


outfile = fullfile(pwd,'supple_fig_4.png');  % 路径可自改
exportgraphics(f1, outfile, 'Resolution', 600, 'BackgroundColor','white');
