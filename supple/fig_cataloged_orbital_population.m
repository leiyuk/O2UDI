clc;
clear;

% Folder containing the ESA CSV files.
folderPath = 'ESA';

filePattern = fullfile(folderPath, '*.csv');
csvFiles = dir(filePattern);
numFiles = length(csvFiles);

if numFiles == 0
    error('未在指定文件夹中找到CSV文件，请检查路径。');
end

% Read column F from each CSV file and create a variable named after the file.
for i = 1:numFiles
    currentFileName = csvFiles(i).name;
    fullFilePath = fullfile(folderPath, currentFileName);
    [~, nameWithoutExt, ~] = fileparts(currentFileName);
    validVarName = matlab.lang.makeValidName(nameWithoutExt);

    try
        columnData = readmatrix(fullFilePath, 'Range', 'F:F');
        columnData(isnan(columnData)) = [];

        if length(columnData) >= 47
            finalData = columnData(1:47);
        else
            warning('文件 %s 的有效数据不足 47 个，已用 NaN 填充。', ...
                currentFileName);
            finalData = NaN(47, 1);
            finalData(1:length(columnData)) = columnData;
        end

        assignin('caller', validVarName, finalData);
        fprintf('已成功读取并创建变量: %s\n', validVarName);
    catch ME
        warning('读取文件 %s 时发生错误: %s', ...
            currentFileName, ME.message);
    end
end

disp('所有文件读取并转换为独立变量完成！');

set(groot, 'defaultAxesFontName', 'Times New Roman');
f1 = figure('Color', 'w', 'Position', [50 50 1400 1300], ...
    'WindowStyle', 'normal');
t = tiledlayout(f1, 2, 1, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

X = (1980:2026).';
allObjects = round(all_up(:));
payloads = round(PF_down_and_PL_up(:));
leoObjects = round(LEO_up(:));

% (a) Comparison by object type.
ax1 = nexttile(t, 1);
b1 = bar(ax1, X, [allObjects, payloads], 'grouped', 'BarWidth', 0.9);
b1(1).FaceColor = [230, 90, 13]/255;
b1(1).EdgeColor = 'none';
b1(2).FaceColor = [0, 82, 155]/255;
b1(2).EdgeColor = 'none';

legend(ax1, b1, {'All cataloged objects', 'Payloads'}, ...
    'Location', 'northwest');
title(ax1, '(a) By object type', 'FontWeight', 'normal');
xlabel(ax1, 'Reference year');
ylabel(ax1, 'Cataloged object count');
xlim(ax1, [1979.5 2026.5]);
ylim(ax1, [0 50000]);

% (b) Comparison by orbital region.
ax2 = nexttile(t, 2);
b2 = bar(ax2, X, [allObjects, leoObjects], 'grouped', 'BarWidth', 0.9);
b2(1).FaceColor = [230, 90, 13]/255;
b2(1).EdgeColor = 'none';
b2(2).FaceColor = [0, 82, 155]/255;
b2(2).EdgeColor = 'none';

legend(ax2, b2, {'All cataloged objects', 'Objects in LEO'}, ...
    'Location', 'northwest');
title(ax2, '(b) By orbital region', 'FontWeight', 'normal');
xlabel(ax2, 'Reference year');
ylabel(ax2, 'Cataloged object count');
xlim(ax2, [1979.5 2026.5]);
ylim(ax2, [0 50000]);

% Apply consistent formatting to both panels.
for ax = [ax1, ax2]
    ax.YAxis.Exponent = 0;
    ax.XTick = [1980:5:2025];
    ax.YGrid = 'on';
    ax.XGrid = 'off';
    ax.GridAlpha = 0.15;
    set(ax, 'FontSize', 20, 'LineWidth', 1, 'Box', 'on');
end

linkaxes([ax1, ax2], 'x');
drawnow;

outfile = fullfile(pwd, 'fig_cataloged_orbital_population.png');
exportgraphics(f1, outfile, 'Resolution', 600, ...
    'BackgroundColor', 'white');
