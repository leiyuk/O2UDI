



rng(3);
% singlemain
% clf;

set(groot, 'defaultAxesFontName', 'Times New Roman');
% ==== 上三下⼀布局 ====
figure('Color','w','Position',[50 50 1200 1000],'WindowStyle','normal');
tl = tiledlayout(3,3,'Padding','compact','TileSpacing','compact');

nums = randperm(length(linexy_record), 9);

% ===== 上排三张子图：三个时间快照 =====
for k = 1:9

    ax = nexttile(k);

    scatter(pointxy_record{nums(k)}(1,:),pointxy_record{nums(k)}(2,:),8, 'green', 'filled')
    hold on;
    
    plot(linexy_record{nums(k)}(1,:),linexy_record{nums(k)}(2,:),'Color','magenta','LineWidth',1);
    hold off;

    grid(ax, 'on');
    box(ax, 'on');
    

    text(0.05, 0.9, sprintf('$R^2 = %.4f$', R2_record{nums(k)}), ...
    'Units', 'normalized', ...
    'Interpreter', 'latex', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 12, ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k', ...
    'LineWidth', 0.6, ...
    'Margin', 3, ...
    'Clipping', 'on');

    xlabel('$B^*\;(R_{\mathrm{E}}^{-1})$', 'Interpreter', 'latex');
    ylabel('$\overline{\dot{a}}\;(\mathrm{km/s})$', ...
       'Interpreter', 'latex');

    xlim([0 1.05*max(pointxy_record{nums(k)}(1,:))]);  % 将X轴下限设为0，上限自动
    ylim([0 1.05*max(pointxy_record{nums(k)}(2,:))]);  % 将Y轴下限设为0，上限自动
   

  
    set(ax,'FontSize',10,'LineWidth',1);

end

drawnow;  % 刷新图像

% ==== 精确导出为 A4 尺寸、400 DPI ====
% A4: 8.27 × 11.69 inch


outfile = fullfile(pwd,'fig2.png');  % 路径可自改
exportgraphics(tl, outfile, 'Resolution', 600, 'BackgroundColor','white');












