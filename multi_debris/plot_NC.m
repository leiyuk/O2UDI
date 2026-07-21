% load('true_10_year.mat');
KL_sanduPQ=zeros(1,DAY);
KL_sanduQP=zeros(1,DAY);
JS_sandu_smooth=zeros(1,DAY);


%计算JS散度
for i=1:DAY
    P_KL=latter_real(1,:,i);
    Q_KL=movmean(latter_tuiyan(1,:,i),10);
    for j=1:(hi_fenbu-lo_fenbu)/dR_fenbu
        if(P_KL(j)~=0)
            KL_sanduPQ(i)=KL_sanduPQ(i)+P_KL(j)*log(P_KL(j)/Q_KL(j));
            JS_sandu_smooth(i)=JS_sandu_smooth(i)+0.5*P_KL(j)*log(P_KL(j)*2/(P_KL(j)+Q_KL(j)));
        end
        if(Q_KL(j)~=0)
            KL_sanduQP(i)=KL_sanduQP(i)+Q_KL(j)*log(Q_KL(j)/P_KL(j));
            JS_sandu_smooth(i)=JS_sandu_smooth(i)+0.5*Q_KL(j)*log(Q_KL(j)*2/(P_KL(j)+Q_KL(j)));
        end
    end
    if(Px_real(i,1)~=0)
        KL_sanduPQ(i)=KL_sanduPQ(i)+Px_real(i,1)*log(Px_real(i,1)/Px_tuiyan(i,1));
        JS_sandu_smooth(i)=JS_sandu_smooth(i)+0.5*Px_real(i,1)*log(Px_real(i,1)*2/(Px_real(i,1)+Px_tuiyan(i,1)));
    end
    if(Px_real(i,2)~=0)
        KL_sanduPQ(i)=KL_sanduPQ(i)+Px_real(i,2)*log(Px_real(i,2)/Px_tuiyan(i,2));
        JS_sandu_smooth(i)=JS_sandu_smooth(i)+0.5*Px_real(i,2)*log(Px_real(i,2)*2/(Px_real(i,2)+Px_tuiyan(i,2)));
    end
    if(Px_tuiyan(i,1)~=0)
        KL_sanduQP(i)=KL_sanduQP(i)+Px_tuiyan(i,1)*log(Px_tuiyan(i,1)/Px_real(i,1));
        JS_sandu_smooth(i)=JS_sandu_smooth(i)+0.5*Px_tuiyan(i,1)*log(Px_tuiyan(i,1)*2/(Px_real(i,1)+Px_tuiyan(i,1)));
    end
    if(Px_tuiyan(i,2)~=0)
        KL_sanduQP(i)=KL_sanduQP(i)+Px_tuiyan(i,2)*log(Px_tuiyan(i,2)/Px_real(i,2));
        JS_sandu_smooth(i)=JS_sandu_smooth(i)+0.5*Px_tuiyan(i,2)*log(Px_tuiyan(i,2)*2/(Px_real(i,2)+Px_tuiyan(i,2)));
    end
end





set(groot, 'defaultAxesFontName', 'Times New Roman');
% ==== 上三下⼀布局 ====
figure('Color','w','Position',[50 50 1200 1000],'WindowStyle','normal');
tl = tiledlayout(3,3,'Padding','compact','TileSpacing','compact');


snap_idx = [round(DAY/3),round(DAY/3*2),DAY];

labels = {'A','B','C','D','E','F','G'};  % 子图标签（无括号）

% ===== 上排三张子图：三个时间快照 =====
for k = 1:3
    pic_num = snap_idx(k);
    ax = nexttile(k);

    plot(latter_real(1,:,pic_num), latter_real(2,:,pic_num), 'green', 'LineWidth', 1); hold on;
    plot(latter_tuiyan(1,:,pic_num), latter_tuiyan(2,:,pic_num), 'magenta', 'LineWidth', 1); hold off;

    ylabel('Orbital radius (km)');
    xlabel('Probability per bin');
    legend('True distribution','Inferred distribution','Location','northeast');

    % 标题：进化时间 + JS散度
    title(['Date: ', get_target_date(pic_num)]);

    set(ax,'FontSize',10,'LineWidth',1);

    % 轴范围：第一个子图更宽，其余更窄（沿用你的逻辑）

    axis([0 0.015 lo_fenbu/1e3 hi_fenbu/1e3]); 

    text(ax, -0.06, 1.12, labels{k}, 'Units','normalized', ...
         'HorizontalAlignment','left', 'VerticalAlignment','top', ...
         'FontName','Times New Roman','FontSize',15,'FontWeight','bold');

end

for k = 1:3
    pic_num = snap_idx(k);
    ax = nexttile(k+3);

    plot(latter_real(1,:,pic_num), latter_real(2,:,pic_num), 'green', 'LineWidth', 1); hold on;
    plot(movmean(latter_tuiyan(1,:,pic_num),10), latter_tuiyan(2,:,pic_num), 'magenta', 'LineWidth', 1); hold off;

    ylabel('Orbital radius (km)');
    xlabel('Probability per bin');
    legend('True distribution','Smoothed inferred distribution','Location','northeast');

    % 标题：进化时间 + JS散度
    title(['Date: ', get_target_date(pic_num)]);

    set(ax,'FontSize',10,'LineWidth',1);

    % 轴范围：第一个子图更宽，其余更窄（沿用你的逻辑）

    axis([0 0.015 lo_fenbu/1e3 hi_fenbu/1e3]); 

    text(ax, -0.06, 1.12, labels{k+3}, 'Units','normalized', ...
         'HorizontalAlignment','left', 'VerticalAlignment','top', ...
         'FontName','Times New Roman','FontSize',15,'FontWeight','bold');

end


% ===== 下排一张子图（横跨整行）：JS随时间 =====
ax4 = nexttile(7,[1 3]);
t_days = 1:DAY;
plot(t_days, JS_sandu,'Color', [0, 82, 155]/255, 'LineWidth', 1); hold on;
plot(t_days, JS_sandu_smooth,'Color', [230, 90, 13]/255, 'LineWidth', 1); hold on;
xlabel('Time elapsed since 2009-05-10 (days)');
ylabel('JS divergence');
legend('Original inferred','Smoothed inferred','Location','southeast');
ylim([0 0.07]);
set(ax4,'FontSize',12,'LineWidth',1);

text(ax4, -0.02, 1.12, labels{7}, 'Units','normalized', ...
     'HorizontalAlignment','left', 'VerticalAlignment','top', ...
     'FontName','Times New Roman','FontSize',15,'FontWeight','bold');




% 输出平均JS
disp(mean(JS_sandu_smooth));
disp(mean(JS_sandu_smooth-JS_sandu))

drawnow;  % 刷新图像

% ==== 精确导出为 A4 尺寸、400 DPI ====
% A4: 8.27 × 11.69 inch


outfile = fullfile(pwd,'fig3.png');  % 路径可自改
exportgraphics(tl, outfile, 'Resolution', 600, 'BackgroundColor','white');












