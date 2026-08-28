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
snap_idx = [round(DAY/3),round(DAY/3*2),DAY];

% ==== 第一张图：六个径向分布子图 ====
fig_dist = figure('Color','w','Position',[50 50 1200 500],'WindowStyle','normal');
tl_dist = tiledlayout(fig_dist,2,3,'Padding','loose','TileSpacing','loose');
top_axes = gobjects(1,3);
bottom_axes = gobjects(1,3);

% ===== 上排三张子图：三个时间快照 =====
for k = 1:3
    pic_num = snap_idx(k);
    ax = nexttile(tl_dist,k);
    top_axes(k) = ax;

    plot(latter_real(1,:,pic_num), latter_real(2,:,pic_num), 'green', 'LineWidth', 1); hold on;
    plot(latter_tuiyan(1,:,pic_num), latter_tuiyan(2,:,pic_num), 'magenta', 'LineWidth', 1); hold off;

    ylabel('Orbital radius (km)');
    xlabel('Probability per bin');
    legend('True distribution','Inferred distribution','Location','northeast');

    % 标题：进化时间 + JS散度
    title(['Date: ', get_target_date(pic_num)],'FontSize',12,'FontWeight','normal');

    set(ax,'FontSize',11,'LineWidth',1);
    ax.XGrid = 'on';
    ax.YGrid = 'on';

    % 轴范围：第一个子图更宽，其余更窄（沿用你的逻辑）

    axis([0 0.015 lo_fenbu/1e3 hi_fenbu/1e3]); 

end

for k = 1:3
    pic_num = snap_idx(k);
    ax = nexttile(tl_dist,k+3);
    bottom_axes(k) = ax;

    plot(latter_real(1,:,pic_num), latter_real(2,:,pic_num), 'green', 'LineWidth', 1); hold on;
    plot(movmean(latter_tuiyan(1,:,pic_num),10), latter_tuiyan(2,:,pic_num), 'magenta', 'LineWidth', 1); hold off;

    ylabel('Orbital radius (km)');
    xlabel('Probability per bin');
    legend('True distribution','Smoothed inferred distribution','Location','northeast');

    % 标题：进化时间 + JS散度
    title(['Date: ', get_target_date(pic_num)],'FontSize',12,'FontWeight','normal');

    set(ax,'FontSize',12,'LineWidth',1);
    ax.XGrid = 'on';
    ax.YGrid = 'on';

    % 轴范围：第一个子图更宽，其余更窄（沿用你的逻辑）

    axis([0 0.015 lo_fenbu/1e3 hi_fenbu/1e3]); 

end

drawnow;

top_pos = top_axes(2).Position;
bottom_pos = bottom_axes(2).Position;

% annotation(fig_dist,'textbox', ...
%     [0.25,min(top_pos(2)+top_pos(4)+0.062,0.955),0.5,0.03], ...
%     'String','(a) Original inferred distributions', ...
%     'HorizontalAlignment','center','VerticalAlignment','middle', ...
%     'FontName','Times New Roman','FontSize',15, ...
%     'EdgeColor','none');
% 
% annotation(fig_dist,'textbox', ...
%     [0.25,bottom_pos(2)+bottom_pos(4)+0.032,0.5,0.03], ...
%     'String','(b) Smoothed inferred distributions', ...
%     'HorizontalAlignment','center','VerticalAlignment','middle', ...
%     'FontName','Times New Roman','FontSize',15, ...
%     'EdgeColor','none');

outfile_dist = fullfile(pwd,'fig_simulation_radial_distribution.png');
exportgraphics(fig_dist, outfile_dist, 'Resolution', 600, 'BackgroundColor','white');

% ==== 第二张图：JS散度随时间变化 ====
fig_js = figure('Color','w','Position',[100 100 1000 400],'WindowStyle','normal');
ax_js = axes(fig_js);
t_days = 1:DAY;
plot(ax_js,t_days, JS_sandu,'Color', [0, 82, 155]/255, 'LineWidth', 1); hold(ax_js,'on');
plot(ax_js,t_days, JS_sandu_smooth,'Color', [230, 90, 13]/255, 'LineWidth', 1); hold(ax_js,'off');
xlabel('Time elapsed since 2009-05-10 (days)');
ylabel('JS divergence');
legend('Original inferred','Smoothed inferred','Location','southeast');
ylim([0 0.07]);
set(ax_js,'FontSize',16,'LineWidth',1);
ax_js.XGrid = 'off';
ax_js.YGrid = 'on';

% text(ax_js, -0.04, 1.08, 'G', 'Units','normalized', ...
%      'HorizontalAlignment','left', 'VerticalAlignment','top', ...
%      'FontName','Times New Roman','FontSize',16,'FontWeight','bold');

% 输出平均JS
disp(mean(JS_sandu_smooth));
disp(mean(JS_sandu_smooth-JS_sandu))

drawnow;  % 刷新图像

outfile_js = fullfile(pwd,'fig_simulation_js_divergence.png');
exportgraphics(fig_js, outfile_js, 'Resolution', 600, 'BackgroundColor','white');











