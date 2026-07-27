set(groot, 'defaultAxesFontName', 'Times New Roman');
% ==== 上三下⼀布局 ====
figure('Color','w','Position',[50 50 1200 1000],'WindowStyle','normal');
tl = tiledlayout(4,2,'Padding','compact','TileSpacing','compact');
miu=398600e9;
R=6378e3;

% 计算三个快照在原循环中的索引（跳过第一个=1）
snap_idx = 1:ns/4:ns+1;   % [1, 1+ns/3, 1+2ns/3, ns+1]
snap_idx = snap_idx(2:5); % 只取后三个

labels = {'A','B','C','D','E','F','G','H'};  % 子图标签（无括号）

% ===== 上排三张子图：三个时间快照 =====
for k = 1:4
    pic_num = snap_idx(k);
    X=zeros(3,10000);
    aX=zeros(1,10000);
    eX=zeros(1,10000);
    % 尝试一下画第10000个周期的gabbard图
    Y=zeros(3,10000);
    aY=zeros(1,10000);
    eY=zeros(1,10000);
    
    index=pic_num;
    a_true=[ae_cankao_L2(index,1:k2(index),1),ae_cankao_L3(index,1:k3(index),1),ae_cankao_L4(index,1:k4(index),1)];
    e_true=[ae_cankao_L2(index,1:k2(index),2),ae_cankao_L3(index,1:k3(index),2),ae_cankao_L4(index,1:k4(index),2)];
    
    num_true=0;
    for j=1:(k2(index)+k3(index)+k4(index))
        if (a_true(j)>6378e3+200e3)
            num_true=num_true+1;
            a=a_true(j);
            e=e_true(j);
            aX(num_true)=a;
            eX(num_true)=e;
            X(1,num_true)=2*pi/sqrt(miu/a^3)/60;
            X(2,num_true)=(a*(1+e)-R)/1e3;
            X(3,num_true)=(a*(1-e)-R)/1e3;
        end
    end
    disp(num_true)
    disp(sum(aX)/num_true)
    disp(sum(eX)/num_true)
    
    num_tuiyan=0;
    for j=1:no_x
        zuijin=erfen(tae_xiao_tuiyan(1,1:k5(j),j),1,k5(j),(index-1)*Ts);
        if(tae_xiao_tuiyan(2,zuijin,j)>200e3+6378e3)
            num_tuiyan=num_tuiyan+1;
            a=tae_xiao_tuiyan(2,zuijin,j);
            e=tae_xiao_tuiyan(3,zuijin,j);
            aY(num_tuiyan)=a;
            eY(num_tuiyan)=e;
            Y(1,num_tuiyan)=2*pi/sqrt(miu/a^3)/60;
            Y(2,num_tuiyan)=(a*(1+e)-R)/1e3;
            Y(3,num_tuiyan)=(a*(1-e)-R)/1e3;
        end
    end
    disp(num_tuiyan)
    disp(sum(aY)/num_tuiyan)
    disp(sum(eY)/num_tuiyan)
    
    
    % --- 1. 提取并清洗数据 (确保列向量) ---
    x1 = X(1, 1:num_true)';  y1 = X(2, 1:num_true)';
    x2 = Y( 1,1:num_tuiyan)'; y2 = Y(2, 1:num_tuiyan)';
    x3 = X(1, 1:num_true)';  y3 = X(3, 1:num_true)';
    x4 = Y( 1,1:num_tuiyan)'; y4 = Y(3, 1:num_tuiyan)';
    
    % --- 2. 手动计算全局最大密度 (解决输出参数报错) ---
    % 我们不需要每个点的密度，只需要知道两组数据中“最高峰”是多少
    [f1, ~] = ksdensity([x1, y1], [x1, y1]);
    [f2, ~] = ksdensity([x2, y2], [x2, y2]);
    [f3, ~] = ksdensity([x3, y3], [x3, y3]);
    [f4, ~] = ksdensity([x4, y4], [x4, y4]);
    max_density = max([max(f1), max(f2),max(f3),max(f4)]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 真实值
    ax = nexttile(tl,2*k-1);
    axis(ax, 'off'); % 隐藏容器本身的坐标轴

    tl_child = tiledlayout(tl, 1, 1);
    tl_child.Layout.Tile = 2*k-1; % 指定子布局占用父布局的第 i 个格子
    tl_child.TileSpacing = 'compact'; % 子图之间紧凑一点

    ax_merge = nexttile(tl_child,1);

    % 远地点：单独计算 KDE 密度
    x_plot_apo = X(1, 1:num_true)';
    y_plot_apo = X(2, 1:num_true)';
    [f_apo, ~] = ksdensity([x_plot_apo, y_plot_apo], [x_plot_apo, y_plot_apo]);

    % 近地点：单独计算 KDE 密度
    x_plot_per = X(1, 1:num_true)';
    y_plot_per = X(3, 1:num_true)';
    [f_per, ~] = ksdensity([x_plot_per, y_plot_per], [x_plot_per, y_plot_per]);

    % 合并到同一张 Gabbard 图中：上三角表示远地点，下三角表示近地点
    scatter(ax_merge,x_plot_apo, y_plot_apo, 2, f_apo, '^', 'filled');
    hold(ax_merge,'on');
    scatter(ax_merge,x_plot_per, y_plot_per, 2, f_per, 'v', 'filled');
    hold(ax_merge,'off');
    colorbar(ax_merge);
    colormap(viridis);
    clim(ax_merge,[0, max_density]); % 统一颜色范围
    xlabel(ax_merge,'Orbital period (minutes)');
    ylabel(ax_merge,'Altitude (km)');
    set(ax_merge,'FontSize',8,'LineWidth',1);
    axis(ax_merge,[90 120 400 1500]);

    text(ax, -0.075,1.3, labels{2*k-1}, 'Units','normalized', ...
         'HorizontalAlignment','left', 'VerticalAlignment','top', ...
         'FontName','Times New Roman','FontSize',12,'FontWeight','bold');
    subtitle(ax,['$t = $', num2str(round((pic_num-1)*Ts/86400)), ' days'],"FontSize",10,'Interpreter', 'latex');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 推演值
    ax = nexttile(tl,2*k);
    axis(ax, 'off'); % 隐藏容器本身的坐标轴

    tl_child = tiledlayout(tl, 1, 1);
    tl_child.Layout.Tile = 2*k; % 指定子布局占用父布局的第 i 个格子
    tl_child.TileSpacing = 'compact'; % 子图之间紧凑一点

    ax_merge = nexttile(tl_child,1);

    % 远地点：单独计算 KDE 密度
    x_plot_apo = Y(1, 1:num_tuiyan)';
    y_plot_apo = Y(2, 1:num_tuiyan)';
    [f_apo, ~] = ksdensity([x_plot_apo, y_plot_apo], [x_plot_apo, y_plot_apo]);

    % 近地点：单独计算 KDE 密度
    x_plot_per = Y(1, 1:num_tuiyan)';
    y_plot_per = Y(3, 1:num_tuiyan)';
    [f_per, ~] = ksdensity([x_plot_per, y_plot_per], [x_plot_per, y_plot_per]);

    % 合并到同一张 Gabbard 图中：上三角表示远地点，下三角表示近地点
    scatter(ax_merge,x_plot_apo, y_plot_apo, 2, f_apo, '^', 'filled');
    hold(ax_merge,'on');
    scatter(ax_merge,x_plot_per, y_plot_per, 2, f_per, 'v', 'filled');
    hold(ax_merge,'off');
    colormap(viridis);
    clim(ax_merge,[0, max_density]); % 统一颜色范围
    xlabel(ax_merge,'Orbital period (minutes)');
    ylabel(ax_merge,'Altitude (km)');
    set(ax_merge,'FontSize',8,'LineWidth',1);
    axis(ax_merge,[90 120 400 1500]);

    text(ax, -0.075,1.3, labels{2*k}, 'Units','normalized', ...
         'HorizontalAlignment','left', 'VerticalAlignment','top', ...
         'FontName','Times New Roman','FontSize',12,'FontWeight','bold');
    subtitle(ax,['$t = $', num2str(round((pic_num-1)*Ts/86400)), ' days'],"FontSize",10,'Interpreter', 'latex');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end



ax = nexttile(tl,1);
title("Reference population","FontSize",14)

ax = nexttile(tl,2);
title("Inferred population","FontSize",14)

drawnow;  % 刷新图像

% ==== 精确导出为 A4 尺寸、400 DPI ====
% A4: 8.27 × 11.69 inch


outfile = fullfile(pwd,'fig1_5.png');  % 路径可自改
exportgraphics(tl, outfile, 'Resolution', 600, 'BackgroundColor','white');





function xuhao=erfen(xulie,lo,hi,zhi)%返回序号对应的值<=查找值，序号+1对应的值大于查找值
if (zhi<=xulie(lo))
    xuhao=lo;
    return
end

if (zhi>=xulie(hi))
    xuhao=hi;
    return
end

while 1
    if(xulie(floor((lo+hi)/2))<=zhi)
        lo=floor((lo+hi)/2);
    else
        hi=floor((lo+hi)/2);
    end
    if(hi-lo<=1)
        xuhao=lo;
        return;
    end
end
end






