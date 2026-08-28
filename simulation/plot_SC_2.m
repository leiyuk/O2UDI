set(groot, 'defaultAxesFontName', 'Times New Roman');
% ==== 上三下⼀布局 ====
figure('Color','w','Position',[50 50 800 1000],'WindowStyle','normal');
tl = tiledlayout(4,2,'Padding','compact','TileSpacing','compact');
miu=398600e9;
R=6378e3;

% 计算三个快照在原循环中的索引（跳过第一个=1）
snap_idx = 1:ns/4:ns+1;   % [1, 1+ns/3, 1+2ns/3, ns+1]
snap_idx = snap_idx(2:5); % 只取后三个

% labels = {'A','B','C','D','E','F','G','H'};  % 子图标签（无括号）

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
    num_true=(k2(index)+k3(index)+k4(index));
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
    x1 = aX( 1:num_true)';  y1 = eX( 1:num_true)';
    x2 = aY( 1:num_tuiyan)'; y2 = eY( 1:num_tuiyan)';
    
    % --- 2. 手动计算全局最大密度 (解决输出参数报错) ---
    % 我们不需要每个点的密度，只需要知道两组数据中“最高峰”是多少
    [f1, ~] = ksdensity([x1, y1], [x1, y1]);
    [f2, ~] = ksdensity([x2, y2], [x2, y2]);
    max_density = max([max(f1), max(f2)]);


    % 真实值
    ax = nexttile(2*k-1);
%     x_plot = X(1, 1:num_true)'; % 转置为列向量
%     y_plot = X(2, 1:num_true)'; % 转置为列向量
    x_plot = aX(1:num_true)'; % 转置为列向量
    y_plot = eX(1:num_true)'; % 转置为列向量    
    % --- 计算密度 ---
    % 这里的输入必须是 [num_true x 2] 的矩阵
    [f, ~] = ksdensity([x_plot, y_plot], [x_plot, y_plot]);    
    % --- 绘图 ---
    % 此时 x_plot, y_plot, f 均为 num_true x 1，维度完全匹配
    scatter(x_plot/1e3, y_plot, 2, f, 'filled'); 
    cb = colorbar; % 每一行四张子图共用该颜色范围，colorbar 保持在中间
    cb.FontName = 'Times New Roman';
    cb.FontSize = 8;
    colormap(viridis);
    clim([0, max_density]); % 统一颜色范围
    if k==4
    xlabel('Semimajor axis (km)');
%     ylabel('Eccentricity');
    end
    ylabel('Eccentricity');
    set(ax,'FontSize',10,'LineWidth',1,'XGrid','on','YGrid','on','Box','on');
    axis([6.6e3 7.8e3 0 0.1]);

%     text(ax, -0.075,1.25, labels{2*k-1}, 'Units','normalized', ...
%          'HorizontalAlignment','left', 'VerticalAlignment','top', ...
%          'FontName','Times New Roman','FontSize',12,'FontWeight','bold');
    subtitle(ax,['$t = $', num2str(round((pic_num-1)*Ts/86400)), ' days'],"FontSize",10,'Interpreter', 'latex');
    % 推演值
    ax = nexttile(2*k);
%     x_plot = Y(1,1:num_tuiyan)'; % 转置为列向量
%     y_plot = Y(2,1:num_tuiyan)'; % 转置为列向量
    x_plot = aY(1:num_tuiyan)'; % 转置为列向量
    y_plot = eY(1:num_tuiyan)'; % 转置为列向量
    % --- 计算密度 ---
    % 这里的输入必须是 [num_true x 2] 的矩阵
    [f, ~] = ksdensity([x_plot, y_plot], [x_plot, y_plot]);    
    % --- 绘图 ---
    % 此时 x_plot, y_plot, f 均为 num_true x 1，维度完全匹配
    scatter(x_plot/1e3, y_plot, 2, f, 'filled'); 
%     colorbar;
    colormap(viridis);
    clim([0, max_density]); % 统一颜色范围
    if k==4
    xlabel('Semimajor axis (km)');
%     ylabel('Eccentricity');
    end
%     subtitle(['evolution time: ', num2str(round((pic_num-1)*Ts/86400)), ' days']);
    set(ax,'FontSize',10,'LineWidth',1,'XGrid','on','YGrid','on','Box','on');
    axis([6.6e3 7.8e3 0 0.1]);

%     text(ax, -0.075,1.25, labels{2*k}, 'Units','normalized', ...
%          'HorizontalAlignment','left', 'VerticalAlignment','top', ...
%          'FontName','Times New Roman','FontSize',12,'FontWeight','bold');
    subtitle(ax,['$t = $', num2str(round((pic_num-1)*Ts/86400)), ' days'],"FontSize",10,'Interpreter', 'latex');

end

ax = nexttile(1);
title("(a) Reference population","FontSize",14)

ax = nexttile(2);
title("(b) Inferred population","FontSize",14)

drawnow;  % 刷新图像

% ==== 精确导出为 A4 尺寸、400 DPI ====
% A4: 8.27 × 11.69 inch


outfile = fullfile(pwd,'fig1_2.png');  % 路径可自改
exportgraphics(tl, outfile, 'Resolution', 1200, 'BackgroundColor','white');





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






