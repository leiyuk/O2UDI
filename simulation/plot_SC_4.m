% 原来.mat文件中计算的JS散度有问题，因此需要重新计算这部分
KL_sanduPQ=zeros(1,ns+1);
KL_sanduQP=zeros(1,ns+1);
JS_sandu=zeros(1,ns+1);
for i=1:ns+1
    P_KL=xiao_cankao(1,:,i);
    Q_KL=xiao_tuiyan(1,:,i);
    for j=1:(hi_fenbu-lo_fenbu)/dR_fenbu
        if(P_KL(j)~=0)
            KL_sanduPQ(i)=KL_sanduPQ(i)+P_KL(j)*log(P_KL(j)/Q_KL(j));
            JS_sandu(i)=JS_sandu(i)+0.5*P_KL(j)*log(P_KL(j)*2/(P_KL(j)+Q_KL(j)));
        end
        if(Q_KL(j)~=0)
            KL_sanduQP(i)=KL_sanduQP(i)+Q_KL(j)*log(Q_KL(j)/P_KL(j));
            JS_sandu(i)=JS_sandu(i)+0.5*Q_KL(j)*log(Q_KL(j)*2/(P_KL(j)+Q_KL(j)));
        end
    end
    if(Px_cankao(i,1)~=0)
        KL_sanduPQ(i)=KL_sanduPQ(i)+Px_cankao(i,1)*log(Px_cankao(i,1)/Px_tuiyan(i,1));
        JS_sandu(i)=JS_sandu(i)+0.5*Px_cankao(i,1)*log(Px_cankao(i,1)*2/(Px_cankao(i,1)+Px_tuiyan(i,1)));
    end
    if(Px_cankao(i,2)~=0)
        KL_sanduPQ(i)=KL_sanduPQ(i)+Px_cankao(i,2)*log(Px_cankao(i,2)/Px_tuiyan(i,2));
        JS_sandu(i)=JS_sandu(i)+0.5*Px_cankao(i,2)*log(Px_cankao(i,2)*2/(Px_cankao(i,2)+Px_tuiyan(i,2)));
    end
    if(Px_tuiyan(i,1)~=0)
        KL_sanduQP(i)=KL_sanduQP(i)+Px_tuiyan(i,1)*log(Px_tuiyan(i,1)/Px_cankao(i,1));
        JS_sandu(i)=JS_sandu(i)+0.5*Px_tuiyan(i,1)*log(Px_tuiyan(i,1)*2/(Px_cankao(i,1)+Px_tuiyan(i,1)));
    end
    if(Px_tuiyan(i,2)~=0)
        KL_sanduQP(i)=KL_sanduQP(i)+Px_tuiyan(i,2)*log(Px_tuiyan(i,2)/Px_cankao(i,2));
        JS_sandu(i)=JS_sandu(i)+0.5*Px_tuiyan(i,2)*log(Px_tuiyan(i,2)*2/(Px_cankao(i,2)+Px_tuiyan(i,2)));
    end
end





set(groot, 'defaultAxesFontName', 'Times New Roman');
% ==== 图1：三个时间快照 ====
figure('Color','w','Position',[50 50 1200 420],'WindowStyle','normal');
tl_dist = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

% 计算三个快照在原循环中的索引（跳过第一个=1）
snap_idx = 1:ns/3:ns+1;   % [1, 1+ns/3, 1+2ns/3, ns+1]
snap_idx = snap_idx(2:4); % 只取后三个

% labels = {'A','B','C'};  % 子图标签（无括号）

% ===== 三张子图：三个时间快照 =====
for k = 1:3
    pic_num = snap_idx(k);
    ax = nexttile(tl_dist,k);

    plot(xiao_cankao(1,:,pic_num), xiao_cankao(2,:,pic_num), 'green', 'LineWidth', 1); hold on;
    plot(xiao_tuiyan(1,:,pic_num), xiao_tuiyan(2,:,pic_num), 'magenta', 'LineWidth', 1); hold off;

    ylabel('Orbital radius (km)');
    xlabel('Probability per bin');
    legend('Reference distribution','Inferred distribution','Location','northeast');

    % 标题：进化时间 + JS散度
    title(['$t = $', num2str(round((pic_num-1)*Ts/86400)), ' days'], ...
        'Interpreter', 'latex','FontSize',14);

    set(ax,'FontSize',12,'LineWidth',1,'XGrid','on','YGrid','on');

    % 轴范围：第一个子图更宽，其余更窄（沿用你的逻辑）

    axis([0 0.005 lo_fenbu/1e3  hi_fenbu/1e3]);

%     text(ax, -0.08, 1.12, labels{k}, 'Units','normalized', ...
%          'HorizontalAlignment','left', 'VerticalAlignment','top', ...
%          'FontName','Times New Roman','FontSize',16,'FontWeight','bold');

end




% ===== 图2：JS散度随时间 =====
figure('Color','w','Position',[100 100 720 450],'WindowStyle','normal');
tl_js = tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
ax4 = nexttile(tl_js,1);
t_days = (0:ns)*Ts/86400;
plot(t_days, JS_sandu,'Color', [0, 82, 155]/255, 'LineWidth', 1);
xlabel('Time (days)');
ylabel('JS divergence');
ylim([0 0.05]);
set(ax4,'FontSize',15,'LineWidth',1,'XGrid','off','YGrid','on');




% 输出平均JS
disp(mean(JS_sandu));

drawnow;  % 刷新图像

% ==== 精确导出为 A4 尺寸、400 DPI ====
% A4: 8.27 × 11.69 inch


outfile_dist = fullfile(pwd,'fig1_4_distribution.png');  % 路径可自改
exportgraphics(tl_dist, outfile_dist, 'Resolution', 600, 'BackgroundColor','white');

outfile_js = fullfile(pwd,'fig1_4_js.png');  % 路径可自改
exportgraphics(tl_js, outfile_js, 'Resolution', 600, 'BackgroundColor','white');












