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
% ==== 上三下⼀布局 ====
tl=figure('Color','w','Position',[50 50 1200 1000],'WindowStyle','normal');



t_days = (0:ns)*Ts/86400;
plot(t_days, JS_sandu,'Color', [0, 82, 155]/255, 'LineWidth', 1); 
xlabel('evolution time (day)');
ylabel('JS divergence');
ylim([0 0.05]);
% 设置所有坐标轴的默认字体大小
set(0, 'DefaultAxesFontSize', 20);


% 设置所有绘图线条（plot 出来的线）的默认粗细
set(0, 'DefaultLineLineWidth', 1);






drawnow;  % 刷新图像

% ==== 精确导出为 A4 尺寸、400 DPI ====
% A4: 8.27 × 11.69 inch


outfile = fullfile(pwd,'fig1_1.png');  % 路径可自改
exportgraphics(tl, outfile, 'Resolution', 600, 'BackgroundColor','white');












