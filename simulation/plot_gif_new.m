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






clf;
flag_tuiyan = 0;
pic_num = 1;

% 创建一个VideoWriter对象来保存视频
videoFileName = 'simulation.mp4';
videoObj = VideoWriter(videoFileName, 'MPEG-4');  % 使用MPEG-4格式并应用H.264编码
videoObj.Quality = 100;   % 设置质量为100（最高质量）
videoObj.FrameRate = 20;  % 设置帧率（可调）

% 打开视频文件以开始写入
open(videoObj);
f_tu = figure(1);  % 创建一个图形句柄

while 1
    if pic_num > ns + 1
        break;
    end
    

    % 绘制数据，xiao_cankao为参考值，xiao_tuiyan为外推值
    plot(xiao_cankao(1,:,pic_num), xiao_cankao(2,:,pic_num), 'green', 'LineWidth', 1);hold on;
    plot(xiao_tuiyan(1,:,pic_num), xiao_tuiyan(2,:,pic_num), 'magenta', 'LineWidth', 1);hold off;
    ylabel('Orbital radius (km)');  % 设置y轴标签
    xlabel('Probability per bin');  % 设置x轴标签
    legend('Reference distribution', 'Inferred distribution');  % 添加图例
    title({['$t = $', num2str(round((pic_num-1)*Ts/86400)), ' day']; ...
           ['JS divergence: ', num2str(JS_sandu(pic_num), '%.4f')]},'Interpreter', 'latex');  % 设置标题
    set(gca, 'FontSize', 20,'LineWidth',1);  % 设置坐标轴字体大小
    
    % 设置坐标轴范围
    axis([0 0.01 lo_fenbu/1e3 hi_fenbu/1e3]);
    set(f_tu, 'unit', 'normalized', 'position', [0, 0, 16*100/1920, 9*100/1080]);  % 使用16:9的比例
    
    drawnow;  % 更新图形
    
    % 获取当前图形帧
    F_tu = getframe(gcf);  
    I_tu = frame2im(F_tu);  % 将帧转换为图像
    
    % 将图像写入视频
    writeVideo(videoObj, I_tu);  
    

    
    % 更新图片编号
    pic_num = pic_num + huitu_jiange;
    clf;
end

% 完成所有帧写入后，关闭视频文件
close(videoObj);
