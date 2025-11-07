clf;
huitu_jiange = 9;
pic_num = 1;
ns = DAY;

% 创建一个VideoWriter对象来保存视频
videoFileName = 'cosmos_2251.mp4';
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
    plot(latter_real(1,:,pic_num), latter_real(2,:,pic_num), 'green', 'LineWidth', 0.5);hold on;
    plot(latter_tuiyan(1,:,pic_num), latter_tuiyan(2,:,pic_num), 'magenta', 'LineWidth', 0.5);hold off;
    ylabel('orbital radius (km)');  % 设置y轴标签
    xlabel('probability density');  % 设置x轴标签
    legend('reference distribution', 'inferred distribution','Location','northeast');  % 添加图例
    title({['evolution time: ', num2str(round(pic_num)), ' day']; ...
           ['JS divergence: ', num2str(JS_sandu(pic_num), '%.4f')]});  % 设置标题
    set(gca, 'FontSize', 20);  % 设置坐标轴字体大小
    
    % 设置坐标轴范围
    axis([0 0.015 lo_fenbu/1e3 hi_fenbu/1e3]);
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
