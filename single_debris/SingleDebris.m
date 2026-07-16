% 验证方法1：单碎片验证法程序


% clear;

% 读取文件
% 该颗卫星数据所在文件夹

function [R2,all,use,pointxy,linexy]=SingleDebris(basefolder,fileName,dt_use,num,plot_flag)

% basefolder = input('输入所在文件夹地址：', 's');
% % 打开数据txt文件
% fileName = input('输入文件名：', 's');
fullPath = fullfile(basefolder, fileName);
fidin = fopen(fullPath, 'r');
if fidin == -1
    error('数据文件打开失败：%s', fullPath);
else
    data='';
    N=0;
    miu=3.986e5;

    while ~feof(fidin)        %判断是否文件读完
        tline=fgetl(fidin);  %逐行读入数据
        if isempty(tline)
            continue
        end


        if mod(N,3)==0
            tline = pad(tline, 24, 'right');   % 在右侧补空格到 24 个字符
        end        
        data=[data tline];
        N=N+1;  %统计文件的行数
    end

    satlitNum=N/3;      %卫星数目
    X=zeros(7,satlitNum);
    fclose(fidin);
end

% satlitNum=100;
%提取数据
t0=str2double(data(-120+162*1:-106+162*1));
t0_string=num2str(floor(t0));
year_start=str2double(t0_string(1:2));



for i=1:satlitNum
%     if str2double(data(-84+162*i:-79+162*i))==0
%         satlitNum=satlitNum-1;
%         continue;
%     end
    bianhao=str2double(data(-136+162*i:-131+162*i));
    X(1,i)=bianhao;        %将卫星的编号存到矩阵第一行中

    t=str2double(data(-120+162*i:-106+162*i));
    X(2,i)=t;   %将数据的时间存到矩阵第二行中

    t_string=num2str(floor(t));
    year_now=str2double(t_string(1:2));


    
    dt=t-t0;

    dt=dt-635*(year_now-year_start);


    if dt>dt_use
        break
    end

    


    B1=str2double(data(-84+162*i:-79+162*i));
    if(B1<0)
        B1=-B1;
    end
    B2=str2double(data(-78+162*i:-77+162*i));
    X(3,i)=B1*10^(-5)*10^B2;  %将弹道系数存到矩阵第三行

    e=str2double(data(-42+162*i:-36+162*i));
    e=e*10^(-7);
    X(4,i)=e;         %将偏心率保存到矩阵第四行中

    omiga=str2double(data(-34+162*i:-27+162*i));
    X(5,i)=omiga;          %将近地点辐角存入矩阵的第五行

    n=str2double(data(-16+162*i:-6+162*i));
    X(6,i)=n;          %将转速存入矩阵的第六行
    a=(86400^2*miu/(4*pi^2*n^2))^(1/3);
    X(7,i)=a;    %将半长轴存入矩阵的第七行
end
satlitNum=i-1;
if dt<dt_use
    disp(["error!",fileName]);
end

i=1;
while i<satlitNum
    if(abs(X(2,i+1)-X(2,i))<0.001)
        for j=i:(satlitNum-1)
            for m=1:7
                X(m,j)=X(m,j+1);
            end
        end
        satlitNum=satlitNum-1;
    else
        i=i+1;
    end
end

% %画散点图
% %偏心率与碎片序号的关系
% figure(1);
% x1=1:satlitNum;
% y1=X(4,x1);
% gscatter(x1,y1);
% xlabel('碎片序号');
% ylabel('偏心率');
% title('偏心率与时间序号关系图');

% %半长轴与碎片序号的关系
% figure(num);
% x2=1:satlitNum;
% y2=X(7,x2);
% gscatter(x2,y2);
% xlabel('碎片序号');
% ylabel('半长轴/km');
% title('半长轴与时间序号关系图');

% %采样时刻在轨道周期中所处相位的相对百分比
% figure(3);
% x3=1:satlitNum;
% y3=X(2,x3)-X(2,1);
% for i=1:satlitNum
%     if y3(i)>600
%         y3(i)=y3(i)-635;
%     end
% 
%     y3(i)=X(6,1)*y3(i);
% 
% end

% gscatter(x3,y3-floor(y3));
% xlabel('时间序号');
% ylabel('相位的相对位置');
% title('采样时刻在轨道周期中所处相位的相对百分比');

% y4=y3-floor(y3);
% edges = 0:0.1:1;
% [Nbox, ~, binIdx] = histcounts(y4, edges);
% maxCount = max(Nbox);
% idx_max_bins = find(Nbox == maxCount);   % 可能是多个
% chosen_bin = idx_max_bins(1);
% idx_in_bin = find(binIdx == chosen_bin);

idx_in_bin=1:satlitNum;
%储存数据
[~,n]=size(idx_in_bin); %将数据的个数记为n
for i=1:n
    t(i)=X(2,idx_in_bin(i));
    B(i)=X(3,idx_in_bin(i));
    e(i)=X(4,idx_in_bin(i));
    w(i)=X(5,idx_in_bin(i));
    N(i)=X(6,idx_in_bin(i));
    a(i)=X(7,idx_in_bin(i));
end



%计算两次数据间偏心率的变化Δe、时间的变化Δt、半长轴的变化Δa
%把两次测量的面质比的均值作为这个时间段的面质比（后续可进行改进，如多项式拟合后积分，测试是否会增加准确度）
%把两次测量的转速均值作为这个时间段的转速
for i=1:(n-1)
    e_delta(i)=e(i+1)-e(i);
    t_delta(i)=t(i+1)-t(i);
    %     disp(t_delta(i))
    if(t_delta(i))>600
        t_delta(i)=t_delta(i)-635;
    end
    t_delta(i)=t_delta(i)*86400;

    a_delta(i)=a(i+1)-a(i);
    B_mean(i)=(B(i)+B(i+1))/2;
    N_mean(i)=(N(i)+N(i+1))/2;

    N_mean(i)=N_mean(i)/86400;


end

%计算偏心率和半长轴在单位时间的变化（对于同一个碎片一段较短时间内，周期可以视为不变）
e_delta_t=abs(e_delta./t_delta);
a_delta_t=abs(a_delta./t_delta);



%计算一周期内的Δa
% a_delta_t_N=(a_delta_t')./(N_mean')*2*pi;
a_delta_t_N=a_delta_t';

x = B_mean(:); y = a_delta_t_N(:);
Zx = abs((x - mean(x,'omitnan')) ./ std(x,'omitnan'));
Zy = abs((y - mean(y,'omitnan')) ./ std(y,'omitnan'));
mask = (Zx <= 3) & (Zy <= 3);     % 3σ 规则，阈值可调成 2.5 或 2
% disp(length(mask))
% disp(sum(mask))
% mask=true(length(Zx),1);
B_mean = x(mask); a_delta_t_N = y(mask).';
all=length(mask);
use=sum(mask);

%计算拟合优度
p4=polyfit(B_mean,a_delta_t_N,1);
f4=polyval(p4,B_mean');
SSres4 = sum((f4 - a_delta_t_N).^2); % 残差平方和
SStot4 = sum((a_delta_t_N - mean(a_delta_t_N)).^2); % 总平方和
R_squared4 = 1 - (SSres4 / SStot4); % R-squared
% disp(R_squared4)  %计算delta_e的拟合优度
R2=R_squared4;
x_linear=linspace(0,1.05*max(B_mean),100);
%一周期Δa-A/M及回归直线作图
if plot_flag==1
    figure(num+1);
    scatter(B_mean.', a_delta_t_N, 8, 'green', 'filled')
    hold on;
    
    plot(x_linear,x_linear*p4(1)+p4(2),'Color','magenta','LineWidth',1);
    hold off;
    xlabel('B^*');
    ylabel('\Deltaa/\DeltaT');
    title({'\Deltaa单位时间变化率与面质比关系图';['R2= ',num2str(R_squared4)]});
    xlim([0 1.05*max(B_mean)]);  % 将X轴下限设为0，上限自动
    ylim([0 1.05*max(a_delta_t_N)]);  % 将Y轴下限设为0，上限自动
end

pointxy=[B_mean.';a_delta_t_N];
linexy=[x_linear;x_linear*p4(1)+p4(2)];

end
