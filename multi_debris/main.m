% 
% 
% clc;
clear;

%使用的基本参数
year0=2009;%推演的第1天的年份
day0=130;%数据第1天在该年的天数
t0=130;%推演的第1天的时间

miu_earth=398600e9;
h0=790e3;%碰撞高度约为790km
R=6378e3;%地球半径
lo_fenbu=6378e3+200e3;%考虑的最高轨道半径
hi_fenbu=6378e3+2*h0-200e3;%最低轨道半径
dR_fenbu=1e3;%每个区间的宽度
i0=74.035/180*pi; %轨道倾角，认为基本不变

DAY=365*10; %总的预测天数
%DAY=200;

ratio=0.5; %预测采用的大尺寸碎片的占比

% 指定读取的两行根数的文件夹路径
% folderPath = input('输入所在文件夹地址：', 's');
% folderPath = input('输入所在文件夹地址：', 's');
% folderPath = strtrim(folderPath); % 去除首尾空格
folderPath = 'cosmos_2251';
if ~isfolder(folderPath)
    error('The folder does not exist. Please check the path!');
end

% 获取文件夹中所有 .txt 文件的信息
filePattern = fullfile(folderPath, '*.txt');
% txtFiles = dir(filePattern);
txtFiles = dir(filePattern);
txtFiles = txtFiles(~startsWith({txtFiles.name}, '.')); % 排除隐藏文件
if isempty(txtFiles)
    error('The file does not exist. Please check the path!');
end

% 初始化一个存储数据的单元数组（或矩阵，视情况而定）
allData = {};
NUM=length(txtFiles);%总共的碎片个数，其中碰撞后本体的碎片数据不用

% 循环读取每个文件
h=waitbar(0,'Reading progress');
for k = 1:length(txtFiles)
    waitbar(k/NUM,h)
    % 获取文件的完整路径
    baseFileName = txtFiles(k).name;
    fullFileName = fullfile(folderPath, baseFileName);
    
    fidin = fopen(fullFileName, 'r');
    data='';
    N=0;
    miu=3.986e5; 

    while ~feof(fidin)        %判断是否文件读完
        tline=fgetl(fidin);  %逐行读入数据
        if isempty(tline)
            continue
        end
        data=[data tline];  
        N=N+1;  %统计文件的行数           
    end
    satlitNum=N/3;      %记录个数
    fclose(fidin);

    %提取数据
    X=zeros(7,satlitNum);
    for i=1:satlitNum 
        bianhao=str2double(data(-136+155*i:-131+155*i));
        X(1,i)=bianhao;        %将卫星的编号存到矩阵第一行中
    
        t=str2double(data(-120+155*i:-106+155*i));
        year1=2000+floor(t/1000);   %2000年之后的年份
        t1=t-floor(t/1000)*1000;
        for j=year0:(year1-1)
            t1=t1+365+isLeapYear(j);
        end
        X(2,i)=t1-t0;   %将数据相对第一天经过的时间存到矩阵第二行中
    
        B1=str2double(data(-84+155*i:-79+155*i));
%         if(B1<0)
%             B1=-B1;
%         end
        B2=str2double(data(-78+155*i:-77+155*i));
        X(3,i)=B1*10^(-5)*10^B2;  %将弹道系数存到矩阵第三行
    
        e=str2double(data(-42+155*i:-36+155*i));
        e=e*10^(-7);
        X(4,i)=e;         %将偏心率保存到矩阵第四行中
    
        omiga=str2double(data(-34+155*i:-27+155*i));
        X(5,i)=omiga/360*2*pi;          %将近地点辐角存入矩阵的第五行
    
        n=str2double(data(-16+155*i:-6+155*i)); 
        X(6,i)=n;          %将转速存入矩阵的第六行
        a=(86400^2*miu/(4*pi^2*n^2))^(1/3);
        X(7,i)=a;    %将半长轴存入矩阵的第七行
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
    
    for i=1:satlitNum
        M(i,5*k-4)=X(2,i);   %将数据的时间
        M(i,5*k-1)=X(3,i);  %将弹道系数
        M(i,5*k-2)=X(4,i);         %将偏心率
        M(i,5*k)=X(5,i);          %将近地点辐角
        M(i,5*k-3)=X(7,i);    %将半长轴存入
    end

end
delete(h);
disp('all file reading finished');

N=NUM; %总共的碎片个数
big=floor(N*ratio); %推演利用的大碎片的数量
small=N-big; %需要推演的小碎片的数量
fprintf('num of large：%d\n', big);
fprintf('num of small：%d\n', small);

real=zeros(5,DAY,N)-1; %每天采样一次，存储碎片的实际值(a/e/B)；第一个维度代表a\e\w\b\t四个数据，
%第二个维度代表DAY天每天采样时的数据，第三个维度代表N个碎片
%每天取一次,取当天的第一个数据，把数据存入real
for i=1:N
    for j=1:DAY
        k=1;
        while (M(k,5*i-4)-(j-1))<0
            k=k+1;
            if(M(k,5*i-4)==0)
                k=k-1;
                break
            end
        end
        real(1,j,i)=M(k,5*i-3)*1e3;%a，单位是m
        real(2,j,i)=M(k,5*i-2);%e
        real(3,j,i)=M(k,5*i);%w
        real(4,j,i)=M(k,5*i-1);%b
        real(5,j,i)=M(k,5*i-4);%t
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% load('data_lyk.mat')

%按照特征长度从大到小排序（特征长度近似与面质比成反比），存储在real数组
SM=zeros(1,N)-1;
for i=1:N
    SM(i)=sum(real(4,:,i));
end

middle=zeros(5,DAY)-1;%存储中间变量
for i=1:(N-1)
    SM1=SM(i);
    r=i;
    for j=(i+1):N
        if(SM(j)<SM1)
            SM1=SM(j);
            r=j;
        end
    end
    if(r==i)
        continue;
    else
        middle=real(:,:,i);
        real(:,:,i)=real(:,:,r);
        real(:,:,r)=middle;
        d=SM(i);
        SM(i)=SM(r);
        SM(r)=d;
    end
end

tic;
%根据前半部分碎片运动情况推演得出后半部分碎片运动情况
tae_latter_tuiyan=zeros(3,DAY,small)-1;%存储推演出的小碎片的tae
k5=zeros(small,1);

par=parpool('local',parcluster('local').NumWorkers, 'IdleTimeout', 60*24);
% parfor_progress(small); %计算进度
parfor i=1:small
    u=zeros(3,1);%递推过程的中间变量
    v=zeros(3,1);%递推过程的中间变量,存储下一天的数据
    u(1)=real(1,1,big+i);
    u(2)=real(2,1,big+i);
    u(3)=real(3,1,big+i);
    %后半碎片的初值，取出a/e/w存储在u
    t_latter_tuiyan=1;%推演第t_latter_tuiyan天的数据；直接从第一天开始推演
    tae_latter_tuiyan_linshi=zeros(3,DAY);%临时存储每天的tae_latter_tuiyan
    while 1
        tae_latter_tuiyan_linshi(1,t_latter_tuiyan)=t_latter_tuiyan;%存储第n天
        tae_latter_tuiyan_linshi(2,t_latter_tuiyan)=u(1);%存储a
        tae_latter_tuiyan_linshi(3,t_latter_tuiyan)=u(2);%存储e
        t_latter_tuiyan=t_latter_tuiyan+1;
        if(t_latter_tuiyan>DAY)
            break;
        end
        t1_jz=round(max(1,t_latter_tuiyan-10));
        t2_jz=t_latter_tuiyan-1;
        %t2_jz_no=erfen(t_da_caiyang,1,t_da_caiyang_lie,t2_jz);
        %t1_jz_no=erfen(t_da_caiyang,1,t2_jz_no,t1_jz);
        [day,num]=search(real(1:3,t1_jz:t2_jz,1:big),u(1),u(2),u(3));
        %day和num是查找到的最近的前半部分碎片的天数和序号
        %搜索：从第1天为止到现在的数据全搜一遍
        %v(1)=u(1)+aeda_sm(5,lie_a_choose)*2*pi*sqrt(u(1)^3/miu_earth)*s_mx(i);
        %v(2)=u(2)+aeda_sm(6,lie_e_choose)*2*pi*sqrt(u(1)^3/miu_earth)*s_mx(i);
        k=1;
%         while(real(5,day+k,num)==real(5,day,num))
%             k=k+1;
%         end

        v(1)=u(1)+(real(1,day+k,num)-real(1,day,num))/(real(5,day+k,num)-real(5,day,num)+1e-10)...
            /SM(num)*SM(big+i)...
            *(real(5,t_latter_tuiyan,big+i)-real(5,t_latter_tuiyan-1,big+i));
        %a
        v(2)=u(2)+(real(2,day+k,num)-real(2,day,num))/(real(5,day+k,num)-real(5,day,num)+1e-10)...
            /SM(num)*SM(big+i)...
            *(real(5,t_latter_tuiyan,big+i)-real(5,t_latter_tuiyan-1,big+i));
        %e
     
        %v(3)=u(3)+3/4*1.083e-3*(6378e3/u(1)/(1-u(2))^2)^2*sqrt(miu_earth/u(1)^3)*(5*dot(cross(r0,v0)/norm(cross(r0,v0)),[0;0;1])^2-1)*2*pi*sqrt(u(1)^3/miu_earth);
        v(3)=u(3)+3/4*1.083e-3*(6378e3/u(1)/(1-u(2))^2)^2*(5*cos(i0)^2-1)*sqrt(miu_earth/u(1)^3)*(real(5,t_latter_tuiyan,big+i)-real(5,t_latter_tuiyan-1,big+i));
        if(v(2)<0)
            v(2)=0;
        else
            if(v(2)>1)
                v(2)=1;
            end
        end
        u=v;
    end
   
    k5(i)=t_latter_tuiyan;
    for j=1:DAY
        if(j>=t_latter_tuiyan)
            break;
        end
        tae_latter_tuiyan(:,j,i)=tae_latter_tuiyan_linshi(:,j);
    end
% parfor_progress;
end
t4=toc;
% parfor_progress(0); 


% %将碎片按照特征长度分成4组
% 
% %绘制不同长度区间现实值和推演值
% tic;
% deltat=1;
% L1_fenbu_cankao=zeros(2,(hi_fenbu-lo_fenbu)/dR_fenbu,DAY);
% L2_fenbu_cankao=zeros(2,(hi_fenbu-lo_fenbu)/dR_fenbu,DAY);
% L3_fenbu_cankao=zeros(2,(hi_fenbu-lo_fenbu)/dR_fenbu,DAY);
% L4_fenbu_cankao=zeros(2,(hi_fenbu-lo_fenbu)/dR_fenbu,DAY);
% P4_shang_cankao=zeros(4,DAY);
% P4_xia_cankao=zeros(4,DAY);
% for i=1:DAY
%     [L1_fenbu_cankao(1,:,i),L1_fenbu_cankao(2,:,i),P4_shang_cankao(1,i),P4_xia_cankao(1,i)]=spaceDistribution(lo_fenbu,hi_fenbu,dR_fenbu,ae_cankao_L1(i,1:k1(i),1),ae_cankao_L1(i,1:k1(i),2));
%     L1_fenbu_cankao(1,:,i)=L1_fenbu_cankao(1,:,i)/no_da;
%     P4_shang_cankao(1,i)=P4_shang_cankao(1,i)/no_da;
%     P4_xia_cankao(1,i)=(P4_xia_cankao(1,i)+no_da-k1(i))/no_da;
%     [L2_fenbu_cankao(1,:,i),L2_fenbu_cankao(2,:,i),P4_shang_cankao(2,i),P4_xia_cankao(2,i)]=spaceDistribution(lo_fenbu,hi_fenbu,dR_fenbu,ae_cankao_L2(i,1:k2(i),1),ae_cankao_L2(i,1:k2(i),2));
%     L2_fenbu_cankao(1,:,i)=L2_fenbu_cankao(1,:,i)/Lcx_1_no_chu;
%     P4_shang_cankao(2,i)=P4_shang_cankao(2,i)/Lcx_1_no_chu;
%     P4_xia_cankao(2,i)=(P4_xia_cankao(2,i)+Lcx_1_no_chu-k2(i))/Lcx_1_no_chu;
%     [L3_fenbu_cankao(1,:,i),L3_fenbu_cankao(2,:,i),P4_shang_cankao(3,i),P4_xia_cankao(3,i)]=spaceDistribution(lo_fenbu,hi_fenbu,dR_fenbu,ae_cankao_L3(i,1:k3(i),1),ae_cankao_L3(i,1:k3(i),2));
%     L3_fenbu_cankao(1,:,i)=L3_fenbu_cankao(1,:,i)/(Lcx_2_no_chu-Lcx_1_no_chu);
%     P4_shang_cankao(3,i)=P4_shang_cankao(3,i)/(Lcx_2_no_chu-Lcx_1_no_chu);
%     P4_xia_cankao(3,i)=(P4_xia_cankao(3,i)+Lcx_2_no_chu-Lcx_1_no_chu-k3(i))/(Lcx_2_no_chu-Lcx_1_no_chu);
%     [L4_fenbu_cankao(1,:,i),L4_fenbu_cankao(2,:,i),P4_shang_cankao(4,i),P4_xia_cankao(4,i)]=spaceDistribution(lo_fenbu,hi_fenbu,dR_fenbu,ae_cankao_L4(i,1:k4(i),1),ae_cankao_L4(i,1:k4(i),2));
%     L4_fenbu_cankao(1,:,i)=L4_fenbu_cankao(1,:,i)/(no_x-Lcx_2_no_chu);
%     P4_shang_cankao(4,i)=P4_shang_cankao(4,i)/(no_x-Lcx_2_no_chu);
%     P4_xia_cankao(4,i)=(P4_xia_cankao(4,i)+no_x-Lcx_2_no_chu-k4(i))/(no_x-Lcx_2_no_chu);
% end
% flag_cankao=0;
% pic_num=1;
% while 1
%     if(pic_num>ns+1)
%         break;
%     end
%     f_tu=figure(1);
%     if(pic_num==1)
%         clf;
%     end
%     plot(L1_fenbu_cankao(1,:,pic_num),L1_fenbu_cankao(2,:,pic_num),L2_fenbu_cankao(1,:,pic_num),...
%         L2_fenbu_cankao(2,:,pic_num),L3_fenbu_cankao(1,:,pic_num),L3_fenbu_cankao(2,:,pic_num),L4_fenbu_cankao(1,:,pic_num),L4_fenbu_cankao(2,:,pic_num));
%     ylabel('轨道半径（km）');
%     xlabel('位于对应轨道半径区间的概率');
%     legend('碎片特征长度大于10cm的参考值',...
%         '碎片特征长度在3cm到10cm之间的参考值',...
%         '碎片特征长度在1.8cm到3cm之间的参考值',...
%         '碎片特征长度在1cm到1.8cm之间的参考值');
%     title({['演化时间=',num2str((pic_num-1)*Ts/86400),'天'];...
%         ['碎片特征长度大于10cm的参考值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(P4_shang_cankao(1,pic_num)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(P4_xia_cankao(1,pic_num))];...
%         ['碎片特征长度在3cm到10cm之间的参考值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(P4_shang_cankao(2,pic_num)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(P4_xia_cankao(2,pic_num))];...
%         ['碎片特征长度在1.8cm到3cm之间的参考值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(P4_shang_cankao(3,pic_num)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(P4_xia_cankao(3,pic_num))];...
%         ['碎片特征长度在1cm到1.8cm之间的参考值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(P4_shang_cankao(4,pic_num)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(P4_xia_cankao(4,pic_num))]});
%     set(gca,'FontSize',16);
%     if(flag_cankao==0)
%         axis([0 0.03 lo_fenbu/1e3 hi_fenbu/1e3]);
%     else
%         axis([0 0.01 lo_fenbu/1e3 hi_fenbu/1e3]);
%     end
%     set(f_tu, 'unit', 'normalized', 'position', [0,0,1,1]);
%     drawnow;
%     F_tu=getframe(gcf);
%     I_tu=frame2im(F_tu);
%     [I_tu,map]=rgb2ind(I_tu,256);
%     if(pic_num == 1)
%         imwrite(I_tu,map,['zhongzi=',num2str(zhongzi),',ns=',num2str(ns),',cankao_duo.gif'],'gif', 'Loopcount',inf,'DelayTime',0.05);
%         clf;
%     else
%         if(pic_num+huitu_jiange>ns+1)
%             if(flag_cankao==0)
%                 imwrite(I_tu,map,['zhongzi=',num2str(zhongzi),',ns=',num2str(ns),',cankao_duo.gif'],'gif','WriteMode','append','DelayTime',2);
%                 clf;
%             else
%                 imwrite(I_tu,map,['zhongzi=',num2str(zhongzi),',ns=',num2str(ns),',cankao_duo.gif'],'gif','WriteMode','append','DelayTime',5);
%             end
%         else
%             imwrite(I_tu,map,['zhongzi=',num2str(zhongzi),',ns=',num2str(ns),',cankao_duo.gif'],'gif','WriteMode','append','DelayTime',0.05);
%             clf;
%         end
%     end
%     if(flag_cankao==0&&pic_num+huitu_jiange>ns+1)
%         flag_cankao=1;
%         pic_num=pic_num-huitu_jiange;%特殊情况，不是错误，是因为需要重复绘制最后一幅图
%     end
%     pic_num=pic_num+huitu_jiange;
% end
% t5=toc;


%验证：先采用前一半碎片数据验证后一半的方法
%绘制后一半碎片参考值和真实值随时间变化的对比图
%计算后半碎片实际的空间分布情况（根据采样的实际数据）,存储在latter_real(概率的形式）
latter_real=zeros(2,(hi_fenbu-lo_fenbu)/dR_fenbu,DAY); %每1km一个高度区间；每个碎片第一个数据是高度区间，第二个数据是在该区间的概率
Px_real=zeros(DAY,2); %存储每个碎片在高度限制之外（hi_fenbu之上或lo_fenbu之下）的概率（一般是0）
for i=1:DAY
    [latter_real(1,:,i),latter_real(2,:,i),Px_real(i,1),Px_real(i,2)]=...
        spaceDistribution(lo_fenbu,hi_fenbu,dR_fenbu,...
        real(1,i,big+1:N),real(2,i,big+1:N));
    latter_real(1,:,i)=latter_real(1,:,i)/small;
    Px_real(i,1)=Px_real(i,1)/small;
    Px_real(i,2)=Px_real(i,2)/small;
end

%计算推演的空间分布
latter_tuiyan=zeros(2,(hi_fenbu-lo_fenbu)/dR_fenbu,DAY);%存储推演得到的后半碎片空间分布情况
Px_tuiyan=zeros(DAY,2);
latter_tuiyan_shuliang=zeros(N,1);
latter_tuiyan_a=zeros(1,small);
latter_tuiyan_e=zeros(1,small);
for i=1:DAY
    [latter_tuiyan(1,:,i),latter_tuiyan(2,:,i),Px_tuiyan(i,1),Px_tuiyan(i,2)]=...
        spaceDistribution(lo_fenbu,hi_fenbu,dR_fenbu,...
        tae_latter_tuiyan(2,i,1:small),tae_latter_tuiyan(3,i,1:small));
    latter_tuiyan(1,:,i)=latter_tuiyan(1,:,i)/small;
    Px_tuiyan(i,1)=Px_tuiyan(i,1)/small;
    Px_tuiyan(i,2)=Px_tuiyan(i,2)/small;
end

KL_sanduPQ=zeros(1,DAY);
KL_sanduQP=zeros(1,DAY);
JS_sandu=zeros(1,DAY);

%计算JS散度
for i=1:DAY
    P_KL=latter_real(1,:,i);
    Q_KL=latter_tuiyan(1,:,i);
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
    if(Px_real(i,1)~=0)
        KL_sanduPQ(i)=KL_sanduPQ(i)+Px_real(i,1)*log(Px_real(i,1)/Px_tuiyan(i,1));
        JS_sandu(i)=JS_sandu(i)+0.5*Px_real(i,1)*log(Px_real(i,1)*2/(Px_real(i,1)+Px_tuiyan(i,1)));
    end
    if(Px_real(i,2)~=0)
        KL_sanduPQ(i)=KL_sanduPQ(i)+Px_real(i,2)*log(Px_real(i,2)/Px_tuiyan(i,2));
        JS_sandu(i)=JS_sandu(i)+0.5*Px_real(i,2)*log(Px_real(i,2)*2/(Px_real(i,2)+Px_tuiyan(i,2)));
    end
    if(Px_tuiyan(i,1)~=0)
        KL_sanduQP(i)=KL_sanduQP(i)+Px_tuiyan(i,1)*log(Px_tuiyan(i,1)/Px_real(i,1));
        JS_sandu(i)=JS_sandu(i)+0.5*Px_tuiyan(i,1)*log(Px_tuiyan(i,1)*2/(Px_real(i,1)+Px_tuiyan(i,1)));
    end
    if(Px_tuiyan(i,2)~=0)
        KL_sanduQP(i)=KL_sanduQP(i)+Px_tuiyan(i,2)*log(Px_tuiyan(i,2)/Px_real(i,2));
        JS_sandu(i)=JS_sandu(i)+0.5*Px_tuiyan(i,2)*log(Px_tuiyan(i,2)*2/(Px_real(i,2)+Px_tuiyan(i,2)));
    end
end
% 
% %第1天的JS散度
% figure;
% pic_num=1;
% plot(latter_real(1,:,pic_num),latter_real(2,:,pic_num),'blue',latter_tuiyan(1,:,pic_num),latter_tuiyan(2,:,pic_num),'magenta');
% ylabel('轨道半径（km）');
% xlabel('位于对应轨道半径区间的概率');
% legend('小碎片的真实值','小碎片的推演值');
% title({['演化时间=',num2str(pic_num),'天'];['JS散度=',num2str(JS_sandu(pic_num))];...
%     ['真实值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(Px_real(pic_num,1)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(Px_real(pic_num,2))];...
%     ['推演值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(Px_tuiyan(pic_num,1)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(Px_tuiyan(pic_num,2))]});
% axis([0 0.02 lo_fenbu/1e3 hi_fenbu/1e3]); 
% %横纵坐标范围
% 
% %中间时刻的JS散度
% figure;
% pic_num=round(DAY/2);
% plot(latter_real(1,:,pic_num),latter_real(2,:,pic_num),'blue',latter_tuiyan(1,:,pic_num),latter_tuiyan(2,:,pic_num),'magenta');
% ylabel('轨道半径（km）');
% xlabel('位于对应轨道半径区间的概率');
% legend('小碎片的真实值','小碎片的推演值');
% title({['演化时间=',num2str(pic_num),'天'];['JS散度=',num2str(JS_sandu(pic_num))];...
%     ['真实值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(Px_real(pic_num,1)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(Px_real(pic_num,2))];...
%     ['推演值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(Px_tuiyan(pic_num,1)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(Px_tuiyan(pic_num,2))]});
% axis([0 0.02 lo_fenbu/1e3 hi_fenbu/1e3]);
% 
% %最终时刻的JS散度
% figure;
% pic_num=DAY;
% plot(latter_real(1,:,pic_num),latter_real(2,:,pic_num),'blue',latter_tuiyan(1,:,pic_num),latter_tuiyan(2,:,pic_num),'magenta');
% ylabel('轨道半径（km）');
% xlabel('位于对应轨道半径区间的概率');
% legend('小碎片的真实值','小碎片的推演值');
% title({['演化时间=',num2str(pic_num),'天'];['JS散度=',num2str(JS_sandu(pic_num))];...
%     ['真实值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(Px_real(pic_num,1)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(Px_real(pic_num,2))];...
%     ['推演值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(Px_tuiyan(pic_num,1)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(Px_tuiyan(pic_num,2))]});
% axis([0 0.02 lo_fenbu/1e3 hi_fenbu/1e3]);
% 
% %整个过程的JS散度变化
% figure;
% plot(JS_sandu);
% xlabel('Prediction Time[day]');
% ylabel('JS divergence');

vars = whos;  % 获取所有变量信息
figVars = {vars(strcmp({vars.class}, 'matlab.ui.Figure')).name};  % 找出figure变量名

% 要保存的变量列表（排除figure）
saveVars = setdiff({vars.name}, figVars);

% 构造文件名，注意加上扩展名.mat
filename = sprintf('true_10_year.mat');

% 保存变量，使用 -v7.3 支持大变量
save(filename, saveVars{:}, '-v7.3');

delete(par);

plot_SC;
