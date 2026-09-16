function [JS_sandu_R1,JS_pingjun_R1,JS_zuida_R1]=zhuchengxu_R1(zhongzi,ns)
% Comment 5: matching-metric ablation
% 1: baseline (r + w)
% 2: baseline + i
% 3: baseline + Omega
% 4: baseline + i + Omega
%
% 本程序直接读取 zhuchengxu.m 已经保存的计算结果，不重新进行大、小碎片的参考轨道积分。
% 其中 baseline 直接由原有的 xiao_tuiyan/Px_tuiyan 重新计算 JS 散度；
% 其余三组仅从“大碎片轨道数据归并排序”之后重新进行匹配和小碎片推演。

data_dir='data';
filename=fullfile(data_dir,sprintf('zhongzi=%d_ns=%d.mat',zhongzi,ns));
if(~isfile(filename))
    error('Cannot find saved data file: %s',filename);
end


%% 先读取原有空间分布结果，重新计算 baseline 的 JS 散度
% 这里不直接使用原文件中保存的 JS_sandu，是为了避免旧服务器版本中
% JS 边界项计算代码可能存在差异。
load(filename,'Ts','xiao_cankao','Px_cankao','xiao_tuiyan','Px_tuiyan',...
    'lo_fenbu','hi_fenbu','dR_fenbu');

JS_sandu_R1=zeros(4,ns+1);
JS_sandu_R1(1,:)=jisuan_JS(xiao_cankao,Px_cankao,xiao_tuiyan,Px_tuiyan,...
    hi_fenbu,lo_fenbu,dR_fenbu,ns);

clear xiao_tuiyan Px_tuiyan;

fprintf('baseline finished: mean JS = %.8f, max JS = %.8f\n',...
    mean(JS_sandu_R1(1,:)),max(JS_sandu_R1(1,:)));


%% 读取从大碎片轨道数据归并排序之后所需要的变量
load(filename,'taewda','taewda_shu','yda_ode','s_m','s_mx','no_da','no_x',...
    'ts','miu_earth','rcx','vcx','r0','v0');


% 暂时不去考虑推演大碎片面质比的问题，认为大碎片面质比为已知量
% 如果说想改变采样周期的话，无需调整这一部分的代码，只需改变所使用的的数据即可，即改变ncaiyang
% 对大碎片每隔ncaiyang个周期采一次样
% 处理得出的大碎片的数据
% 由于从实际出发，小碎片推演时只能够使用当前时刻之前的数据，因此需要对大碎片的轨道数据进行归并排序
% 使用多次二路归并，多路归并效率远低于多次二路归并
% 但似乎影响不大，因此为了代码的简便这里暂时使用多路归并
tic;
ncaiyang=5;

% 相比原 aeda_sm 增加两行：
% 第8行：上一次采样时的大碎片轨道倾角 i
% 第9行：上一次采样时的大碎片升交点赤经 Omega
% 原始计算没有保存“加噪后的”i和Omega，因此这里直接由已经保存的大碎片
% 笛卡尔状态 yda_ode 提取，不重新积分参考轨道。
aeda_sm_R1=zeros(9,sum(ceil((taewda_shu-1)/ncaiyang))+no_da)-2;
k=1;%排列后的数据的列号

dangqian=zeros(no_da,1)+1+ncaiyang;
while(sum(dangqian>taewda_shu)~=no_da)
    tda_min=ts+1;%比所有的时刻值都大
    i_choose=1;
    for i=1:no_da
        if(dangqian(i)<=taewda_shu(i)&&taewda(1,dangqian(i),i)<tda_min)
            i_choose=i;
            tda_min=taewda(1,dangqian(i),i);
        end
    end

    aeda_sm_R1(1,k)=taewda(1,dangqian(i_choose),i_choose);%本次采样的时刻值
    aeda_sm_R1(2,k)=taewda(2,dangqian(i_choose)-ncaiyang,i_choose);%上一次采样时的a
    aeda_sm_R1(3,k)=taewda(3,dangqian(i_choose)-ncaiyang,i_choose);%上一次采样时的e
    aeda_sm_R1(4,k)=taewda(4,dangqian(i_choose)-ncaiyang,i_choose);%上一次采样时的w
    aeda_sm_R1(5,k)=(taewda(2,dangqian(i_choose),i_choose)-taewda(2,dangqian(i_choose)-ncaiyang,i_choose))...
        /(taewda(1,dangqian(i_choose),i_choose)-taewda(1,dangqian(i_choose)-ncaiyang,i_choose))/s_m(i_choose);%delta_a/delta_t/s_m
    aeda_sm_R1(6,k)=(taewda(3,dangqian(i_choose),i_choose)-taewda(3,dangqian(i_choose)-ncaiyang,i_choose))...
        /(taewda(1,dangqian(i_choose),i_choose)-taewda(1,dangqian(i_choose)-ncaiyang,i_choose))/s_m(i_choose);%delta_e/delta_t/s_m
    aeda_sm_R1(7,k)=s_m(i_choose);%大碎片面质比

    [aeda_sm_R1(8,k),aeda_sm_R1(9,k)]=r0v0_i_Omega(...
        yda_ode(dangqian(i_choose)-ncaiyang,1:3,i_choose).',...
        yda_ode(dangqian(i_choose)-ncaiyang,4:6,i_choose).');

    dangqian(i_choose)=dangqian(i_choose)+ncaiyang;
    k=k+1;
end

[~,aeda_sm_lie]=size(aeda_sm_R1);
for i=aeda_sm_lie:-1:1
    if(aeda_sm_R1(1,i)<0)
        aeda_sm_R1(:,i)=[];
    end
end

t_guibing=toc;
fprintf('large-debris record merging finished: %.2f s\n',t_guibing);

% 后续不再需要大碎片笛卡尔状态，释放内存
clear taewda taewda_shu yda_ode s_m dangqian;


%% 分别计算 +i、+Omega、+i+Omega 三组结果
% baseline 已经由原有推演结果获得，因此不重复计算。
mode_name={'baseline','+i','+Omega','+i+Omega'};
baoliu=2;%使用t1_jz~t2_jz时刻的大碎片轨道数据进行推演（与原程序一致）
t_da_caiyang=aeda_sm_R1(1,:);%采样时刻的升序序列
[~,t_da_caiyang_lie]=size(t_da_caiyang);
if(t_da_caiyang(end)<0)
    fprintf('error!caiyang_guibing');
    pause;
end

max_tuiyan_shu=ceil(ts/(2*pi*sqrt(6378000^3/miu_earth)))+10;

% 本消融实验中 tae_xiao_tuiyan_R1 为大体量三维数组。
% 如果JS阶段在process-based pool中按时间parfor，会把整个数组复制到各worker，
% 容易出现“反序列化期间内存不足”。因此优先使用thread-based pool。
% thread pool共享同一进程内存，可同时用于下面的小碎片推演和JS计算。
use_thread_pool_R1=prepare_parallel_pool_R1();

for mode=2:4
    tic;
    fprintf('\nstart mode %d: %s\n',mode,mode_name{mode});

    % 根据大碎片运动情况推演得出的小碎片的运动情况
    tae_xiao_tuiyan_R1=zeros(3,max_tuiyan_shu,no_x);
    k5_R1=zeros(no_x,1);

    parfor i=1:no_x
        % u(1:5)分别对应a、e、w、i、Omega
        u=zeros(5,1);%递推过程的中间变量
        v=zeros(5,1);%递推过程的中间变量
        [u(1),u(2),u(3)]=r0v0_genshu(rcx(:,i),vcx(:,i),miu_earth);
        [u(4),u(5)]=r0v0_i_Omega(rcx(:,i),vcx(:,i));

        t_xiao_tuiyan=0;
        k5_linshi=0;
        tae_xiao_tuiyan_linshi=zeros(3,max_tuiyan_shu);

        while 1
            k5_linshi=k5_linshi+1;
            tae_xiao_tuiyan_linshi(1,k5_linshi)=t_xiao_tuiyan;
            tae_xiao_tuiyan_linshi(2,k5_linshi)=u(1);
            tae_xiao_tuiyan_linshi(3,k5_linshi)=u(2);

            if(t_xiao_tuiyan==0)%由于第一个周期数据量较少，因此t_xiaotuiyan需要多加一个轨道周期的时间,以此确保可以访问到大碎片的轨道数据
                dt_tuiyan=(ncaiyang+1)*2*pi*sqrt(u(1)^3/miu_earth);
                t_xiao_tuiyan=t_xiao_tuiyan+dt_tuiyan;
                flag_1=1;
            else%否则的话，一个轨道周期一个轨道周期地进行递推
                dt_tuiyan=2*pi*sqrt(u(1)^3/miu_earth);
                t_xiao_tuiyan=t_xiao_tuiyan+dt_tuiyan;
                flag_1=0;
            end

            if(t_xiao_tuiyan>ts)
                break;
            end
            if(t_xiao_tuiyan<t_da_caiyang(1))
                fprintf('error!t_xiao_tuiyan');
                pause;
            end

            t1_jz=max(0,t_xiao_tuiyan-baoliu*ncaiyang*Ts);
            t2_jz=t_xiao_tuiyan;
            t2_jz_no=erfen_R1(t_da_caiyang,1,t_da_caiyang_lie,t2_jz);
            t1_jz_no=erfen_R1(t_da_caiyang,1,t2_jz_no,t1_jz);

            % 这里与原程序相同：首先按照近地点/远地点筛选2%的径向相似候选，
            % 然后在候选中进行第二阶段匹配。mode决定是否在w之外加入i和Omega。
            [lie_a_choose,lie_e_choose]=sousuo_R1(...
                aeda_sm_R1([2,3,4,8,9],t1_jz_no:t2_jz_no),...
                u(1),u(2),u(3),u(4),u(5),mode);%这里的lie_*_choose是部分矩阵的列号

            lie_a_choose=lie_a_choose+t1_jz_no-1;
            lie_e_choose=lie_e_choose+t1_jz_no-1;

            if(flag_1==1)
                v(1)=u(1)+aeda_sm_R1(5,lie_a_choose)*dt_tuiyan*s_mx(i);
                v(2)=u(2)+aeda_sm_R1(6,lie_e_choose)*dt_tuiyan*s_mx(i);
                v(3)=u(3)+3/4*1.083e-3*(6378e3/u(1)/(1-u(2))^2)^2*sqrt(miu_earth/u(1)^3)...
                    *(5*dot(cross(r0,v0)/norm(cross(r0,v0)),[0;0;1])^2-1)*dt_tuiyan;
            else
                v(1)=u(1)+aeda_sm_R1(5,lie_a_choose)*dt_tuiyan*s_mx(i);
                v(2)=u(2)+aeda_sm_R1(6,lie_e_choose)*dt_tuiyan*s_mx(i);
                v(3)=u(3)+3/4*1.083e-3*(6378e3/u(1)/(1-u(2))^2)^2*sqrt(miu_earth/u(1)^3)...
                    *(5*dot(cross(r0,v0)/norm(cross(r0,v0)),[0;0;1])^2-1)*dt_tuiyan;
            end

            if(v(2)<0)
                v(2)=0;
            else
                if(v(2)>1)
                    v(2)=1;
                end
            end

            % 一阶J2下倾角无长期漂移，因此i保持初始值。
            % Omega仅用于本次消融实验的匹配，不使用小碎片reference未来信息，
            % 而是根据一阶J2长期项自行递推。
            v(4)=u(4);
            p_R1=u(1)*(1-u(2)^2);
            if(p_R1<=0)
                p_R1=eps;
            end
            Omega_dot=-3/2*1.083e-3*sqrt(miu_earth/u(1)^3)*(6378e3/p_R1)^2*cos(u(4));
            v(5)=mod(u(5)+Omega_dot*dt_tuiyan,2*pi);

            u=v;
            if(u(1)*(1-u(2))<=6378e3)%坠落，之后的数据无需记载
                k5_linshi=k5_linshi+1;
                tae_xiao_tuiyan_linshi(1,k5_linshi)=t_xiao_tuiyan;
                tae_xiao_tuiyan_linshi(2,k5_linshi)=-1;
                tae_xiao_tuiyan_linshi(3,k5_linshi)=-1;
                break;
            end
        end

        k5_R1(i)=k5_linshi;
        % parfor切片变量要求除循环索引i外，其余维度使用固定索引或冒号。
        % 因此这里写回完整的固定尺寸临时数组；实际有效长度仍由k5_R1(i)记录。
        tae_xiao_tuiyan_R1(:,:,i)=tae_xiao_tuiyan_linshi;
    end

    if(max(k5_R1)>max_tuiyan_shu)
        fprintf('error!max_tuiyan_shu');
        pause;
    end

    t_tuiyan=toc;
    fprintf('%s inference finished: %.2f h\n',mode_name{mode},t_tuiyan/3600);


    % 计算该组推演结果与reference之间的JS散度
    tic;
    JS_sandu_R1(mode,:)=jisuan_JS_from_tae(tae_xiao_tuiyan_R1,k5_R1,...
        xiao_cankao,Px_cankao,lo_fenbu,hi_fenbu,dR_fenbu,no_x,ns,Ts,use_thread_pool_R1);
    t_JS=toc;

    fprintf('%s: mean JS = %.8f, max JS = %.8f, JS calculation = %.2f s\n',...
        mode_name{mode},mean(JS_sandu_R1(mode,:)),max(JS_sandu_R1(mode,:)),t_JS);

    clear tae_xiao_tuiyan_R1 k5_R1;
end


%% 汇总结果并绘制四种matching metric的JS散度变化图
JS_pingjun_R1=mean(JS_sandu_R1,2);
JS_zuida_R1=max(JS_sandu_R1,[],2);

fprintf('\n================ matching-metric ablation ================\n');
for mode=1:4
    fprintf('%-10s  mean JS = %.8f   max JS = %.8f\n',...
        mode_name{mode},JS_pingjun_R1(mode),JS_zuida_R1(mode));
end
fprintf('==========================================================\n');


t_days=(0:ns)*Ts/86400;
f_tu=figure;
plot(t_days,JS_sandu_R1(1,:),'LineWidth',1.2);
hold on;
plot(t_days,JS_sandu_R1(2,:),'LineWidth',1.2);
plot(t_days,JS_sandu_R1(3,:),'LineWidth',1.2);
plot(t_days,JS_sandu_R1(4,:),'LineWidth',1.2);
hold off;

xlabel('Time (days)');
ylabel('Jensen-Shannon divergence');
legend('Baseline','+i','+\Omega','+i+\Omega','Location','best');
grid on;
box on;
set(gca,'FontSize',12);
set(f_tu,'Color','w');

fig_name=sprintf('zhongzi=%d_ns=%d_matching_metric_ablation_JS',zhongzi,ns);
savefig(f_tu,fullfile(data_dir,[fig_name,'.fig']));
try
    exportgraphics(f_tu,fullfile(data_dir,[fig_name,'.png']),'Resolution',600);
catch
    print(f_tu,fullfile(data_dir,[fig_name,'.png']),'-dpng','-r600');
end


%% 只保存本次消融实验需要的结果，避免再次保存大体量reference数据
result_filename=fullfile(data_dir,sprintf('zhongzi=%d_ns=%d_R1.mat',zhongzi,ns));
save(result_filename,'JS_sandu_R1','JS_pingjun_R1','JS_zuida_R1','mode_name',...
    'zhongzi','ns','Ts','ncaiyang','baoliu','-v7.3');

fprintf('R1 results saved to %s\n',result_filename);
end



function [lie_a_choose,lie_e_choose]=sousuo_R1(aeda_sm,a,e,w,inc,Omega,mode)
% aeda_sm矩阵中，第1~5行依次为a、e、w、i、Omega
% 与原sousuo函数保持相同的两级匹配结构：
% 1) 先用近地点/远地点差筛选前2%%；
% 2) 再使用角度匹配选择最终记录。
% mode=2: w+i; mode=3: w+Omega; mode=4: w+i+Omega

[~,max_lie]=size(aeda_sm);
chushai=ceil(0.02*max_lie);

% 第一阶段：径向范围匹配
E_r=abs(aeda_sm(1,:).*(1+aeda_sm(2,:))-a*(1+e))+...
    abs(aeda_sm(1,:).*(1-aeda_sm(2,:))-a*(1-e));

% 使用升序排列得到与原程序相同含义的前2%候选
[~,xuhao]=sort(E_r,'ascend');
cs_xuhao=xuhao(1:chushai);

% 第二阶段：角度匹配
E_w=zeros(1,chushai);
E_i=zeros(1,chushai);
E_Omega=zeros(1,chushai);
for j=1:chushai
    E_w(j)=min(abs(aeda_sm(3,cs_xuhao(j))-w),2*pi-abs(aeda_sm(3,cs_xuhao(j))-w));
    E_i(j)=abs(aeda_sm(4,cs_xuhao(j))-inc);
    E_Omega(j)=min(abs(aeda_sm(5,cs_xuhao(j))-Omega),2*pi-abs(aeda_sm(5,cs_xuhao(j))-Omega));
end

if(mode==2)
    E_angle=E_w+E_i;
else
    if(mode==3)
        E_angle=E_w+E_Omega;
    else
        if(mode==4)
            E_angle=E_w+E_i+E_Omega;
        else
            error('Unsupported matching mode: %d',mode);
        end
    end
end

[~,lie_min]=min(E_angle);

% 保留原程序针对近圆轨道的处理：
% 当e较小时，半长轴变化只采用第一阶段中径向最接近的候选；
% 偏心率变化仍采用第二阶段角度匹配结果。
if(e>3e-3)
    lie_a_choose=cs_xuhao(lie_min);
    lie_e_choose=cs_xuhao(lie_min);
else
    lie_a_choose=cs_xuhao(1);
    lie_e_choose=cs_xuhao(lie_min);
end
end



function [inc,Omega]=r0v0_i_Omega(r0,v0)
% 由笛卡尔状态计算轨道倾角和升交点赤经
h=cross(r0,v0);
h_norm=norm(h);

cosi=h(3)/h_norm;
cosi=max(-1,min(1,cosi));
inc=acos(cosi);

% n=k×h=[-hy,hx,0]，因此Omega=atan2(hx,-hy)
if(abs(sin(inc))<1e-12)
    Omega=0;
else
    Omega=mod(atan2(h(1),-h(2)),2*pi);
end
end



function xuhao=erfen_R1(xulie,lo,hi,zhi)%返回序号对应的值<=查找值，序号+1对应的值大于查找值
if(hi<=lo)
    xuhao=lo;
    return;
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



function JS_sandu=jisuan_JS(xiao_cankao,Px_cankao,xiao_tuiyan,Px_tuiyan,hi_fenbu,lo_fenbu,dR_fenbu,ns)
% 根据已经计算好的reference和inferred空间分布重新计算JS散度
JS_sandu=zeros(1,ns+1);
no_bin=(hi_fenbu-lo_fenbu)/dR_fenbu;

for i=1:ns+1
    P_KL=xiao_cankao(1,:,i);
    Q_KL=xiao_tuiyan(1,:,i);

    for j=1:no_bin
        if(P_KL(j)~=0)
            JS_sandu(i)=JS_sandu(i)+0.5*P_KL(j)*log(P_KL(j)*2/(P_KL(j)+Q_KL(j)));
        end
        if(Q_KL(j)~=0)
            JS_sandu(i)=JS_sandu(i)+0.5*Q_KL(j)*log(Q_KL(j)*2/(P_KL(j)+Q_KL(j)));
        end
    end

    if(Px_cankao(i,1)~=0)
        JS_sandu(i)=JS_sandu(i)+0.5*Px_cankao(i,1)*log(Px_cankao(i,1)*2/(Px_cankao(i,1)+Px_tuiyan(i,1)));
    end
    if(Px_cankao(i,2)~=0)
        JS_sandu(i)=JS_sandu(i)+0.5*Px_cankao(i,2)*log(Px_cankao(i,2)*2/(Px_cankao(i,2)+Px_tuiyan(i,2)));
    end
    if(Px_tuiyan(i,1)~=0)
        JS_sandu(i)=JS_sandu(i)+0.5*Px_tuiyan(i,1)*log(Px_tuiyan(i,1)*2/(Px_cankao(i,1)+Px_tuiyan(i,1)));
    end
    if(Px_tuiyan(i,2)~=0)
        JS_sandu(i)=JS_sandu(i)+0.5*Px_tuiyan(i,2)*log(Px_tuiyan(i,2)*2/(Px_cankao(i,2)+Px_tuiyan(i,2)));
    end
end
end



function JS_sandu=jisuan_JS_from_tae(tae_xiao_tuiyan,k5,xiao_cankao,Px_cankao,...
    lo_fenbu,hi_fenbu,dR_fenbu,no_x,ns,Ts,use_thread_pool)
% 由新推演轨迹直接计算每个evaluation epoch的JS散度。
%
% 重要：这里不再对evaluation epoch直接做parfor并广播整个
% tae_xiao_tuiyan三维数组。即使使用thread-based pool，MATLAB在parfor
% 调度和函数调用过程中仍可能为大数组建立worker局部副本或临时副本，
% 从而在16个worker下造成较大的瞬时内存占用。
%
% 当前实现采用“分批提取 + 小数组并行”：
% 1) 主线程从tae_xiao_tuiyan中按时间批次提取a、e；
% 2) parfor仅接收当前批次的小型a/e矩阵和reference分布；
% 3) 每个worker独立计算一个epoch的空间分布和JS散度。
% 这样既保留并行计算，又避免将完整三维轨迹数组广播到worker。

JS_sandu=zeros(1,ns+1);

% 每批处理的evaluation epoch数量。
% 64个时刻对应的a/e矩阵通常只有几MB，远小于完整轨迹数组。
batch_size=64;
no_epoch=ns+1;
no_batch=ceil(no_epoch/batch_size);

for ibatch=1:no_batch
    i_begin=(ibatch-1)*batch_size+1;
    i_end=min(ibatch*batch_size,no_epoch);
    epoch_no=i_begin:i_end;
    nb=length(epoch_no);
    t_target=(epoch_no-1)*Ts;

    % 当前批次所有碎片在各evaluation epoch对应的a、e。
    % a<=0表示该碎片在该时刻已经再入。
    a_batch=zeros(no_x,nb);
    e_batch=zeros(no_x,nb);

    % 这一阶段保持串行，避免把完整tae_xiao_tuiyan送入parfor。
    % 对每个碎片只对当前批次第一个时刻做一次二分查找，后续时刻
    % 利用时间单调性向后移动索引，比每个时刻重新二分查找更省时间。
    for j=1:no_x
        kj=k5(j);
        if(kj<1)
            continue;
        end

        k_now=erfen_tae_R1(tae_xiao_tuiyan,j,kj,t_target(1));

        for b=1:nb
            tt=t_target(b);
            while(k_now<kj-1 && tae_xiao_tuiyan(1,k_now+1,j)<=tt)
                k_now=k_now+1;
            end

            a_now=tae_xiao_tuiyan(2,k_now,j);
            if(a_now>0)
                a_batch(j,b)=a_now;
                e_batch(j,b)=tae_xiao_tuiyan(3,k_now,j);
            end
        end
    end

    % 只截取当前批次所需的reference分布，避免在parfor中广播整个
    % xiao_cankao三维数组。
    P_ref_batch=permute(xiao_cankao(1,:,epoch_no),[3,2,1]);
    Px_ref_batch=Px_cankao(epoch_no,:);

    JS_batch=zeros(1,nb);

    if(use_thread_pool)
        parfor b=1:nb
            JS_batch(b)=jisuan_JS_one_epoch_from_ae_R1(...
                a_batch(:,b),e_batch(:,b),P_ref_batch(b,:),Px_ref_batch(b,:),...
                lo_fenbu,hi_fenbu,dR_fenbu,no_x);
        end
    else
        for b=1:nb
            JS_batch(b)=jisuan_JS_one_epoch_from_ae_R1(...
                a_batch(:,b),e_batch(:,b),P_ref_batch(b,:),Px_ref_batch(b,:),...
                lo_fenbu,hi_fenbu,dR_fenbu,no_x);
        end
    end

    JS_sandu(epoch_no)=JS_batch;

    % 输出进度，便于判断长时间计算是否正常推进。
    if(ibatch==1 || ibatch==no_batch || mod(ibatch,max(1,floor(no_batch/20)))==0)
        fprintf('    JS progress: %d/%d epochs (%.1f%%)\n',...
            i_end,no_epoch,100*i_end/no_epoch);
    end
end
end


function xuhao=erfen_tae_R1(tae_xiao_tuiyan,j,hi,zhi)
% 在第j个碎片的时间序列tae_xiao_tuiyan(1,:,j)中进行二分查找。
% 与erfen_R1含义相同，但不再构造tae_xiao_tuiyan(1,1:k5(j),j)的临时切片，
% 从而降低大规模重复调用时的内存分配和复制开销。

lo=1;
if(hi<=lo)
    xuhao=lo;
    return;
end

while 1
    mid=floor((lo+hi)/2);
    if(tae_xiao_tuiyan(1,mid,j)<=zhi)
        lo=mid;
    else
        hi=mid;
    end

    if(hi-lo<=1)
        xuhao=lo;
        return;
    end
end
end


function JS_i=jisuan_JS_one_epoch_from_ae_R1(a_all,e_all,P_KL,Px_ref,...
    lo_fenbu,hi_fenbu,dR_fenbu,no_x)
% 已知单个evaluation epoch全部碎片的a/e后，计算该时刻的JS散度。
% 本函数只接收当前epoch的小数组，不接收完整三维推演轨迹。

alive=(a_all>0);
xiao_tuiyan_shuliang=sum(alive);

% kongjianfenbu_R1原函数按照行向量读取，因此转置为行向量。
a_now=a_all(alive).';
e_now=e_all(alive).';

[Q_KL,~,P_shang,P_xia]=kongjianfenbu_R1(lo_fenbu,hi_fenbu,dR_fenbu,...
    a_now,e_now);

Q_KL=Q_KL/no_x;
Px_tuiyan_1=P_shang/no_x;
Px_tuiyan_2=(P_xia+no_x-xiao_tuiyan_shuliang)/no_x;

JS_i=0;
for jj=1:length(P_KL)
    if(P_KL(jj)~=0)
        JS_i=JS_i+0.5*P_KL(jj)*log(P_KL(jj)*2/(P_KL(jj)+Q_KL(jj)));
    end
    if(Q_KL(jj)~=0)
        JS_i=JS_i+0.5*Q_KL(jj)*log(Q_KL(jj)*2/(P_KL(jj)+Q_KL(jj)));
    end
end

if(Px_ref(1)~=0)
    JS_i=JS_i+0.5*Px_ref(1)*log(Px_ref(1)*2/(Px_ref(1)+Px_tuiyan_1));
end
if(Px_ref(2)~=0)
    JS_i=JS_i+0.5*Px_ref(2)*log(Px_ref(2)*2/(Px_ref(2)+Px_tuiyan_2));
end
if(Px_tuiyan_1~=0)
    JS_i=JS_i+0.5*Px_tuiyan_1*log(Px_tuiyan_1*2/(Px_ref(1)+Px_tuiyan_1));
end
if(Px_tuiyan_2~=0)
    JS_i=JS_i+0.5*Px_tuiyan_2*log(Px_tuiyan_2*2/(Px_ref(2)+Px_tuiyan_2));
end
end


function use_thread_pool=prepare_parallel_pool_R1()
% 优先建立thread-based pool，避免大数组在process worker之间复制。
% 若当前版本不支持thread pool，则建立普通process pool供小碎片推演使用，
% 同时令JS计算自动退回串行模式。

p=gcp('nocreate');
if(~isempty(p) && isa(p,'parallel.ThreadPool'))
    use_thread_pool=true;
    fprintf('Existing thread-based parallel pool will be used.\n');
    return;
end

if(~isempty(p))
    fprintf('Closing existing process-based parallel pool...\n');
    delete(p);
end

try
    parpool('Threads');
    use_thread_pool=true;
    fprintf('Thread-based parallel pool started. JS calculation will run in parallel without copying the full trajectory array.\n');
catch ME
    warning('Thread-based pool is unavailable: %s',ME.message);
    fprintf('Starting process-based pool for inference; JS calculation will use serial mode to avoid large-array broadcasting.\n');
    parpool('local');
    use_thread_pool=false;
end
end



function [yspace,xspace,P_shang,P_xia]=kongjianfenbu_R1(lo,hi,dR,a,e)
% 与原zhuchengxu.m中的kongjianfenbu保持一致
xspace=linspace(lo,hi-dR,(hi-lo)/dR)+dR/2;
yspace=zeros(1,(hi-lo)/dR);
[~,lie_a]=size(a);
[~,lie_e]=size(e);
if(lie_a~=lie_e)
    fprintf('error!kongjianfenbu');
    pause;
end
P_shang=0;
P_xia=0;

for i=1:lie_a
    % 与原程序保持一致。对于e极小但非零的情况直接使用原解析表达式。
    for j=1:(hi-lo)/dR
        r_down=lo+(j-1)*dR;
        r_up=lo+j*dR;
        if(a(i)*(1+e(i))>r_down&&a(i)*(1-e(i))<r_up)
            if(e(i)==0)
                if(a(i)>=r_down&&a(i)<r_up)
                    yspace(j)=yspace(j)+1;
                end
            else
                jifen_down=asin(max((r_down-a(i))/a(i)/e(i),-1));
                jifen_shang=asin(min((r_up-a(i))/a(i)/e(i),1));
                yspace(j)=yspace(j)+((jifen_shang-jifen_down)-e(i)*(cos(jifen_shang)-cos(jifen_down)))/pi;
            end
        end
    end

    if(a(i)*(1+e(i))>hi)
        if(e(i)==0)
            P_shang=P_shang+1;
        else
            jifen_down=asin(max((hi-a(i))/a(i)/e(i),-1));
            jifen_shang=asin(1);
            P_shang=P_shang+((jifen_shang-jifen_down)-e(i)*(cos(jifen_shang)-cos(jifen_down)))/pi;
        end
    end

    if(a(i)*(1-e(i))<lo)
        if(e(i)==0)
            P_xia=P_xia+1;
        else
            jifen_down=asin(-1);
            jifen_shang=asin(min((lo-a(i))/a(i)/e(i),1));
            P_xia=P_xia+((jifen_shang-jifen_down)-e(i)*(cos(jifen_shang)-cos(jifen_down)))/pi;
        end
    end
end

xspace=xspace/1e3;%单位变为1km
end
