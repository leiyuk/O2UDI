function [L1_fenbu_cankao,L2_fenbu_cankao,L3_fenbu_cankao,L4_fenbu_cankao,P4_shang_cankao,P4_xia_cankao,xiao_cankao,xiao_tuiyan,Px_cankao,Px_tuiyan]...
    =zhuchengxu(zhongzi,ns,h0)
rng(zhongzi);%使用了随机数种子


%生成大碎片（特征长度大于10cm）的特征长度
s=1;
i=1;
Lc=zeros(1,1);
while 1
    Lc(i)=(i/6/s)^(-1/1.6);
    if(Lc(i)<0.1)
        no_da=i-1;
        break;
    end
    i=i+1;
end
Lc(no_da+1)=[];


%根据概率密度函数得出各个大碎片的面质比（原则上这一概率密度函数适用于特征长度大于11cm的情况，但这里近似认为大于10cm都服从这一分布）
nc=log10(Lc);
s_m=zeros(1,no_da);
opt_solve=optimset('Display','off');
parfor i=1:no_da
    r=rand;
    s_m(i)=fsolve(@(x1) integral(@(x) alfa(nc(i))*logn(miu1(nc(i)),sigma1(nc(i)),x),0,x1)+integral(@(x) (1-alfa(nc(i)))*logn(miu2(nc(i)),sigma2(nc(i)),x),0,x1)-r,0.12,opt_solve);
    while 1
        if(imag(s_m(i))==0&&s_m(i)>0&&integral(@(x) alfa(nc(i))*logn(miu1(nc(i)),sigma1(nc(i)),x),0,s_m(i))+integral(@(x) (1-alfa(nc(i)))*logn(miu2(nc(i)),sigma2(nc(i)),x),0,s_m(i))>0.005&&integral(@(x) alfa(nc(i))*logn(miu1(nc(i)),sigma1(nc(i)),x),0,s_m(i))+integral(@(x) (1-alfa(nc(i)))*logn(miu2(nc(i)),sigma2(nc(i)),x),0,s_m(i))<0.995)
            break;
        end
        r=rand;
        s_m(i)=fsolve(@(x1) integral(@(x) alfa(nc(i))*logn(miu1(nc(i)),sigma1(nc(i)),x),0,x1)+integral(@(x) (1-alfa(nc(i)))*logn(miu2(nc(i)),sigma2(nc(i)),x),0,x1)-r,0.12,opt_solve);
    end
end


%航天器解体后大碎片的速度分布（原则上这一概率密度函数适用于特征长度大于11cm的情况，但这里近似认为大于10cm都服从这一分布）
miu_earth=398600e9;
r0=[0;(6378e3+h0)*cos(pi/6);(6378e3+h0)*sin(pi/6)];
v0=[-sqrt(miu_earth/(6378e3+h0));0;0];
dv=zeros(1,no_da);
parfor i=1:no_da
    r=rand;
    dv(i)=fsolve(@(x1) integral(@(x) alfa(nc(i))*logn(0.2*miu1(nc(i))+1.85,sqrt(0.04*sigma1(nc(i))^2+0.16),x),0,x1)+integral(@(x) (1-alfa(nc(i)))*logn(0.2*miu2(nc(i))+1.85,sqrt(0.04*sigma2(nc(i))^2+0.16),x),0,x1)-r,40,opt_solve);
    while 1
        if(imag(dv(i))==0&&dv(i)>=0&&integral(@(x) alfa(nc(i))*logn(0.2*miu1(nc(i))+1.85,sqrt(0.04*sigma1(nc(i))^2+0.16),x),0,dv(i))+integral(@(x) (1-alfa(nc(i)))*logn(0.2*miu2(nc(i))+1.85,sqrt(0.04*sigma2(nc(i))^2+0.16),x),0,dv(i))>0.005&&integral(@(x) alfa(nc(i))*logn(0.2*miu1(nc(i))+1.85,sqrt(0.04*sigma1(nc(i))^2+0.16),x),0,dv(i))+integral(@(x) (1-alfa(nc(i)))*logn(0.2*miu2(nc(i))+1.85,sqrt(0.04*sigma2(nc(i))^2+0.16),x),0,dv(i))<0.995)
            break;
        end
        r=rand;
        dv(i)=fsolve(@(x1) integral(@(x) alfa(nc(i))*logn(0.2*miu1(nc(i))+1.85,sqrt(0.04*sigma1(nc(i))^2+0.16),x),0,x1)+integral(@(x) (1-alfa(nc(i)))*logn(0.2*miu2(nc(i))+1.85,sqrt(0.04*sigma2(nc(i))^2+0.16),x),0,x1)-r,40,opt_solve);
    end
end
rc=zeros(3,no_da);
vc=zeros(3,no_da);
parfor i=1:no_da
    sita=rand*pi;
    fai=rand*2*pi;
    rc(:,i)=r0;
    vc(:,i)=v0+[sin(sita)*cos(fai);sin(sita)*sin(fai);cos(sita)]*dv(i);
end
for i=1:no_da
    if(~(imag(s_m(i))==0&&s_m(i)>0&&imag(dv(i))==0&&dv(i)>=0))
        fprintf("error!s_m/dv");
        pause;
    end
end






%生成小碎片（特征长度在1~10cm之间）的特征长度
i=no_da+1;
j=1;%小碎片编号从1开始
Lcx=zeros(1,1);
while 1
    Lcx(j)=(i/6/s)^(-1/1.6);
    if(Lcx(j)<0.01)
        no_x=j-1;
        break;
    end
    i=i+1;
    j=j+1;
end
Lcx(no_x+1)=[];


%根据概率密度函数得出各小碎片的面质比（原则上这一概率密度函数适用于特征长度小于8cm的情况，但这里近似认为小于10cm都服从这一分布）
ncx=log10(Lcx);
s_mx=zeros(1,no_x);
parfor i=1:no_x
    r=rand;
    s_mx(i)=fsolve(@(x1) integral(@(x) logn(miu3(ncx(i)),sigma3(ncx(i)),x),0,x1)-r,0.4,opt_solve);
    while 1
        if(imag(s_mx(i))==0&&s_mx(i)>0&&integral(@(x) logn(miu3(ncx(i)),sigma3(ncx(i)),x),0,s_mx(i))>0.005&&integral(@(x) logn(miu3(ncx(i)),sigma3(ncx(i)),x),0,s_mx(i))<0.995)
            break;
        end
        r=rand;
        s_mx(i)=fsolve(@(x1) integral(@(x) logn(miu3(ncx(i)),sigma3(ncx(i)),x),0,x1)-r,0.4,opt_solve);
    end
end



%航天器解体后小碎片的速度分布（原则上这一概率密度函数适用于特征长度小于8cm的情况，但这里近似认为小于10cm都服从这一分布）
dvx=zeros(1,no_x);
parfor i=1:no_x
    r=rand;
    dvx(i)=fsolve(@(x1) integral(@(x) logn(0.2*miu3(ncx(i))+1.85,sqrt(0.04*sigma3(ncx(i))^2+0.16),x),0,x1)-r,60,opt_solve);
    while 1
        if(imag(dvx(i))==0&&dvx(i)>=0&&integral(@(x) logn(0.2*miu3(ncx(i))+1.85,sqrt(0.04*sigma3(ncx(i))^2+0.16),x),0,dvx(i))>0.005&&integral(@(x) logn(0.2*miu3(ncx(i))+1.85,sqrt(0.04*sigma3(ncx(i))^2+0.16),x),0,dvx(i))<0.995)
            break;
        end
        r=rand;
        dvx(i)=fsolve(@(x1) integral(@(x) logn(0.2*miu3(ncx(i))+1.85,sqrt(0.04*sigma3(ncx(i))^2+0.16),x),0,x1)-r,60,opt_solve);
    end
end
rcx=zeros(3,no_x);
vcx=zeros(3,no_x);
parfor i=1:no_x
    sita=rand*pi;
    fai=rand*2*pi;
    rcx(:,i)=r0;
    vcx(:,i)=v0+[sin(sita)*cos(fai);sin(sita)*sin(fai);cos(sita)]*dvx(i);
end
for i=1:no_x
    if(~(imag(s_mx(i))==0&&s_mx(i)>0&&imag(dvx(i))==0&&dvx(i)>=0))
        fprintf("error!s_mx/dvx");
        pause;
    end
end




Ts=2*pi*sqrt((6378e3+h0)^3/miu_earth);
ts=ns*Ts;%总演化时间
%大碎片作为一个区间，小碎片分为三个长度区间
Lcx_1=0.03;
Lcx_2=0.018;
for i=1:no_x
    if(Lcx(i)<Lcx_1)
        Lcx_1_no_chu=i;
        break;
    end
end
for i=1:no_x
    if(Lcx(i)<Lcx_2)
        Lcx_2_no_chu=i;
        break;
    end
end
Lcx_1_no_chu=Lcx_1_no_chu-1;
Lcx_2_no_chu=Lcx_2_no_chu-1;




disp('generated finished')

tic;
%小碎片直接积分演化（用于检验，从应用的角度来讲不能用这种方法计算）
ae_cankao_L2=zeros(ns+1,Lcx_1_no_chu,2)-1;
k2=zeros(ns+1,1);%对应的是第二个维度的实际长度，因为考虑到碎片存在坠落判断，因此不是每个地方都可以排满
ae_cankao_L3=zeros(ns+1,Lcx_2_no_chu-Lcx_1_no_chu,2)-1;
k3=zeros(ns+1,1);
ae_cankao_L4=zeros(ns+1,no_x-Lcx_2_no_chu,2)-1;
k4=zeros(ns+1,1);
yxiao_cankao=zeros(ns+1,6,no_x)-1;
hang_yxiao_cankao=zeros(no_x,1);
nxh=ceil(1.3*(ns+1));
parfor i=1:no_x
    opt=odeset('Events',@myEventsFcn_xiao,'RelTol', 1e-10, 'AbsTol', 1e-10);
    [~,yxiao_cankao_linshi,~,~,~]=ode45(@(t,xx) df(t,xx,s_mx(i)),linspace(0,ts,ns+1),[rcx(:,i);vcx(:,i)],opt);%之所以设置三个采样点只是为了格式要求
    [hang_yxiao_cankao_linshi,~]=size(yxiao_cankao_linshi);
    if(hang_yxiao_cankao_linshi~=ns+1)
        hang_yxiao_cankao_linshi=hang_yxiao_cankao_linshi-1;
    end
    hang_yxiao_cankao(i)=hang_yxiao_cankao_linshi;
    for j=1:nxh
        if(j>hang_yxiao_cankao_linshi)
            break;
        end
        yxiao_cankao(j,:,i)=yxiao_cankao_linshi(j,:);
    end
    if mod(i,100)==0
        disp(['index small debris = ',num2str(i)])
    end
end
for j=1:ns+1
    for i=1:no_x
        if(j<=hang_yxiao_cankao(i))
            if(i<=Lcx_1_no_chu)
                k2(j)=k2(j)+1;
                [ae_cankao_L2(j,k2(j),1),ae_cankao_L2(j,k2(j),2),~]=r0v0_genshu(yxiao_cankao(j,1:3,i).',yxiao_cankao(j,4:6,i).',miu_earth);
            else
                if(i<=Lcx_2_no_chu)
                    k3(j)=k3(j)+1;
                    [ae_cankao_L3(j,k3(j),1),ae_cankao_L3(j,k3(j),2),~]=r0v0_genshu(yxiao_cankao(j,1:3,i).',yxiao_cankao(j,4:6,i).',miu_earth);
                else
                    k4(j)=k4(j)+1;
                    [ae_cankao_L4(j,k4(j),1),ae_cankao_L4(j,k4(j),2),~]=r0v0_genshu(yxiao_cankao(j,1:3,i).',yxiao_cankao(j,4:6,i).',miu_earth);
                end
            end
        end
    end
end
if(sum(hang_yxiao_cankao>=nxh)~=0)
    fprintf('error!nxh');
    pause;
end
t1=toc;

disp('small debris finished')



tic;
%大碎片此后的演化状态
ae_cankao_L1=zeros(ns+1,no_da,2)-1;
k1=zeros(ns+1,1);
taewda=zeros(4,ceil(ts/(2*pi*sqrt(6378000^3/miu_earth))),no_da)-2;%第一个维度的1,2,3,4分别对应采样时间、半长轴、偏心率、近地点角
taewda_shu=zeros(no_da,1);%对应第二个维度上的实际有效数据
yda_cankao=zeros(ns+1,6,no_da)-1;
hang_yda_cankao=zeros(no_da,1);
yda_ode=zeros(ns+1,6,no_da)-1;
tda_ode=zeros(ns+1,no_da)-1;
parfor i=1:no_da
    opt=odeset('Events',@myEventsFcn_da,'RelTol', 1e-10, 'AbsTol', 1e-10);
    [~,yda_cankao_linshi,tda_ode_linshi,yda_ode_linshi,itda_ode_linshi]=ode45(@(t,xx) df(t,xx,s_m(i)),linspace(0,ts,ns+1),[rc(:,i);vc(:,i)],opt);%之所以设置三个采样点只是为了格式要求
    %处理参考值
    [hang_yda_cankao_linshi,~]=size(yda_cankao_linshi);
    if(hang_yda_cankao_linshi~=ns+1)
        hang_yda_cankao_linshi=hang_yda_cankao_linshi-1;
    end
    hang_yda_cankao(i)=hang_yda_cankao_linshi;
    for j=1:nxh
        if(j>hang_yda_cankao_linshi)
            break;
        end
        yda_cankao(j,:,i)=yda_cankao_linshi(j,:);
    end
    %处理采样值
    [hang_itda_ode_linshi,~]=size(itda_ode_linshi);
    if(itda_ode_linshi(end)==1)
        hang_itda_ode_linshi=hang_itda_ode_linshi-1;
    end
    taewda_shu(i)=hang_itda_ode_linshi;
    for j=1:nxh
        if(j>hang_itda_ode_linshi)
            break;
        end
        yda_ode(j,:,i)=yda_ode_linshi(j,:);
        tda_ode(j,i)=tda_ode_linshi(j);
    end
end
for j=1:ns+1
    for i=1:no_da
        if(j<=hang_yda_cankao(i))
            k1(j)=k1(j)+1;
            [ae_cankao_L1(j,k1(j),1),ae_cankao_L1(j,k1(j),2),~]=r0v0_genshu(yda_cankao(j,1:3,i).',yda_cankao(j,4:6,i).',miu_earth);
        end
    end
end
for i=1:no_da
    for j=1:taewda_shu(i)
        taewda(1,j,i)=tda_ode(j,i);
        [taewda(2,j,i),taewda(3,j,i),taewda(4,j,i)]=r0v0_genshu(yda_ode(j,1:3,i).'+zao_sheng(10),yda_ode(j,4:6,i).'+zao_sheng(0.2),miu_earth);%如有噪声在这里加到位置和速度上即可
    end
end
if(sum(taewda_shu>=nxh)~=0||sum(hang_yda_cankao>=nxh)~=0)
    fprintf('error!nxh');
    pause;
end
t2=toc;

disp('large debris finished')



%暂时不去考虑推演大碎片面质比的问题，认为大碎片面质比为已知量
%如果说想改变采样周期的话，无需调整这一部分的代码，只需改变所使用的的数据即可，即改变ncaiyang
%对大碎片每隔ncaiyang个周期采一次样
%处理得出的大碎片的数据
%由于从实际出发，小碎片推演时只能够使用当前时刻之前的数据，因此需要对大碎片的轨道数据进行归并排序
%使用多次二路归并，多路归并效率远低于多次二路归并
%但似乎影响不大，因此为了代码的简便这里暂时使用多路归并
tic;
ncaiyang=5;
aeda_sm=zeros(7,sum(ceil((taewda_shu-1)/ncaiyang))+no_da)-2;
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
    aeda_sm(1,k)=taewda(1,dangqian(i_choose),i_choose);%本次采样的时刻值
    aeda_sm(2,k)=taewda(2,dangqian(i_choose)-ncaiyang,i_choose);%上一次采样时的a
    aeda_sm(3,k)=taewda(3,dangqian(i_choose)-ncaiyang,i_choose);%上一次采样时的e
    aeda_sm(4,k)=taewda(4,dangqian(i_choose)-ncaiyang,i_choose);%上一次采样时的w
    aeda_sm(5,k)=(taewda(2,dangqian(i_choose),i_choose)-taewda(2,dangqian(i_choose)-ncaiyang,i_choose))...
        /(taewda(1,dangqian(i_choose),i_choose)-taewda(1,dangqian(i_choose)-ncaiyang,i_choose))/s_m(i_choose);            %delta_a/delta_t/s_m
    aeda_sm(6,k)=(taewda(3,dangqian(i_choose),i_choose)-taewda(3,dangqian(i_choose)-ncaiyang,i_choose))...
        /(taewda(1,dangqian(i_choose),i_choose)-taewda(1,dangqian(i_choose)-ncaiyang,i_choose))/s_m(i_choose);            %delta_e/delta_t/s_m
    aeda_sm(7,k)=s_m(i_choose);%大碎片面质比
    dangqian(i_choose)=dangqian(i_choose)+ncaiyang;
    k=k+1;
end
[~,aeda_sm_lie]=size(aeda_sm);
for i=aeda_sm_lie:-1:1
    if(aeda_sm(1,i)<0)
        aeda_sm(:,i)=[];
    end
end
t3=toc;




tic;
%根据大碎片运动情况推演得出的小碎片的运动情况
tae_xiao_tuiyan=zeros(3,ceil(ts/(2*pi*sqrt(6378000^3/miu_earth)))+10,no_x);
k5=zeros(no_x,1);
t_da_caiyang=aeda_sm(1,:);%采样时刻的升序序列
[~,t_da_caiyang_lie]=size(t_da_caiyang);
if(t_da_caiyang(end)<0)
    fprintf("error!caiyang_guibing");
    pause;
end
baoliu=2;%使用t1_jz~t2_jz时刻的大碎片轨道数据进行推演（这里的'2'这个值可以进行调整）
parfor i=1:no_x
    u=zeros(3,1);%递推过程的中间变量
    v=zeros(3,1);%递推过程的中间变量
    [u(1),u(2),u(3)]=r0v0_genshu(rcx(:,i),vcx(:,i),miu_earth);
    t_xiao_tuiyan=0;
    k5_linshi=0;
    tae_xiao_tuiyan_linshi=zeros(3,ceil(ts/(2*pi*sqrt(6378000^3/miu_earth)))+10);
    while 1
        k5_linshi=k5_linshi+1;
        tae_xiao_tuiyan_linshi(1,k5_linshi)=t_xiao_tuiyan;
        tae_xiao_tuiyan_linshi(2,k5_linshi)=u(1);
        tae_xiao_tuiyan_linshi(3,k5_linshi)=u(2);
        if(t_xiao_tuiyan==0)%由于第一个周期数据量较少，因此t_xiaotuiyan需要多加一个轨道周期的时间,以此确保可以访问到大碎片的轨道数据
            t_xiao_tuiyan=t_xiao_tuiyan+(ncaiyang+1)*2*pi*sqrt(u(1)^3/miu_earth);
            flag_1=1;
        else%否则的话，一个轨道周期一个轨道周期地进行递推
            t_xiao_tuiyan=t_xiao_tuiyan+2*pi*sqrt(u(1)^3/miu_earth);
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
        t2_jz_no=erfen(t_da_caiyang,1,t_da_caiyang_lie,t2_jz);
        t1_jz_no=erfen(t_da_caiyang,1,t2_jz_no,t1_jz);
        [lie_a_choose,lie_e_choose]=sousuo(aeda_sm(2:4,t1_jz_no:t2_jz_no),u(1),u(2),u(3));%这里的lie_*_choose是部分矩阵的列号
        lie_a_choose=lie_a_choose+t1_jz_no-1;
        lie_e_choose=lie_e_choose+t1_jz_no-1;
        if(flag_1==1)
            v(1)=u(1)+aeda_sm(5,lie_a_choose)*(ncaiyang+1)*2*pi*sqrt(u(1)^3/miu_earth)*s_mx(i);
            v(2)=u(2)+aeda_sm(6,lie_e_choose)*(ncaiyang+1)*2*pi*sqrt(u(1)^3/miu_earth)*s_mx(i);
            v(3)=u(3)+3/4*1.083e-3*(6378e3/u(1)/(1-u(2))^2)^2*sqrt(miu_earth/u(1)^3)*(5*dot(cross(r0,v0)/norm(cross(r0,v0)),[0;0;1])^2-1)*(ncaiyang+1)*2*pi*sqrt(u(1)^3/miu_earth);
            if(v(2)<0)
                v(2)=0;
            else
                if(v(2)>1)
                    v(2)=1;
                end
            end
        else
            v(1)=u(1)+aeda_sm(5,lie_a_choose)*2*pi*sqrt(u(1)^3/miu_earth)*s_mx(i);
            v(2)=u(2)+aeda_sm(6,lie_e_choose)*2*pi*sqrt(u(1)^3/miu_earth)*s_mx(i);
            v(3)=u(3)+3/4*1.083e-3*(6378e3/u(1)/(1-u(2))^2)^2*sqrt(miu_earth/u(1)^3)*(5*dot(cross(r0,v0)/norm(cross(r0,v0)),[0;0;1])^2-1)*2*pi*sqrt(u(1)^3/miu_earth);
            if(v(2)<0)
                v(2)=0;
            else
                if(v(2)>1)
                    v(2)=1;
                end
            end
        end
        u=v;
        if(u(1)*(1-u(2))<=6378e3)%坠落，之后的数据无需记载
            k5_linshi=k5_linshi+1;
            tae_xiao_tuiyan_linshi(1,k5_linshi)=t_xiao_tuiyan;
            tae_xiao_tuiyan_linshi(2,k5_linshi)=-1;
            tae_xiao_tuiyan_linshi(3,k5_linshi)=-1;
            break;
        end
    end
    k5(i)=k5_linshi;
    for j=1:nxh
        if(j>k5_linshi)
            break;
        end
        tae_xiao_tuiyan(:,j,i)=tae_xiao_tuiyan_linshi(:,j);
    end
end
if(sum(k5>=nxh)~=0)
    fprintf('error!nxh');
    pause;
end
t4=toc;






%计算空间分布，并按照不同特征尺度绘图
tic;
%首先绘制四个长度区间参考值随时间的变化
huitu_jiange=max(1,floor(ns/400));
lo_fenbu=6378e3+200e3;
hi_fenbu=6378e3+2*h0-200e3;
dR_fenbu=1e3;
L1_fenbu_cankao=zeros(2,(hi_fenbu-lo_fenbu)/dR_fenbu,ns+1);
L2_fenbu_cankao=zeros(2,(hi_fenbu-lo_fenbu)/dR_fenbu,ns+1);
L3_fenbu_cankao=zeros(2,(hi_fenbu-lo_fenbu)/dR_fenbu,ns+1);
L4_fenbu_cankao=zeros(2,(hi_fenbu-lo_fenbu)/dR_fenbu,ns+1);
P4_shang_cankao=zeros(4,ns+1);
P4_xia_cankao=zeros(4,ns+1);
for i=1:ns+1
    [L1_fenbu_cankao(1,:,i),L1_fenbu_cankao(2,:,i),P4_shang_cankao(1,i),P4_xia_cankao(1,i)]=kongjianfenbu(lo_fenbu,hi_fenbu,dR_fenbu,ae_cankao_L1(i,1:k1(i),1),ae_cankao_L1(i,1:k1(i),2));
    L1_fenbu_cankao(1,:,i)=L1_fenbu_cankao(1,:,i)/no_da;
    P4_shang_cankao(1,i)=P4_shang_cankao(1,i)/no_da;
    P4_xia_cankao(1,i)=(P4_xia_cankao(1,i)+no_da-k1(i))/no_da;
    [L2_fenbu_cankao(1,:,i),L2_fenbu_cankao(2,:,i),P4_shang_cankao(2,i),P4_xia_cankao(2,i)]=kongjianfenbu(lo_fenbu,hi_fenbu,dR_fenbu,ae_cankao_L2(i,1:k2(i),1),ae_cankao_L2(i,1:k2(i),2));
    L2_fenbu_cankao(1,:,i)=L2_fenbu_cankao(1,:,i)/Lcx_1_no_chu;
    P4_shang_cankao(2,i)=P4_shang_cankao(2,i)/Lcx_1_no_chu;
    P4_xia_cankao(2,i)=(P4_xia_cankao(2,i)+Lcx_1_no_chu-k2(i))/Lcx_1_no_chu;
    [L3_fenbu_cankao(1,:,i),L3_fenbu_cankao(2,:,i),P4_shang_cankao(3,i),P4_xia_cankao(3,i)]=kongjianfenbu(lo_fenbu,hi_fenbu,dR_fenbu,ae_cankao_L3(i,1:k3(i),1),ae_cankao_L3(i,1:k3(i),2));
    L3_fenbu_cankao(1,:,i)=L3_fenbu_cankao(1,:,i)/(Lcx_2_no_chu-Lcx_1_no_chu);
    P4_shang_cankao(3,i)=P4_shang_cankao(3,i)/(Lcx_2_no_chu-Lcx_1_no_chu);
    P4_xia_cankao(3,i)=(P4_xia_cankao(3,i)+Lcx_2_no_chu-Lcx_1_no_chu-k3(i))/(Lcx_2_no_chu-Lcx_1_no_chu);
    [L4_fenbu_cankao(1,:,i),L4_fenbu_cankao(2,:,i),P4_shang_cankao(4,i),P4_xia_cankao(4,i)]=kongjianfenbu(lo_fenbu,hi_fenbu,dR_fenbu,ae_cankao_L4(i,1:k4(i),1),ae_cankao_L4(i,1:k4(i),2));
    L4_fenbu_cankao(1,:,i)=L4_fenbu_cankao(1,:,i)/(no_x-Lcx_2_no_chu);
    P4_shang_cankao(4,i)=P4_shang_cankao(4,i)/(no_x-Lcx_2_no_chu);
    P4_xia_cankao(4,i)=(P4_xia_cankao(4,i)+no_x-Lcx_2_no_chu-k4(i))/(no_x-Lcx_2_no_chu);
end
flag_cankao=0;
pic_num=1;
while 1
    if(pic_num>ns+1)
        break;
    end
    f_tu=figure(1);
    if(pic_num==1)
        clf;
    end
    plot(L1_fenbu_cankao(1,:,pic_num),L1_fenbu_cankao(2,:,pic_num),L2_fenbu_cankao(1,:,pic_num),...
        L2_fenbu_cankao(2,:,pic_num),L3_fenbu_cankao(1,:,pic_num),L3_fenbu_cankao(2,:,pic_num),L4_fenbu_cankao(1,:,pic_num),L4_fenbu_cankao(2,:,pic_num));
    ylabel('轨道半径（km）');
    xlabel('位于对应轨道半径区间的概率');
    legend('碎片特征长度大于10cm的参考值',...
        '碎片特征长度在3cm到10cm之间的参考值',...
        '碎片特征长度在1.8cm到3cm之间的参考值',...
        '碎片特征长度在1cm到1.8cm之间的参考值');
    title({['演化时间=',num2str((pic_num-1)*Ts/86400),'天'];...
        ['碎片特征长度大于10cm的参考值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(P4_shang_cankao(1,pic_num)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(P4_xia_cankao(1,pic_num))];...
        ['碎片特征长度在3cm到10cm之间的参考值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(P4_shang_cankao(2,pic_num)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(P4_xia_cankao(2,pic_num))];...
        ['碎片特征长度在1.8cm到3cm之间的参考值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(P4_shang_cankao(3,pic_num)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(P4_xia_cankao(3,pic_num))];...
        ['碎片特征长度在1cm到1.8cm之间的参考值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(P4_shang_cankao(4,pic_num)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(P4_xia_cankao(4,pic_num))]});
    set(gca,'FontSize',16);
    % if(flag_cankao==0)
    %     axis([0 0.03 lo_fenbu/1e3 hi_fenbu/1e3]);
    % else
    %     axis([0 0.01 lo_fenbu/1e3 hi_fenbu/1e3]);
    % end
    axis([0 0.01 lo_fenbu/1e3 hi_fenbu/1e3]);
    set(f_tu, 'unit', 'normalized', 'position', [0,0,1,1]);
    drawnow;
    F_tu=getframe(gcf);
    I_tu=frame2im(F_tu);
    [I_tu,map]=rgb2ind(I_tu,256);
    if(pic_num == 1)
        imwrite(I_tu,map,['zhongzi=',num2str(zhongzi),',ns=',num2str(ns),',cankao_duo.gif'],'gif', 'Loopcount',inf,'DelayTime',0.05);
        clf;
    else
        if(pic_num+huitu_jiange>ns+1)
            if(flag_cankao==0)
                imwrite(I_tu,map,['zhongzi=',num2str(zhongzi),',ns=',num2str(ns),',cankao_duo.gif'],'gif','WriteMode','append','DelayTime',2);
                clf;
            else
                imwrite(I_tu,map,['zhongzi=',num2str(zhongzi),',ns=',num2str(ns),',cankao_duo.gif'],'gif','WriteMode','append','DelayTime',5);
            end
        else
            imwrite(I_tu,map,['zhongzi=',num2str(zhongzi),',ns=',num2str(ns),',cankao_duo.gif'],'gif','WriteMode','append','DelayTime',0.05);
            clf;
        end
    end
    if(flag_cankao==0&&pic_num+huitu_jiange>ns+1)
        flag_cankao=1;
        pic_num=pic_num-huitu_jiange;%特殊情况，不是错误，是因为需要重复绘制最后一幅图
    end
    pic_num=pic_num+huitu_jiange;
end
t5=toc;





%绘制小碎片参考值和推演值随时间变化的对比图
tic;
xiao_cankao=zeros(2,(hi_fenbu-lo_fenbu)/dR_fenbu,ns+1);
Px_cankao=zeros(ns+1,2);
for i=1:ns+1
    [xiao_cankao(1,:,i),xiao_cankao(2,:,i),Px_cankao(i,1),Px_cankao(i,2)]=kongjianfenbu(lo_fenbu,hi_fenbu,dR_fenbu,...
        [ae_cankao_L2(i,1:k2(i),1),ae_cankao_L3(i,1:k3(i),1),ae_cankao_L4(i,1:k4(i),1)],[ae_cankao_L2(i,1:k2(i),2),ae_cankao_L3(i,1:k3(i),2),ae_cankao_L4(i,1:k4(i),2)]);
    xiao_cankao(1,:,i)=xiao_cankao(1,:,i)/no_x;
    Px_cankao(i,1)=Px_cankao(i,1)/no_x;
    Px_cankao(i,2)=(Px_cankao(i,2)+no_x-k2(i)-k3(i)-k4(i))/no_x;
end
xiao_tuiyan=zeros(2,(hi_fenbu-lo_fenbu)/dR_fenbu,ns+1);
Px_tuiyan=zeros(ns+1,2);
xiao_tuiyan_shuliang=zeros(ns+1,1);
xiao_tuiyan_a=zeros(1,no_x);
xiao_tuiyan_e=zeros(1,no_x);
for i=1:ns+1
    for j=1:no_x
        zuijin=erfen(tae_xiao_tuiyan(1,1:k5(j),j),1,k5(j),(i-1)*Ts);
        if(tae_xiao_tuiyan(2,zuijin,j)>0)
            xiao_tuiyan_shuliang(i)=xiao_tuiyan_shuliang(i)+1;
            xiao_tuiyan_a(xiao_tuiyan_shuliang(i))=tae_xiao_tuiyan(2,zuijin,j);
            xiao_tuiyan_e(xiao_tuiyan_shuliang(i))=tae_xiao_tuiyan(3,zuijin,j);
        end
    end
    [xiao_tuiyan(1,:,i),xiao_tuiyan(2,:,i),Px_tuiyan(i,1),Px_tuiyan(i,2)]=kongjianfenbu(lo_fenbu,hi_fenbu,dR_fenbu,xiao_tuiyan_a(1:xiao_tuiyan_shuliang(i)),xiao_tuiyan_e(1:xiao_tuiyan_shuliang(i)));
    xiao_tuiyan(1,:,i)=xiao_tuiyan(1,:,i)/no_x;
    Px_tuiyan(i,1)=Px_tuiyan(i,1)/no_x;
    Px_tuiyan(i,2)=(Px_tuiyan(i,2)+no_x-xiao_tuiyan_shuliang(i))/no_x;
end
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
        %在服务器上面跑的时候，加了注释后面那一段，但这是不对的，所以载入数据之后需要重新计算
        JS_sandu(i)=JS_sandu(i)+0.5*Px_tuiyan(i,2)*log(Px_tuiyan(i,2)*2/(Px_cankao(i,2)+Px_tuiyan(i,2)));%+0.5*Px_cankao(i,2)*log(Px_cankao(i,2)*2/(Px_cankao(i,2)+Px_tuiyan(i,2)))
    end
end
flag_tuiyan=0;
pic_num=1;
while 1
    if(pic_num>ns+1)
        break;
    end
    f_tu=figure(2);
    if(pic_num==1)
        clf;
    end
    plot(xiao_cankao(1,:,pic_num),xiao_cankao(2,:,pic_num),'blue',xiao_tuiyan(1,:,pic_num),xiao_tuiyan(2,:,pic_num),'magenta');
    ylabel('轨道半径（km）');
    xlabel('位于对应轨道半径区间的概率');
    legend('小碎片（特征长度小于10cm）的参考值','小碎片（特征长度小于10cm）的推演值');
    title({['演化时间=',num2str((pic_num-1)*Ts/86400),'天'];['JS散度=',num2str(JS_sandu(pic_num))];...
        ['参考值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(Px_cankao(pic_num,1)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(Px_cankao(pic_num,2))];...
        ['推演值：',num2str(hi_fenbu/1e3),'km之上的概率=',num2str(Px_tuiyan(pic_num,1)),'；',num2str(lo_fenbu/1e3),'km之下的概率=',num2str(Px_tuiyan(pic_num,2))]});
    set(gca,'FontSize',16);
    % if(flag_tuiyan==0)
    %     axis([0 0.03 lo_fenbu/1e3 hi_fenbu/1e3]);
    % else
    %     axis([0 0.01 lo_fenbu/1e3 hi_fenbu/1e3]);
    % end
    axis([0 0.01 lo_fenbu/1e3 hi_fenbu/1e3]);
    set(f_tu, 'unit', 'normalized', 'position', [0,0,1,1]);
    drawnow;
    F_tu=getframe(gcf);
    I_tu=frame2im(F_tu);
    [I_tu,map]=rgb2ind(I_tu,256);
    if(pic_num == 1)
        imwrite(I_tu,map,['zhongzi=',num2str(zhongzi),',ns=',num2str(ns),',xiao_cankao_tuiyan_ncaiyang=',num2str(ncaiyang),'_baoliu=',num2str(baoliu),'.gif'],'gif', 'Loopcount',inf,'DelayTime',0.05);
        clf;
    else
        if(pic_num+huitu_jiange>ns+1)
            if(flag_tuiyan==0)
                imwrite(I_tu,map,['zhongzi=',num2str(zhongzi),',ns=',num2str(ns),',xiao_cankao_tuiyan_ncaiyang=',num2str(ncaiyang),'_baoliu=',num2str(baoliu),'.gif'],'gif','WriteMode','append','DelayTime',2);
                clf;
            else
                imwrite(I_tu,map,['zhongzi=',num2str(zhongzi),',ns=',num2str(ns),',xiao_cankao_tuiyan_ncaiyang=',num2str(ncaiyang),'_baoliu=',num2str(baoliu),'.gif'],'gif','WriteMode','append','DelayTime',5);
            end
        else
            imwrite(I_tu,map,['zhongzi=',num2str(zhongzi),',ns=',num2str(ns),',xiao_cankao_tuiyan_ncaiyang=',num2str(ncaiyang),'_baoliu=',num2str(baoliu),'.gif'],'WriteMode','append','DelayTime',0.05);
            clf;
        end
    end
    if(flag_tuiyan==0&&pic_num+huitu_jiange>ns+1)
        flag_tuiyan=1;
        pic_num=pic_num-huitu_jiange;%特殊情况，不是错误，是因为需要重复绘制最后一幅图
    end
    pic_num=pic_num+huitu_jiange;
end
t6=toc;

vars = whos;  % 获取所有变量信息
figVars = {vars(strcmp({vars.class}, 'matlab.ui.Figure')).name};  % 找出figure变量名

% 要保存的变量列表（排除figure）
saveVars = setdiff({vars.name}, figVars);

% 构造文件名，注意加上扩展名.mat
filename = sprintf('zhongzi=%d_ns=%d.mat', zhongzi, ns);

% 保存变量，使用 -v7.3 支持大变量
save(filename, saveVars{:}, '-v7.3');
end


function [value,isterminal,direction]=myEventsFcn_xiao(~,y)
value=norm(y(1:3))-6378e3;
isterminal=1;
direction=-1;
end

function [value,isterminal,direction]=myEventsFcn_da(~,y)
value=[norm(y(1:3))-6378e3...
    ,y(1)+1];
isterminal=[1,0];
direction=[-1,-1];
end

function xuhao=erfen(xulie,lo,hi,zhi)%返回序号对应的值<=查找值，序号+1对应的值大于查找值
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

function [yspace,xspace,P_shang,P_xia] = kongjianfenbu(lo,hi,dR,a,e)
%lo和hi分别为半径的上下限，dR为间隔，a、e为1维行向量，no为向量维度
%yspace大致代表着空间碎片在半径方向的一个分布
xspace=linspace(lo,hi-dR,(hi-lo)/dR)+dR/2;
yspace=zeros(1,(hi-lo)/dR);
[~,lie_a]=size(a);
[~,lie_e]=size(e);
if(lie_a~=lie_e)
    fprintf("error!kongjianfenbu");
    pause;
end
P_shang=0;
P_xia=0;
for i=1:lie_a
    for j=1:(hi-lo)/dR
        r_down=lo+(j-1)*dR;
        r_up=lo+j*dR;
        if(a(i)*(1+e(i))>r_down&&a(i)*(1-e(i))<r_up)
            jifen_down=asin(max((r_down-a(i))/a(i)/e(i),-1));
            jifen_shang=asin(min((r_up-a(i))/a(i)/e(i),1));
            yspace(j)=yspace(j)+((jifen_shang-jifen_down)-e(i)*(cos(jifen_shang)-cos(jifen_down)))/pi;
        end
    end
    if(a(i)*(1+e(i))>hi)
        jifen_down=asin(max((hi-a(i))/a(i)/e(i),-1));
        jifen_shang=asin(1);
        P_shang=P_shang+((jifen_shang-jifen_down)-e(i)*(cos(jifen_shang)-cos(jifen_down)))/pi;
    end
    if(a(i)*(1-e(i))<lo)
        jifen_down=asin(-1);
        jifen_shang=asin(min((lo-a(i))/a(i)/e(i),1));
        P_xia=P_xia+((jifen_shang-jifen_down)-e(i)*(cos(jifen_shang)-cos(jifen_down)))/pi;
    end
end
xspace=xspace/1e3;%单位变为1km
end