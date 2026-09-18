function results=zhuchengxu_sensitivity(zhongzi,ns)
% Comment 4: parameter-sensitivity analysis for the closed numerical experiment
%
% 本程序直接读取 zhuchengxu.m 已经保存的计算结果，不重新进行大、小碎片的
% 数值积分。仅从“大碎片轨道数据归并排序”之后重新执行 O2UDI 推演，并改变：
% 1) historical-data window T_hist;
% 2) first-stage screening fraction eta_r;
% 3) inference interval delta_t;
% 4) radial-bin width dR;
% 5) fragment-specific drag-coefficient assumption.
%
% 每一类参数均只改变一个量，其余设置保持 baseline 不变。baseline 直接复用
% 原始保存的 O2UDI 推演结果，不重复计算。输出每组的 time-averaged JS
% divergence 与 maximum JS divergence。
%
% 建议调用：
%   zhuchengxu_sensitivity(2,21000)

%% ========================= 参数设置 =========================
data_dir='data';
filename=fullfile(data_dir,sprintf('zhongzi=%d_ns=%d.mat',zhongzi,ns));
if(~isfile(filename))
    error('Cannot find saved data file: %s',filename);
end

% 五组 sensitivity 参数
% T_hist 采用大碎片采样间隔的 1~5 倍。原程序 baseline 中 baoliu=2，
% 因此 baseline 精确对应 2*ncaiyang*Ts；论文中四舍五入写为 16.8 h。
T_hist_multiple_values=[1,2,3,4,5];
eta_values=[0.01,0.015,0.02,0.03,0.05];
delta_multiple_values=[1,2,3,4,5];
dR_values_km=[0.5,1,2,3,5];
CD_uncertainty_values=[0,0.025,0.05,0.075,0.10];

% baseline 设置
T_hist_baseline_multiple=2;
eta_baseline=0.02;
delta_multiple_baseline=1;
dR_baseline_km=1;
CD_uncertainty_baseline=0;

% 大碎片每隔 ncaiyang 个轨道周期采样一次，与原程序保持一致
ncaiyang=5;

% CD 随机扰动固定种子。每颗碎片只生成一次 z，五种扰动幅度共用同一组 z。
CD_random_seed=20260917;

% 是否保存完整的 JS 时间序列
save_full_JS=true;

%% ========================= 读取原有结果 =========================
load(filename,'Ts','ts','miu_earth','no_da','no_x','taewda','taewda_shu',...
    's_m','s_mx','rcx','vcx','r0','v0',...
    'xiao_cankao','Px_cankao','xiao_tuiyan','Px_tuiyan',...
    'tae_xiao_tuiyan','k5',...
    'lo_fenbu','hi_fenbu','dR_fenbu');

fprintf('\n================ Comment 4 sensitivity =================\n');
fprintf('Saved data: %s\n',filename);
fprintf('Fragments: observable = %d, unobservable = %d\n',no_da,no_x);
fprintf('Ts = %.6f h\n',Ts/3600);
fprintf('Observable sampling interval = %.6f h\n',ncaiyang*Ts/3600);

% T_hist 使用采样间隔的整数倍进行 sensitivity。这样 baseline 与原程序
% baoliu=2 完全一致，而论文中的 16.8 h 只是对精确值的四舍五入。
T_hist_values_s=T_hist_multiple_values*ncaiyang*Ts;
T_hist_values_h=T_hist_values_s/3600;
T_hist_baseline_s=T_hist_baseline_multiple*ncaiyang*Ts;
T_hist_baseline_h=T_hist_baseline_s/3600;
T_hist_original_code_h=T_hist_baseline_h;

fprintf('Historical-window settings used in calculation (h): ');
fprintf('%.6f ',T_hist_values_h);
fprintf('\n');
fprintf('Manuscript-rounded settings (h): ');
fprintf('%.1f ',T_hist_values_h);
fprintf('\n');
fprintf('Baseline T_hist = %.6f h (reported as %.1f h)\n',...
    T_hist_baseline_h,T_hist_baseline_h);

% baseline 直接使用原始 O2UDI 结果，不重复推演；只重新计算 JS，
% 以统一边界概率项的计算方式。
JS_saved_baseline=jisuan_JS_from_distribution(xiao_cankao,Px_cankao,...
    xiao_tuiyan,Px_tuiyan);
JS_baseline=JS_saved_baseline;
fprintf('Saved original baseline: mean JS = %.8f, max JS = %.8f\n',...
    mean(JS_baseline),max(JS_baseline));

clear xiao_tuiyan Px_tuiyan;

%% ========================= 大碎片记录归并 =========================
% 与 zhuchengxu.m 相同，只增加第8行记录 observable fragment ID，
% 用于 CD sensitivity 中读取该匹配大碎片的随机 CD。
tic;
aeda_sm=zeros(8,sum(ceil((taewda_shu-1)/ncaiyang))+no_da)-2;
k=1;
dangqian=zeros(no_da,1)+1+ncaiyang;

while(sum(dangqian>taewda_shu)~=no_da)
    tda_min=ts+1;
    i_choose=1;
    for i=1:no_da
        if(dangqian(i)<=taewda_shu(i)&&taewda(1,dangqian(i),i)<tda_min)
            i_choose=i;
            tda_min=taewda(1,dangqian(i),i);
        end
    end

    aeda_sm(1,k)=taewda(1,dangqian(i_choose),i_choose);
    aeda_sm(2,k)=taewda(2,dangqian(i_choose)-ncaiyang,i_choose);
    aeda_sm(3,k)=taewda(3,dangqian(i_choose)-ncaiyang,i_choose);
    aeda_sm(4,k)=taewda(4,dangqian(i_choose)-ncaiyang,i_choose);
    aeda_sm(5,k)=(taewda(2,dangqian(i_choose),i_choose)-taewda(2,dangqian(i_choose)-ncaiyang,i_choose))...
        /(taewda(1,dangqian(i_choose),i_choose)-taewda(1,dangqian(i_choose)-ncaiyang,i_choose))/s_m(i_choose);
    aeda_sm(6,k)=(taewda(3,dangqian(i_choose),i_choose)-taewda(3,dangqian(i_choose)-ncaiyang,i_choose))...
        /(taewda(1,dangqian(i_choose),i_choose)-taewda(1,dangqian(i_choose)-ncaiyang,i_choose))/s_m(i_choose);
    aeda_sm(7,k)=s_m(i_choose);
    aeda_sm(8,k)=i_choose;

    dangqian(i_choose)=dangqian(i_choose)+ncaiyang;
    k=k+1;
end

[~,aeda_sm_lie]=size(aeda_sm);
for i=aeda_sm_lie:-1:1
    if(aeda_sm(1,i)<0)
        aeda_sm(:,i)=[];
    end
end

fprintf('Large-debris record merging finished: %.2f s, records = %d\n',...
    toc,size(aeda_sm,2));

clear taewda taewda_shu s_m dangqian;

%% ========================= CD 固定随机扰动 =========================
rng(CD_random_seed,'twister');
z_da=2*rand(no_da,1)-1;
z_xiao=2*rand(no_x,1)-1;

%% ========================= 并行池 =========================
use_thread_pool=prepare_parallel_pool_sensitivity();

%% ========================= screening-fraction sensitivity =========================
fprintf('\n---------------- eta_r sensitivity ----------------\n');
JS_eta=zeros(length(eta_values),ns+1);
for icase=1:length(eta_values)
    value=eta_values(icase);

    if(abs(value-eta_baseline)<1e-15)
        JS_eta(icase,:)=JS_baseline;
    else
        [tae_case,k5_case]=run_inference_case(aeda_sm,rcx,vcx,s_mx,...
            no_x,ts,miu_earth,r0,v0,ncaiyang,Ts,...
            T_hist_baseline_s,value,delta_multiple_baseline,...
            CD_uncertainty_baseline,z_da,z_xiao);

        JS_eta(icase,:)=jisuan_JS_from_tae_sensitivity(tae_case,k5_case,...
            xiao_cankao,Px_cankao,lo_fenbu,hi_fenbu,dR_fenbu,no_x,ns,Ts,use_thread_pool);
        clear tae_case k5_case;
    end

    fprintf('eta_r = %6.2f%%: mean JS = %.8f, max JS = %.8f\n',...
        value*100,mean(JS_eta(icase,:)),max(JS_eta(icase,:)));
end

%% ========================= radial-bin width sensitivity =========================
% 为避免对 0.5/1/2/3/5 km 分别重复计算，先在最细的 0.5 km 网格上
% 重构一次 reference 和 inferred 分布，再将相邻 0.5 km bins 相加得到
% 其余四种 bin width。由于各网格边界完全对齐，该处理与直接使用对应
% 较粗网格积分等价。
fprintf('\n---------------- radial-bin width sensitivity ----------------\n');
load(filename,'ae_cankao_L2','ae_cankao_L3','ae_cankao_L4','k2','k3','k4');

fine_dR=min(dR_values_km)*1e3;
[P_ref_fine,P_inf_fine,Px_ref_fine,Px_inf_fine]=...
    compute_fine_distributions(ae_cankao_L2,ae_cankao_L3,ae_cankao_L4,...
    k2,k3,k4,tae_xiao_tuiyan,k5,no_x,ns,Ts,...
    lo_fenbu,hi_fenbu,fine_dR,use_thread_pool);

JS_dR=zeros(length(dR_values_km),ns+1);
for icase=1:length(dR_values_km)
    factor=round(dR_values_km(icase)/(fine_dR/1e3));
    if(abs(factor*fine_dR/1e3-dR_values_km(icase))>1e-12)
        error('Radial-bin width %.6f km is not an integer multiple of %.6f km.',...
            dR_values_km(icase),fine_dR/1e3);
    end

    P_ref=rebin_probability(P_ref_fine,factor);
    P_inf=rebin_probability(P_inf_fine,factor);
    JS_dR(icase,:)=jisuan_JS_matrix(P_ref,Px_ref_fine,P_inf,Px_inf_fine);

    fprintf('dR = %5.2f km: mean JS = %.8f, max JS = %.8f\n',...
        dR_values_km(icase),mean(JS_dR(icase,:)),max(JS_dR(icase,:)));
end

clear ae_cankao_L2 ae_cankao_L3 ae_cankao_L4 k2 k3 k4;
clear P_ref_fine P_inf_fine Px_ref_fine Px_inf_fine;

% baseline 推演轨迹后续不再需要，释放内存。
clear tae_xiao_tuiyan k5;

%% ========================= historical-data window sensitivity =========================
fprintf('\n---------------- T_hist sensitivity ----------------\n');
JS_Thist=zeros(length(T_hist_multiple_values),ns+1);
for icase=1:length(T_hist_multiple_values)
    multiple=T_hist_multiple_values(icase);
    value_s=T_hist_values_s(icase);
    value_h=T_hist_values_h(icase);

    if(multiple==T_hist_baseline_multiple)
        JS_Thist(icase,:)=JS_baseline;
    else
        [tae_case,k5_case]=run_inference_case(aeda_sm,rcx,vcx,s_mx,...
            no_x,ts,miu_earth,r0,v0,ncaiyang,Ts,...
            value_s,eta_baseline,delta_multiple_baseline,...
            CD_uncertainty_baseline,z_da,z_xiao);

        JS_Thist(icase,:)=jisuan_JS_from_tae_sensitivity(tae_case,k5_case,...
            xiao_cankao,Px_cankao,lo_fenbu,hi_fenbu,dR_fenbu,no_x,ns,Ts,use_thread_pool);
        clear tae_case k5_case;
    end

    fprintf('T_hist = %4.1f h (%d sampling interval%s): mean JS = %.8f, max JS = %.8f\n',...
        value_h,multiple,plural_s(multiple),mean(JS_Thist(icase,:)),max(JS_Thist(icase,:)));
end



%% ========================= inference-interval sensitivity =========================
fprintf('\n---------------- delta_t sensitivity ----------------\n');
JS_delta=zeros(length(delta_multiple_values),ns+1);
delta_values_h=delta_multiple_values*Ts/3600;
for icase=1:length(delta_multiple_values)
    value=delta_multiple_values(icase);

    if(value==delta_multiple_baseline)
        JS_delta(icase,:)=JS_baseline;
    else
        [tae_case,k5_case]=run_inference_case(aeda_sm,rcx,vcx,s_mx,...
            no_x,ts,miu_earth,r0,v0,ncaiyang,Ts,...
            T_hist_baseline_s,eta_baseline,value,...
            CD_uncertainty_baseline,z_da,z_xiao);

        JS_delta(icase,:)=jisuan_JS_from_tae_sensitivity(tae_case,k5_case,...
            xiao_cankao,Px_cankao,lo_fenbu,hi_fenbu,dR_fenbu,no_x,ns,Ts,use_thread_pool);
        clear tae_case k5_case;
    end

    fprintf('delta_t ~= %6.3f h (%d orbit%s): mean JS = %.8f, max JS = %.8f\n',...
        delta_values_h(icase),value,plural_s(value),...
        mean(JS_delta(icase,:)),max(JS_delta(icase,:)));
end

%% ========================= drag-coefficient sensitivity =========================
fprintf('\n---------------- C_D assumption sensitivity ----------------\n');
JS_CD=zeros(length(CD_uncertainty_values),ns+1);
for icase=1:length(CD_uncertainty_values)
    value=CD_uncertainty_values(icase);

    if(abs(value-CD_uncertainty_baseline)<1e-15)
        JS_CD(icase,:)=JS_baseline;
    else
        [tae_case,k5_case]=run_inference_case(aeda_sm,rcx,vcx,s_mx,...
            no_x,ts,miu_earth,r0,v0,ncaiyang,Ts,...
            T_hist_baseline_s,eta_baseline,delta_multiple_baseline,...
            value,z_da,z_xiao);

        JS_CD(icase,:)=jisuan_JS_from_tae_sensitivity(tae_case,k5_case,...
            xiao_cankao,Px_cankao,lo_fenbu,hi_fenbu,dR_fenbu,no_x,ns,Ts,use_thread_pool);
        clear tae_case k5_case;
    end

    fprintf('C_D random range = +/- %5.2f%%: mean JS = %.8f, max JS = %.8f\n',...
        value*100,mean(JS_CD(icase,:)),max(JS_CD(icase,:)));
end

%% ========================= 汇总 =========================
results=struct();
results.T_hist_multiple_values=T_hist_multiple_values;
results.T_hist_values_s=T_hist_values_s;
results.T_hist_values_h=T_hist_values_h;
results.eta_values=eta_values;
results.delta_multiple_values=delta_multiple_values;
results.delta_values_h=delta_values_h;
results.dR_values_km=dR_values_km;
results.CD_uncertainty_values=CD_uncertainty_values;
results.JS_saved_baseline=JS_saved_baseline;
results.JS_baseline=JS_baseline;
results.JS_Thist=JS_Thist;
results.JS_eta=JS_eta;
results.JS_delta=JS_delta;
results.JS_dR=JS_dR;
results.JS_CD=JS_CD;

results.mean_Thist=mean(JS_Thist,2);
results.max_Thist=max(JS_Thist,[],2);
results.mean_eta=mean(JS_eta,2);
results.max_eta=max(JS_eta,[],2);
results.mean_delta=mean(JS_delta,2);
results.max_delta=max(JS_delta,[],2);
results.mean_dR=mean(JS_dR,2);
results.max_dR=max(JS_dR,[],2);
results.mean_CD=mean(JS_CD,2);
results.max_CD=max(JS_CD,[],2);

% 生成统一 CSV 表，每类参数5行，共25行。
Parameter=[repmat("Historical-data window T_hist",5,1);...
    repmat("Screening fraction eta_r",5,1);...
    repmat("Inference interval delta_t",5,1);...
    repmat("Radial-bin width",5,1);...
    repmat("C_D assumption",5,1)];

Setting=[compose("%.1f h",T_hist_values_h(:));...
    compose("%.2f%%",100*eta_values(:));...
    compose("%.3f h",delta_values_h(:));...
    compose("%.2f km",dR_values_km(:));...
    compose("+/- %.2f%%",100*CD_uncertainty_values(:))];

Mean_JS=[results.mean_Thist;results.mean_eta;results.mean_delta;...
    results.mean_dR;results.mean_CD];
Maximum_JS=[results.max_Thist;results.max_eta;results.max_delta;...
    results.max_dR;results.max_CD];

summary_table=table(Parameter,Setting,Mean_JS,Maximum_JS);
results.summary_table=summary_table;

disp(summary_table);

csv_filename=fullfile(data_dir,sprintf('zhongzi=%d_ns=%d_comment4_sensitivity.csv',zhongzi,ns));
writetable(summary_table,csv_filename);

mat_filename=fullfile(data_dir,sprintf('zhongzi=%d_ns=%d_comment4_sensitivity.mat',zhongzi,ns));
if(save_full_JS)
    save(mat_filename,'results','zhongzi','ns','Ts','T_hist_original_code_h',...
        'T_hist_baseline_multiple','T_hist_baseline_s','T_hist_baseline_h',...
        'eta_baseline','delta_multiple_baseline',...
        'dR_baseline_km','CD_uncertainty_baseline','CD_random_seed','-v7.3');
else
    results_for_save=rmfield(results,{'JS_saved_baseline','JS_baseline','JS_Thist',...
        'JS_eta','JS_delta','JS_dR','JS_CD'});
    save(mat_filename,'results_for_save','zhongzi','ns','Ts','T_hist_original_code_h',...
        'T_hist_baseline_multiple','T_hist_baseline_s','T_hist_baseline_h',...
        'eta_baseline','delta_multiple_baseline',...
        'dR_baseline_km','CD_uncertainty_baseline','CD_random_seed','-v7.3');
end

fprintf('\nCSV saved to: %s\n',csv_filename);
fprintf('MAT saved to: %s\n',mat_filename);
fprintf('==========================================================\n');
end


%% ========================================================================
function [tae_xiao_tuiyan,k5]=run_inference_case(aeda_sm,rcx,vcx,s_mx,...
    no_x,ts,miu_earth,r0,v0,ncaiyang,Ts,T_hist,eta_r,delta_multiple,...
    CD_uncertainty,z_da,z_xiao)
% 从已经保存的大碎片轨道记录重新执行 O2UDI，不做任何数值积分。
%
% T_hist: historical-data window, s
% eta_r: first-stage screening fraction
% delta_multiple: inference interval expressed as an integer multiple of
%                 the instantaneous orbital period
% CD_uncertainty: 0~0.1, each fragment uses C_D=2.2*(1+alpha*z)

fprintf('\nStart inference case: T_hist=%.3f h, eta_r=%.3f%%, delta=%d orbit(s), CD=+/-%.3f%%\n',...
    T_hist/3600,eta_r*100,delta_multiple,CD_uncertainty*100);

tic;
t_da_caiyang=aeda_sm(1,:);
t_da_caiyang_lie=length(t_da_caiyang);
if(t_da_caiyang(end)<0)
    error('error!caiyang_guibing');
end

max_tuiyan_shu=ceil(ts/(2*pi*sqrt(6378000^3/miu_earth)))+10;
tae_xiao_tuiyan=zeros(3,max_tuiyan_shu,no_x);
k5=zeros(no_x,1);

CD_da=2.2*(1+CD_uncertainty*z_da);
CD_xiao=2.2*(1+CD_uncertainty*z_xiao);

parfor i=1:no_x
    u=zeros(3,1);
    v=zeros(3,1);
    [u(1),u(2),u(3)]=r0v0_genshu(rcx(:,i),vcx(:,i),miu_earth);

    t_xiao_tuiyan=0;
    k5_linshi=0;
    tae_linshi=zeros(3,max_tuiyan_shu);

    while 1
        k5_linshi=k5_linshi+1;
        tae_linshi(1,k5_linshi)=t_xiao_tuiyan;
        tae_linshi(2,k5_linshi)=u(1);
        tae_linshi(3,k5_linshi)=u(2);

        T_orbit=2*pi*sqrt(u(1)^3/miu_earth);
        if(t_xiao_tuiyan==0)
            % 与原程序 baseline 的 6T 首步一致；对于 mT sensitivity，
            % 首步推广为 (ncaiyang+m)T，从而保证第一次匹配时已有 observable record。
            dt_tuiyan=(ncaiyang+delta_multiple)*T_orbit;
            t_xiao_tuiyan=t_xiao_tuiyan+dt_tuiyan;
        else
            dt_tuiyan=delta_multiple*T_orbit;
            t_xiao_tuiyan=t_xiao_tuiyan+dt_tuiyan;
        end

        if(t_xiao_tuiyan>ts)
            break;
        end
        if(t_xiao_tuiyan<t_da_caiyang(1))
            error('error!t_xiao_tuiyan');
        end

        t1_jz=max(0,t_xiao_tuiyan-T_hist);
        t2_jz=t_xiao_tuiyan;
        t2_jz_no=erfen_sensitivity(t_da_caiyang,1,t_da_caiyang_lie,t2_jz);
        t1_jz_no=erfen_sensitivity(t_da_caiyang,1,t2_jz_no,t1_jz);

        [lie_a_choose,lie_e_choose]=sousuo_sensitivity(...
            aeda_sm(2:4,t1_jz_no:t2_jz_no),u(1),u(2),u(3),eta_r);

        lie_a_choose=lie_a_choose+t1_jz_no-1;
        lie_e_choose=lie_e_choose+t1_jz_no-1;

        % Drag-coefficient sensitivity:
        % reference trajectories remain generated with C_D=2.2;
        % only the assumed transfer ratio C_D,u / C_D,o is perturbed.
        id_da_a=round(aeda_sm(8,lie_a_choose));
        id_da_e=round(aeda_sm(8,lie_e_choose));
        CD_ratio_a=CD_xiao(i)/CD_da(id_da_a);
        CD_ratio_e=CD_xiao(i)/CD_da(id_da_e);

        v(1)=u(1)+aeda_sm(5,lie_a_choose)*dt_tuiyan*s_mx(i)*CD_ratio_a;
        v(2)=u(2)+aeda_sm(6,lie_e_choose)*dt_tuiyan*s_mx(i)*CD_ratio_e;
        v(3)=u(3)+3/4*1.083e-3*(6378e3/u(1)/(1-u(2))^2)^2*sqrt(miu_earth/u(1)^3)...
            *(5*dot(cross(r0,v0)/norm(cross(r0,v0)),[0;0;1])^2-1)*dt_tuiyan;

        if(v(2)<0)
            v(2)=0;
        elseif(v(2)>1)
            v(2)=1;
        end

        u=v;
        if(u(1)*(1-u(2))<=6378e3)
            k5_linshi=k5_linshi+1;
            tae_linshi(1,k5_linshi)=t_xiao_tuiyan;
            tae_linshi(2,k5_linshi)=-1;
            tae_linshi(3,k5_linshi)=-1;
            break;
        end
    end

    k5(i)=k5_linshi;
    tae_xiao_tuiyan(:,:,i)=tae_linshi;
end

fprintf('Inference finished: %.2f h\n',toc/3600);
end


%% ========================================================================
function [lie_a_choose,lie_e_choose]=sousuo_sensitivity(aeda_sm,a,e,w,eta_r)
% 与原 sousuo.m 保持相同的两阶段结构，只将固定的 2%% 改为 eta_r。
[~,max_lie]=size(aeda_sm);
chushai=max(1,ceil(eta_r*max_lie));

E_r=abs(aeda_sm(1,:).*(1+aeda_sm(2,:))-a*(1+e))+...
    abs(aeda_sm(1,:).*(1-aeda_sm(2,:))-a*(1-e));
[~,xuhao]=sort(E_r,'ascend');
cs_xuhao=xuhao(1:chushai);

E_w=abs(aeda_sm(3,cs_xuhao)-w);
E_w=min(E_w,2*pi-E_w);
[~,lie_min]=min(E_w);

if(e>3e-3)
    lie_a_choose=cs_xuhao(lie_min);
    lie_e_choose=cs_xuhao(lie_min);
else
    lie_a_choose=cs_xuhao(1);
    lie_e_choose=cs_xuhao(lie_min);
end
end


%% ========================================================================
function xuhao=erfen_sensitivity(xulie,lo,hi,zhi)
% 返回序号对应的值 <= 查找值，序号+1对应的值 > 查找值。
if(hi<=lo)
    xuhao=lo;
    return;
end
while 1
    mid=floor((lo+hi)/2);
    if(xulie(mid)<=zhi)
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


%% ========================================================================
function JS_sandu=jisuan_JS_from_distribution(xiao_cankao,Px_cankao,...
    xiao_tuiyan,Px_tuiyan)
% 根据已经保存的分布重新计算 JS，避免旧服务器版本边界项差异。
P=squeeze(xiao_cankao(1,:,:)).';
Q=squeeze(xiao_tuiyan(1,:,:)).';
JS_sandu=jisuan_JS_matrix(P,Px_cankao,Q,Px_tuiyan);
end


%% ========================================================================
function JS_sandu=jisuan_JS_matrix(P,Px_P,Q,Px_Q)
% P/Q: epoch x bin
% Px_*: epoch x 2, [upper-bound probability, lower-bound probability]
no_epoch=size(P,1);
JS_sandu=zeros(1,no_epoch);
for i=1:no_epoch
    p=[P(i,:),Px_P(i,1),Px_P(i,2)];
    q=[Q(i,:),Px_Q(i,1),Px_Q(i,2)];
    m=0.5*(p+q);

    idp=(p>0);
    idq=(q>0);
    JS_sandu(i)=0.5*sum(p(idp).*log(p(idp)./m(idp)))+...
        0.5*sum(q(idq).*log(q(idq)./m(idq)));
end
end


%% ========================================================================
function JS_sandu=jisuan_JS_from_tae_sensitivity(tae_xiao_tuiyan,k5,...
    xiao_cankao,Px_cankao,lo_fenbu,hi_fenbu,dR_fenbu,no_x,ns,Ts,use_thread_pool)
% 由一组新的 inferred trajectory 计算 JS。
% 采用“分批提取 + 小数组并行”，避免将完整三维轨迹广播到 worker。

JS_sandu=zeros(1,ns+1);
batch_size=64;
no_epoch=ns+1;
no_batch=ceil(no_epoch/batch_size);

for ibatch=1:no_batch
    i_begin=(ibatch-1)*batch_size+1;
    i_end=min(ibatch*batch_size,no_epoch);
    epoch_no=i_begin:i_end;
    nb=length(epoch_no);
    t_target=(epoch_no-1)*Ts;

    a_batch=zeros(no_x,nb);
    e_batch=zeros(no_x,nb);

    for j=1:no_x
        kj=k5(j);
        if(kj<1)
            continue;
        end

        k_now=erfen_tae_sensitivity(tae_xiao_tuiyan,j,kj,t_target(1));
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

    P_ref_batch=permute(xiao_cankao(1,:,epoch_no),[3,2,1]);
    Px_ref_batch=Px_cankao(epoch_no,:);
    JS_batch=zeros(1,nb);

    if(use_thread_pool)
        parfor b=1:nb
            alive=(a_batch(:,b)>0);
            n_alive=sum(alive);
            [P_inf,~,P_up,P_low]=kongjianfenbu_fast(lo_fenbu,hi_fenbu,dR_fenbu,...
                a_batch(alive,b).',e_batch(alive,b).');
            P_inf=P_inf/no_x;
            Px_inf=[P_up/no_x,(P_low+no_x-n_alive)/no_x];
            JS_batch(b)=one_epoch_JS(P_ref_batch(b,:),Px_ref_batch(b,:),P_inf,Px_inf);
        end
    else
        for b=1:nb
            alive=(a_batch(:,b)>0);
            n_alive=sum(alive);
            [P_inf,~,P_up,P_low]=kongjianfenbu_fast(lo_fenbu,hi_fenbu,dR_fenbu,...
                a_batch(alive,b).',e_batch(alive,b).');
            P_inf=P_inf/no_x;
            Px_inf=[P_up/no_x,(P_low+no_x-n_alive)/no_x];
            JS_batch(b)=one_epoch_JS(P_ref_batch(b,:),Px_ref_batch(b,:),P_inf,Px_inf);
        end
    end

    JS_sandu(epoch_no)=JS_batch;

    if(ibatch==1 || ibatch==no_batch || mod(ibatch,max(1,floor(no_batch/20)))==0)
        fprintf('    JS progress: %d/%d epochs (%.1f%%)\n',...
            i_end,no_epoch,100*i_end/no_epoch);
    end
end
end


%% ========================================================================
function xuhao=erfen_tae_sensitivity(tae_xiao_tuiyan,j,hi,zhi)
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


%% ========================================================================
function [P_ref,P_inf,Px_ref,Px_inf]=compute_fine_distributions(...
    ae2,ae3,ae4,k2,k3,k4,tae_inf,k5_inf,no_x,ns,Ts,lo,hi,dR,use_thread_pool)
% 在最细 radial grid 上一次性计算 baseline reference/inferred 分布。
% 后续更粗的 bin width 直接通过相邻 bins 求和获得。

no_epoch=ns+1;
no_bin=round((hi-lo)/dR);
P_ref=zeros(no_epoch,no_bin);
P_inf=zeros(no_epoch,no_bin);
Px_ref=zeros(no_epoch,2);
Px_inf=zeros(no_epoch,2);

batch_size=32;
no_batch=ceil(no_epoch/batch_size);

for ibatch=1:no_batch
    i_begin=(ibatch-1)*batch_size+1;
    i_end=min(ibatch*batch_size,no_epoch);
    epoch_no=i_begin:i_end;
    nb=length(epoch_no);
    t_target=(epoch_no-1)*Ts;

    a_ref_cell=cell(1,nb);
    e_ref_cell=cell(1,nb);
    n_ref=zeros(1,nb);

    for b=1:nb
        ii=epoch_no(b);
        n2=k2(ii);
        n3=k3(ii);
        n4=k4(ii);
        n_ref(b)=n2+n3+n4;

        a_ref_cell{b}=[reshape(ae2(ii,1:n2,1),1,[]),...
            reshape(ae3(ii,1:n3,1),1,[]),...
            reshape(ae4(ii,1:n4,1),1,[])];
        e_ref_cell{b}=[reshape(ae2(ii,1:n2,2),1,[]),...
            reshape(ae3(ii,1:n3,2),1,[]),...
            reshape(ae4(ii,1:n4,2),1,[])];
    end

    a_inf_batch=zeros(no_x,nb);
    e_inf_batch=zeros(no_x,nb);
    for j=1:no_x
        kj=k5_inf(j);
        if(kj<1)
            continue;
        end
        k_now=erfen_tae_sensitivity(tae_inf,j,kj,t_target(1));
        for b=1:nb
            tt=t_target(b);
            while(k_now<kj-1 && tae_inf(1,k_now+1,j)<=tt)
                k_now=k_now+1;
            end
            a_now=tae_inf(2,k_now,j);
            if(a_now>0)
                a_inf_batch(j,b)=a_now;
                e_inf_batch(j,b)=tae_inf(3,k_now,j);
            end
        end
    end

    P_ref_batch=zeros(nb,no_bin);
    P_inf_batch=zeros(nb,no_bin);
    Px_ref_batch=zeros(nb,2);
    Px_inf_batch=zeros(nb,2);

    if(use_thread_pool)
        parfor b=1:nb
            [pr,~,up_r,low_r]=kongjianfenbu_fast(lo,hi,dR,a_ref_cell{b},e_ref_cell{b});
            pr=pr/no_x;
            pxr=[up_r/no_x,(low_r+no_x-n_ref(b))/no_x];

            alive=(a_inf_batch(:,b)>0);
            n_alive=sum(alive);
            [pi,~,up_i,low_i]=kongjianfenbu_fast(lo,hi,dR,...
                a_inf_batch(alive,b).',e_inf_batch(alive,b).');
            pi=pi/no_x;
            pxi=[up_i/no_x,(low_i+no_x-n_alive)/no_x];

            P_ref_batch(b,:)=pr;
            P_inf_batch(b,:)=pi;
            Px_ref_batch(b,:)=pxr;
            Px_inf_batch(b,:)=pxi;
        end
    else
        for b=1:nb
            [pr,~,up_r,low_r]=kongjianfenbu_fast(lo,hi,dR,a_ref_cell{b},e_ref_cell{b});
            pr=pr/no_x;
            pxr=[up_r/no_x,(low_r+no_x-n_ref(b))/no_x];

            alive=(a_inf_batch(:,b)>0);
            n_alive=sum(alive);
            [pi,~,up_i,low_i]=kongjianfenbu_fast(lo,hi,dR,...
                a_inf_batch(alive,b).',e_inf_batch(alive,b).');
            pi=pi/no_x;
            pxi=[up_i/no_x,(low_i+no_x-n_alive)/no_x];

            P_ref_batch(b,:)=pr;
            P_inf_batch(b,:)=pi;
            Px_ref_batch(b,:)=pxr;
            Px_inf_batch(b,:)=pxi;
        end
    end

    P_ref(epoch_no,:)=P_ref_batch;
    P_inf(epoch_no,:)=P_inf_batch;
    Px_ref(epoch_no,:)=Px_ref_batch;
    Px_inf(epoch_no,:)=Px_inf_batch;

    if(ibatch==1 || ibatch==no_batch || mod(ibatch,max(1,floor(no_batch/20)))==0)
        fprintf('    radial distribution progress: %d/%d epochs (%.1f%%)\n',...
            i_end,no_epoch,100*i_end/no_epoch);
    end
end
end


%% ========================================================================
function P_coarse=rebin_probability(P_fine,factor)
% 将连续 factor 个 fine bins 相加。
[no_epoch,no_bin]=size(P_fine);
if(mod(no_bin,factor)~=0)
    error('Fine-bin number %d cannot be evenly grouped by factor %d.',no_bin,factor);
end
P_coarse=sum(reshape(P_fine.',factor,no_bin/factor,no_epoch),1);
P_coarse=permute(P_coarse,[3,2,1]);
P_coarse=reshape(P_coarse,no_epoch,no_bin/factor);
end


%% ========================================================================
function JS=one_epoch_JS(P_ref,Px_ref,P_inf,Px_inf)
p=[P_ref,Px_ref(1),Px_ref(2)];
q=[P_inf,Px_inf(1),Px_inf(2)];
m=0.5*(p+q);
idp=(p>0);
idq=(q>0);
JS=0.5*sum(p(idp).*log(p(idp)./m(idp)))+...
   0.5*sum(q(idq).*log(q(idq)./m(idq)));
end


%% ========================================================================
function [yspace,xspace,P_shang,P_xia]=kongjianfenbu_fast(lo,hi,dR,a,e)
% 与原 kongjianfenbu 使用同一解析表达式，但只遍历每颗碎片实际相交的 bins，
% 避免对所有 radial bins 逐一测试，适合 sensitivity 中的大量重复统计。

no_bin=round((hi-lo)/dR);
xspace=linspace(lo,hi-dR,no_bin)+dR/2;
yspace=zeros(1,no_bin);
P_shang=0;
P_xia=0;

if(isempty(a))
    xspace=xspace/1e3;
    return;
end

for i=1:length(a)
    ai=a(i);
    ei=e(i);
    if(ai<=0)
        continue;
    end

    if(abs(ei)<1e-14)
        if(ai<lo)
            P_xia=P_xia+1;
        elseif(ai>=hi)
            P_shang=P_shang+1;
        else
            j=floor((ai-lo)/dR)+1;
            j=max(1,min(no_bin,j));
            yspace(j)=yspace(j)+1;
        end
        continue;
    end

    rp=ai*(1-ei);
    ra=ai*(1+ei);

    if(ra>lo && rp<hi)
        j_begin=max(1,floor((rp-lo)/dR)+1);
        j_end=min(no_bin,ceil((ra-lo)/dR));

        if(j_end>=j_begin)
            jj=j_begin:j_end;
            r_down=lo+(jj-1)*dR;
            r_up=lo+jj*dR;

            jifen_down=asin(max((r_down-ai)/ai/ei,-1));
            jifen_shang=asin(min((r_up-ai)/ai/ei,1));
            yspace(jj)=yspace(jj)+...
                ((jifen_shang-jifen_down)-ei*(cos(jifen_shang)-cos(jifen_down)))/pi;
        end
    end

    if(ra>hi)
        jifen_down=asin(max((hi-ai)/ai/ei,-1));
        P_shang=P_shang+((pi/2-jifen_down)-ei*(cos(pi/2)-cos(jifen_down)))/pi;
    end

    if(rp<lo)
        jifen_shang=asin(min((lo-ai)/ai/ei,1));
        P_xia=P_xia+((jifen_shang+pi/2)-ei*(cos(jifen_shang)-cos(-pi/2)))/pi;
    end
end

xspace=xspace/1e3;
end


%% ========================================================================
function use_thread_pool=prepare_parallel_pool_sensitivity()
% 优先使用 thread-based pool，避免大数组在 process workers 之间复制。
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
    fprintf('Thread-based parallel pool started.\n');
catch ME
    warning('Thread-based pool is unavailable: %s',ME.message);
    fprintf(['Starting process-based pool for inference. ',...
        'Distribution/JS calculations will fall back to serial mode to avoid large-array broadcasting.\n']);
    parpool('local');
    use_thread_pool=false;
end
end


%% ========================================================================
function s=plural_s(n)
if(n==1)
    s='';
else
    s='s';
end
end
