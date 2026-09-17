%% Cosmos 2251: radial-distribution bands and JS-divergence band
% This script loads previously saved results (true_10_year.mat) and draws:
%   1) radial-distribution snapshots with uncertainty bands;
%   2) JS-divergence curves with uncertainty band.
%
% The plotting style follows plot_SC.m.

clear; clc;
set(groot, 'defaultAxesFontName', 'Times New Roman');

%% ========================= user settings =========================
G = 20;                         % delete-group jackknife groups
rng_seed = 20260917;            % fixed seed for reproducibility
reuse_band_cache = true;        % load cached jackknife results if available
save_cache = true;
show_smoothed_band_in_js = false; % default: only the original-JS curve shows a band

%% ========================= file locations =========================
this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir)
    this_dir = pwd;
end
addpath(this_dir);

if exist(fullfile(this_dir,'data'),'dir')
    data_dir = fullfile(this_dir,'data');
else
    data_dir = this_dir;
end

data_file = fullfile(data_dir, 'true_10_year.mat');
cache_file = fullfile(data_dir, 'true_10_year_jackknife_band_v2.mat');

if ~exist(data_file, 'file')
    error('Data file not found: %s', data_file);
end

%% ========================= load source data =========================
S = load(data_file, ...
    'DAY','big','small','real','tae_latter_tuiyan', ...
    'latter_real','latter_tuiyan','Px_real','Px_tuiyan', ...
    'JS_sandu','lo_fenbu','hi_fenbu','dR_fenbu');

DAY = S.DAY;
big = S.big;
small = S.small;
real_data = S.real;
tae_latter_tuiyan = S.tae_latter_tuiyan;
latter_real = S.latter_real;
latter_tuiyan = S.latter_tuiyan;
Px_real = S.Px_real;
Px_tuiyan = S.Px_tuiyan;
JS_sandu = S.JS_sandu;
lo_fenbu = S.lo_fenbu;
hi_fenbu = S.hi_fenbu;
dR_fenbu = S.dR_fenbu;

n_bins = size(latter_real, 2);
radius_km = latter_real(2,:,1);
snap_idx = [round(DAY/3), round(DAY/3*2), DAY];

%% ========================= smoothed nominal results =========================
JS_sandu_smooth = compute_js_smooth(latter_real, latter_tuiyan, Px_real, Px_tuiyan, lo_fenbu, hi_fenbu, dR_fenbu);

%% ========================= jackknife replicates =========================
if reuse_band_cache && exist(cache_file,'file')
    C = load(cache_file);
    if isfield(C,'group_id') && length(C.group_id)==small && isfield(C,'G') && C.G==G ...
            && isfield(C,'snap_idx') && isequal(C.snap_idx, snap_idx) ...
            && isfield(C,'JS_delete') && all(isfinite(C.JS_delete(:))) && isreal(C.JS_delete)
        fprintf('Loading cached jackknife bands from %s\n', cache_file);
        group_id = C.group_id;
        real_snap_delete = C.real_snap_delete;
        inf_snap_delete = C.inf_snap_delete;
        inf_snap_delete_smooth = C.inf_snap_delete_smooth;
        JS_delete = C.JS_delete;
        JS_delete_smooth = C.JS_delete_smooth;
    else
        warning('Cache exists but settings do not match. Recomputing jackknife bands.');
        [group_id, real_snap_delete, inf_snap_delete, inf_snap_delete_smooth, JS_delete, JS_delete_smooth] = ...
            compute_jackknife_delete_group(G, rng_seed, small, big, DAY, snap_idx, ...
            latter_real, latter_tuiyan, Px_real, Px_tuiyan, real_data, tae_latter_tuiyan, ...
            lo_fenbu, hi_fenbu, dR_fenbu);
        if save_cache
            save(cache_file, 'G', 'rng_seed', 'group_id', 'snap_idx', ...
                'real_snap_delete', 'inf_snap_delete', 'inf_snap_delete_smooth', ...
                'JS_delete', 'JS_delete_smooth', '-v7.3');
        end
    end
else
    [group_id, real_snap_delete, inf_snap_delete, inf_snap_delete_smooth, JS_delete, JS_delete_smooth] = ...
        compute_jackknife_delete_group(G, rng_seed, small, big, DAY, snap_idx, ...
        latter_real, latter_tuiyan, Px_real, Px_tuiyan, real_data, tae_latter_tuiyan, ...
        lo_fenbu, hi_fenbu, dR_fenbu);
    if save_cache
        save(cache_file, 'G', 'rng_seed', 'group_id', 'snap_idx', ...
            'real_snap_delete', 'inf_snap_delete', 'inf_snap_delete_smooth', ...
            'JS_delete', 'JS_delete_smooth', '-v7.3');
    end
end

%% ========================= derive uncertainty bands =========================
real_center = permute(latter_real(1,:,snap_idx), [3 2 1]);                % 3 x n_bins
inf_center = permute(latter_tuiyan(1,:,snap_idx), [3 2 1]);                % 3 x n_bins
inf_center_smooth = movmean(inf_center, 10, 2);

real_mean_delete = squeeze(mean(real_snap_delete, 1));
inf_mean_delete = squeeze(mean(inf_snap_delete, 1));
inf_smooth_mean_delete = squeeze(mean(inf_snap_delete_smooth, 1));

real_se = sqrt((G-1)/G * squeeze(sum((real_snap_delete - reshape(real_mean_delete,[1 size(real_mean_delete)])).^2, 1)));
inf_se = sqrt((G-1)/G * squeeze(sum((inf_snap_delete - reshape(inf_mean_delete,[1 size(inf_mean_delete)])).^2, 1)));
inf_smooth_se = sqrt((G-1)/G * squeeze(sum((inf_snap_delete_smooth - reshape(inf_smooth_mean_delete,[1 size(inf_smooth_mean_delete)])).^2, 1)));

real_low = max(0, real_center - 1.96 * real_se);
real_high = min(1, real_center + 1.96 * real_se);
inf_low = max(0, inf_center - 1.96 * inf_se);
inf_high = min(1, inf_center + 1.96 * inf_se);
inf_smooth_low = max(0, inf_center_smooth - 1.96 * inf_smooth_se);
inf_smooth_high = min(1, inf_center_smooth + 1.96 * inf_smooth_se);

JS_mean_delete = mean(JS_delete, 1);
JS_se = sqrt((G-1)/G * sum((JS_delete - JS_mean_delete).^2, 1));
JS_low = max(0, JS_sandu - 1.96 * JS_se);
JS_high = JS_sandu + 1.96 * JS_se;

JS_smooth_mean_delete = mean(JS_delete_smooth, 1);
JS_smooth_se = sqrt((G-1)/G * sum((JS_delete_smooth - JS_smooth_mean_delete).^2, 1));
JS_smooth_low = max(0, JS_sandu_smooth - 1.96 * JS_smooth_se);
JS_smooth_high = JS_sandu_smooth + 1.96 * JS_smooth_se;

%% ========================= Figure 1: radial distributions with bands =========================
fig_dist = figure('Color','w','Position',[50 50 1200 500],'WindowStyle','normal');
tl_dist = tiledlayout(fig_dist,2,3,'Padding','loose','TileSpacing','loose');
top_axes = gobjects(1,3);
bottom_axes = gobjects(1,3);

for k = 1:3
    pic_num = snap_idx(k);
    ax = nexttile(tl_dist,k);
    top_axes(k) = ax;
    hold(ax,'on');
    h_inf_band = fill(ax, [inf_low(k,:) fliplr(inf_high(k,:))], [radius_km fliplr(radius_km)], ...
        [0.78 0.18 0.68], 'EdgeColor','none', 'FaceAlpha',0.65);
    h_true_line = plot(ax, latter_real(1,:,pic_num), latter_real(2,:,pic_num), 'green', 'LineWidth', 1);
    h_inf_line = plot(ax, latter_tuiyan(1,:,pic_num), latter_tuiyan(2,:,pic_num), 'magenta', 'LineWidth', 1);
    hold(ax,'off');

    ylabel('Orbital radius (km)');
    xlabel('Probability per bin');
    legend(ax, [h_true_line, h_inf_line, h_inf_band], ...
        {'True distribution','Inferred distribution','Inferred 95% band'}, ...
        'Location','northeast');
    title(['Date: ', get_target_date(pic_num)],'FontSize',12,'FontWeight','normal');
    set(ax,'FontSize',11,'LineWidth',1);
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    axis([0 0.02 lo_fenbu/1e3 hi_fenbu/1e3]);
end

for k = 1:3
    pic_num = snap_idx(k);
    ax = nexttile(tl_dist,k+3);
    bottom_axes(k) = ax;
    hold(ax,'on');
    h_inf_s_band = fill(ax, [inf_smooth_low(k,:) fliplr(inf_smooth_high(k,:))], [radius_km fliplr(radius_km)], ...
        [0.78 0.18 0.68], 'EdgeColor','none', 'FaceAlpha',0.65);
    h_true_line = plot(ax, latter_real(1,:,pic_num), latter_real(2,:,pic_num), 'green', 'LineWidth', 1);
    h_inf_s_line = plot(ax, movmean(latter_tuiyan(1,:,pic_num),10), latter_tuiyan(2,:,pic_num), 'magenta', 'LineWidth', 1);
    hold(ax,'off');

    ylabel('Orbital radius (km)');
    xlabel('Probability per bin');
    legend(ax, [h_true_line, h_inf_s_line, h_inf_s_band], ...
        {'True distribution','Smoothed inferred distribution','Smoothed inferred 95% band'}, ...
        'Location','northeast');
    title(['Date: ', get_target_date(pic_num)],'FontSize',12,'FontWeight','normal');
    set(ax,'FontSize',12,'LineWidth',1);
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    axis([0 0.02 lo_fenbu/1e3 hi_fenbu/1e3]);
end

drawnow;
outfile_dist = fullfile(pwd,'fig_cosmos2251_radial_distribution_band.png');
exportgraphics(fig_dist, outfile_dist, 'Resolution', 600, 'BackgroundColor','white');

%% ========================= Figure 2: JS divergence with band =========================
% Two vertically stacked panels are used so that the uncertainty bands of
% the original and smoothed inferred distributions can be read separately.
fig_js = figure('Color','w','Position',[100 80 1000 700],'WindowStyle','normal');
tl_js = tiledlayout(fig_js,2,1,'Padding','compact','TileSpacing','compact');
t_days = 1:DAY;

% ----- (a) Original inferred distribution -----
ax_js1 = nexttile(tl_js,1);
hold(ax_js1,'on');
h_js_band = fill(ax_js1, [t_days fliplr(t_days)], [JS_low fliplr(JS_high)], ...
    [0.32 0.62 0.86], 'EdgeColor','none', 'FaceAlpha',0.62);
h_js_line = plot(ax_js1, t_days, JS_sandu, ...
    'Color', [0, 82, 155]/255, 'LineWidth', 1.0);
hold(ax_js1,'off');
ylabel(ax_js1,'JS divergence');
legend(ax_js1, [h_js_line, h_js_band], ...
    {'Original inferred','95% uncertainty band'}, ...
    'Location','southeast');
title(ax_js1,'(a) Original inferred distribution', ...
    'FontSize',14,'FontWeight','normal');
ylim(ax_js1,[0 0.08]);
set(ax_js1,'FontSize',15,'LineWidth',1,'XGrid','off','YGrid','on');

% ----- (b) Smoothed inferred distribution -----
ax_js2 = nexttile(tl_js,2);
hold(ax_js2,'on');
h_js_s_band = fill(ax_js2, [t_days fliplr(t_days)], [JS_smooth_low fliplr(JS_smooth_high)], ...
    [0.95 0.58 0.28], 'EdgeColor','none', 'FaceAlpha',0.55);
h_js_s_line = plot(ax_js2, t_days, JS_sandu_smooth, ...
    'Color', [230, 90, 13]/255, 'LineWidth', 1.0);
hold(ax_js2,'off');
xlabel(ax_js2,'Time elapsed since 2009-05-10 (days)');
ylabel(ax_js2,'JS divergence');
legend(ax_js2, [h_js_s_line, h_js_s_band], ...
    {'Smoothed inferred','95% uncertainty band'}, ...
    'Location','southeast');
title(ax_js2,'(b) Smoothed inferred distribution', ...
    'FontSize',14,'FontWeight','normal');
ylim(ax_js2,[0 0.08]);
set(ax_js2,'FontSize',15,'LineWidth',1,'XGrid','off','YGrid','on');

% Keep the same x range in both panels.
linkaxes([ax_js1,ax_js2],'x');
xlim(ax_js1,[0 4000]);

drawnow;
outfile_js = fullfile(pwd,'fig_cosmos2251_js_divergence_band_2panel.png');
exportgraphics(tl_js, outfile_js, 'Resolution', 600, 'BackgroundColor','white');

%% ========================= console output =========================
fprintf('Mean JS (original inferred) = %.8f\n', mean(JS_sandu));
fprintf('Mean JS (smoothed inferred) = %.8f\n', mean(JS_sandu_smooth));

fprintf('Mean JS-band half-width (original) = %.8e\n', mean(1.96 * JS_se));
fprintf('Max  JS-band half-width (original) = %.8e\n', max(1.96 * JS_se));
fprintf('Mean JS-band half-width (smoothed) = %.8e\n', mean(1.96 * JS_smooth_se));
fprintf('Max  JS-band half-width (smoothed) = %.8e\n', max(1.96 * JS_smooth_se));
fprintf('Non-finite JS jackknife values = %d\n', sum(~isfinite(JS_delete(:))));
fprintf('Complex JS jackknife flag     = %d\n', ~isreal(JS_delete));

%% ========================= local functions =========================
function JS_sandu_smooth = compute_js_smooth(latter_real, latter_tuiyan, Px_real, Px_tuiyan, lo_fenbu, hi_fenbu, dR_fenbu)
DAY = size(latter_real,3);
JS_sandu_smooth = zeros(1,DAY);
for i = 1:DAY
    P_KL = latter_real(1,:,i);
    Q_KL = movmean(latter_tuiyan(1,:,i),10);
    for j = 1:(hi_fenbu-lo_fenbu)/dR_fenbu
        if(P_KL(j)~=0)
            JS_sandu_smooth(i)=JS_sandu_smooth(i)+0.5*P_KL(j)*log(P_KL(j)*2/(P_KL(j)+Q_KL(j)));
        end
        if(Q_KL(j)~=0)
            JS_sandu_smooth(i)=JS_sandu_smooth(i)+0.5*Q_KL(j)*log(Q_KL(j)*2/(P_KL(j)+Q_KL(j)));
        end
    end
    if(Px_real(i,1)~=0)
        JS_sandu_smooth(i)=JS_sandu_smooth(i)+0.5*Px_real(i,1)*log(Px_real(i,1)*2/(Px_real(i,1)+Px_tuiyan(i,1)));
    end
    if(Px_real(i,2)~=0)
        JS_sandu_smooth(i)=JS_sandu_smooth(i)+0.5*Px_real(i,2)*log(Px_real(i,2)*2/(Px_real(i,2)+Px_tuiyan(i,2)));
    end
    if(Px_tuiyan(i,1)~=0)
        JS_sandu_smooth(i)=JS_sandu_smooth(i)+0.5*Px_tuiyan(i,1)*log(Px_tuiyan(i,1)*2/(Px_real(i,1)+Px_tuiyan(i,1)));
    end
    if(Px_tuiyan(i,2)~=0)
        JS_sandu_smooth(i)=JS_sandu_smooth(i)+0.5*Px_tuiyan(i,2)*log(Px_tuiyan(i,2)*2/(Px_real(i,2)+Px_tuiyan(i,2)));
    end
end
end

function [group_id, real_snap_delete, inf_snap_delete, inf_snap_delete_smooth, JS_delete, JS_delete_smooth] = ...
    compute_jackknife_delete_group(G, rng_seed, small, big, DAY, snap_idx, latter_real, latter_tuiyan, Px_real, Px_tuiyan, real_data, tae_latter_tuiyan, lo_fenbu, hi_fenbu, dR_fenbu)

rng(rng_seed, 'twister');
perm = randperm(small);
group_id = zeros(small,1);
base_group = repelem(1:G, ceil(small/G));
base_group = base_group(1:small);
group_id(perm) = base_group(:);

n_bins = size(latter_real, 2);
real_snap_delete = zeros(G, length(snap_idx), n_bins);
inf_snap_delete = zeros(G, length(snap_idx), n_bins);
inf_snap_delete_smooth = zeros(G, length(snap_idx), n_bins);
JS_delete = zeros(G, DAY);
JS_delete_smooth = zeros(G, DAY);

full_real_sum = squeeze(latter_real(1,:,:)) * small;   % n_bins x DAY
full_inf_sum = squeeze(latter_tuiyan(1,:,:)) * small;   % n_bins x DAY
full_real_up = Px_real(:,1).' * small;
full_real_low = Px_real(:,2).' * small;
full_inf_up = Px_tuiyan(:,1).' * small;
full_inf_low = Px_tuiyan(:,2).' * small;

for g = 1:G
    ids_grp = find(group_id == g);
    n_grp = numel(ids_grp);
    n_remain = small - n_grp;
    fprintf('Jackknife group %d/%d, removed target fragments = %d\n', g, G, n_grp);

    snap_counter = 1;
    for it = 1:DAY
        % True population: use the latter part of real(:,:,big+1:big+small).
        a_true = reshape(real_data(1,it,big+ids_grp), 1, []);
        e_true = reshape(real_data(2,it,big+ids_grp), 1, []);
        [real_bins_removed, ~, real_up_removed, real_low_removed] = ...
            spaceDistribution_local(lo_fenbu, hi_fenbu, dR_fenbu, a_true, e_true);

        % Inferred population.
        a_inf = reshape(tae_latter_tuiyan(2,it,ids_grp), 1, []);
        e_inf = reshape(tae_latter_tuiyan(3,it,ids_grp), 1, []);
        [inf_bins_removed, ~, inf_up_removed, inf_low_removed] = ...
            spaceDistribution_local(lo_fenbu, hi_fenbu, dR_fenbu, a_inf, e_inf);

        % Delete-group distributions.
        P_real_del = (full_real_sum(:,it).' - real_bins_removed) / n_remain;
        P_inf_del = (full_inf_sum(:,it).' - inf_bins_removed) / n_remain;

        Px_real_up_del = (full_real_up(it) - real_up_removed) / n_remain;
        Px_real_low_del = (full_real_low(it) - real_low_removed) / n_remain;
        Px_inf_up_del = (full_inf_up(it) - inf_up_removed) / n_remain;
        Px_inf_low_del = (full_inf_low(it) - inf_low_removed) / n_remain;

        % Independent subtraction can leave tiny negative round-off values.
        % These lead to complex/NaN values in log() during JS evaluation.
        [P_real_del, Px_real_low_del, Px_real_up_del] = ...
            sanitize_probability_components(P_real_del, Px_real_low_del, Px_real_up_del);
        [P_inf_del, Px_inf_low_del, Px_inf_up_del] = ...
            sanitize_probability_components(P_inf_del, Px_inf_low_del, Px_inf_up_del);

        P_inf_del_smooth = movmean(P_inf_del, 10);

        JS_delete(g,it) = one_epoch_js(P_real_del, P_inf_del, Px_real_up_del, Px_real_low_del, Px_inf_up_del, Px_inf_low_del);
        JS_delete_smooth(g,it) = one_epoch_js(P_real_del, P_inf_del_smooth, Px_real_up_del, Px_real_low_del, Px_inf_up_del, Px_inf_low_del);

        if snap_counter <= length(snap_idx) && it == snap_idx(snap_counter)
            real_snap_delete(g, snap_counter, :) = P_real_del;
            inf_snap_delete(g, snap_counter, :) = P_inf_del;
            inf_snap_delete_smooth(g, snap_counter, :) = P_inf_del_smooth;
            snap_counter = snap_counter + 1;
        end
    end
end
end

function JS_one = one_epoch_js(P_ref, P_inf, P_ref_up, P_ref_low, P_inf_up, P_inf_low)
% Robust JS-divergence evaluation.
P = real([P_ref_low, P_ref, P_ref_up]);
Q = real([P_inf_low, P_inf, P_inf_up]);

P(~isfinite(P)) = 0;
Q(~isfinite(Q)) = 0;
P(P < 0) = 0;
Q(Q < 0) = 0;

M = 0.5 * (P + Q);
idxP = (P > 0) & (M > 0);
idxQ = (Q > 0) & (M > 0);

JS_one = 0.5 * sum(P(idxP) .* log(P(idxP) ./ M(idxP))) ...
       + 0.5 * sum(Q(idxQ) .* log(Q(idxQ) ./ M(idxQ)));
JS_one = real(JS_one);

if ~isfinite(JS_one)
    error('Non-finite JS divergence encountered after probability sanitization.');
end
end

function [P_bins, P_low, P_up] = sanitize_probability_components(P_bins, P_low, P_up)
% Internal bins + lower/upper boundary probabilities should sum to one.
% Remove floating-point residuals and renormalize.
P_bins = real(P_bins);
P_low = real(P_low);
P_up = real(P_up);

P_bins(~isfinite(P_bins)) = 0;
if ~isfinite(P_low), P_low = 0; end
if ~isfinite(P_up),  P_up = 0;  end

P_bins(P_bins < 0) = 0;
P_low = max(P_low, 0);
P_up = max(P_up, 0);

s = sum(P_bins) + P_low + P_up;
if s <= 0
    error('Probability normalization failed: total probability is non-positive.');
end

P_bins = P_bins / s;
P_low = P_low / s;
P_up = P_up / s;
end

function [yspace, xspace, P_shang, P_xia] = spaceDistribution_local(lo,hi,dR,a,e)
xspace=linspace(lo,hi-dR,(hi-lo)/dR)+dR/2;
yspace=zeros(1,(hi-lo)/dR);
[~,lie_a]=size(a);
[~,lie_e]=size(e);
if(lie_a~=lie_e)
    fprintf("error!spaceDistribution");
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
xspace=xspace/1e3;
end
