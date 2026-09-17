%% Closed numerical experiment: radial-distribution bands and JS-divergence band
% This script loads the saved .mat file and generates:
% 1) three snapshot radial-distribution plots with uncertainty bands for both
%    the reference and inferred populations;
% 2) one baseline JS-divergence plot with an uncertainty band.
%
% Uncertainty quantification:
% delete-group jackknife with G groups.
%
% Plotting style follows plot_SC_4.m.

clear; clc;
set(groot, 'defaultAxesFontName', 'Times New Roman');

%% ========================= user settings =========================
zhongzi = 2;
ns_in_name = 21000;
G = 20;                         % delete-group jackknife groups
rng_seed = 20260917;            % fixed seed for reproducibility
reuse_band_cache = true;        % if cache exists, load it directly
save_cache = true;

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

data_file = fullfile(data_dir, sprintf('zhongzi=%d_ns=%d.mat', zhongzi, ns_in_name));
cache_file = fullfile(data_dir, sprintf('zhongzi=%d_ns=%d_jackknife_band_v2.mat', zhongzi, ns_in_name));

if ~exist(data_file,'file')
    error('Data file not found: %s', data_file);
end

%% ========================= load source data =========================
S = load(data_file, ...
    'xiao_cankao','xiao_tuiyan','Px_cankao','Px_tuiyan', ...
    'tae_xiao_tuiyan','k5', ...
    'yxiao_cankao','hang_yxiao_cankao', ...
    'Ts','ns','no_x','lo_fenbu','hi_fenbu','dR_fenbu','miu_earth');

xiao_cankao = S.xiao_cankao;
xiao_tuiyan = S.xiao_tuiyan;
Px_cankao = S.Px_cankao;
Px_tuiyan = S.Px_tuiyan;
tae_xiao_tuiyan = S.tae_xiao_tuiyan;
k5 = S.k5;
yxiao_cankao = S.yxiao_cankao;
hang_yxiao_cankao = S.hang_yxiao_cankao;
Ts = S.Ts;
ns = S.ns;
no_x = S.no_x;
lo_fenbu = S.lo_fenbu;
hi_fenbu = S.hi_fenbu;
dR_fenbu = S.dR_fenbu;
miu_earth = S.miu_earth;

n_bins = size(xiao_cankao, 2);
radius_km = xiao_cankao(2,:,1);
t_days = (0:ns) * Ts / 86400;

snap_idx = 1:ns/3:ns+1;
snap_idx = round(snap_idx(2:4));

%% ========================= full-sample JS =========================
if isfield(S, 'JS_sandu')
    JS_sandu = S.JS_sandu;
else
    JS_sandu = compute_js_from_distribution(xiao_cankao, xiao_tuiyan, Px_cankao, Px_tuiyan);
end

%% ========================= jackknife bands =========================
if reuse_band_cache && exist(cache_file,'file')
    C = load(cache_file);
    if isfield(C,'group_id') && length(C.group_id)==no_x && isfield(C,'G') && C.G==G ...
            && isfield(C,'snap_idx') && isequal(C.snap_idx, snap_idx) ...
            && isfield(C,'JS_delete') && all(isfinite(C.JS_delete(:))) && isreal(C.JS_delete)
        fprintf('Loading cached jackknife bands from %s\n', cache_file);
        group_id = C.group_id;
        ref_snap_delete = C.ref_snap_delete;
        inf_snap_delete = C.inf_snap_delete;
        JS_delete = C.JS_delete;
    else
        warning('Cache exists but settings do not match. Recomputing jackknife bands.');
        [group_id, ref_snap_delete, inf_snap_delete, JS_delete] = ...
            compute_jackknife_delete_group(G, rng_seed, no_x, snap_idx, Ts, ...
            xiao_cankao, xiao_tuiyan, Px_cankao, Px_tuiyan, ...
            yxiao_cankao, hang_yxiao_cankao, miu_earth, ...
            tae_xiao_tuiyan, k5, lo_fenbu, hi_fenbu, dR_fenbu);
        if save_cache
            save(cache_file, 'G', 'rng_seed', 'group_id', 'snap_idx', ...
                'ref_snap_delete', 'inf_snap_delete', 'JS_delete', '-v7.3');
        end
    end
else
    [group_id, ref_snap_delete, inf_snap_delete, JS_delete] = ...
        compute_jackknife_delete_group(G, rng_seed, no_x, snap_idx, Ts, ...
        xiao_cankao, xiao_tuiyan, Px_cankao, Px_tuiyan, ...
        yxiao_cankao, hang_yxiao_cankao, miu_earth, ...
        tae_xiao_tuiyan, k5, lo_fenbu, hi_fenbu, dR_fenbu);
    if save_cache
        save(cache_file, 'G', 'rng_seed', 'group_id', 'snap_idx', ...
            'ref_snap_delete', 'inf_snap_delete', 'JS_delete', '-v7.3');
    end
end

%% ========================= derive uncertainty bands =========================
% Distribution bands (reference and inferred) at the three snapshots.
ref_center = permute(xiao_cankao(1,:,snap_idx), [3 2 1]);   % 3 x n_bins
inf_center = permute(xiao_tuiyan(1,:,snap_idx), [3 2 1]);   % 3 x n_bins

ref_mean_delete = squeeze(mean(ref_snap_delete, 1));        % 3 x n_bins
inf_mean_delete = squeeze(mean(inf_snap_delete, 1));        % 3 x n_bins

ref_se = sqrt((G-1)/G * squeeze(sum((ref_snap_delete - reshape(ref_mean_delete,[1 size(ref_mean_delete)])).^2, 1)));
inf_se = sqrt((G-1)/G * squeeze(sum((inf_snap_delete - reshape(inf_mean_delete,[1 size(inf_mean_delete)])).^2, 1)));

ref_low = max(0, ref_center - 1.96 * ref_se);
ref_high = min(1, ref_center + 1.96 * ref_se);
inf_low = max(0, inf_center - 1.96 * inf_se);
inf_high = min(1, inf_center + 1.96 * inf_se);

% JS band.
JS_mean_delete = mean(JS_delete, 1);
JS_se = sqrt((G-1)/G * sum((JS_delete - JS_mean_delete).^2, 1));
JS_low = max(0, JS_sandu - 1.96 * JS_se);
JS_high = JS_sandu + 1.96 * JS_se;

%% ========================= Figure 1: radial distributions with bands =========================
figure('Color','w','Position',[50 50 1200 420],'WindowStyle','normal');
tl_dist = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

for k = 1:3
    pic_num = snap_idx(k);
    ax = nexttile(tl_dist, k);
    hold(ax,'on');

    y_ref = radius_km;
    h_inf_band = fill(ax, [inf_low(k,:) fliplr(inf_high(k,:))], [y_ref fliplr(y_ref)], ...
        [0.78 0.18 0.68], 'EdgeColor','none', 'FaceAlpha',0.65);

    h_ref_line = plot(ax, xiao_cankao(1,:,pic_num), xiao_cankao(2,:,pic_num), ...
        'Color','green','LineWidth',1);
    h_inf_line = plot(ax, xiao_tuiyan(1,:,pic_num), xiao_tuiyan(2,:,pic_num), ...
        'Color','magenta','LineWidth',1);
    hold(ax,'off');

    ylabel('Orbital radius (km)');
    xlabel('Probability per bin');
    title(['$t = $', num2str(round((pic_num-1)*Ts/86400)), ' days'], ...
        'Interpreter','latex','FontSize',14);
    legend(ax, [h_ref_line, h_inf_line, h_inf_band], ...
        {'Reference distribution','Inferred distribution','Inferred 95% band'}, ...
        'Location','northeast');

    set(ax,'FontSize',12,'LineWidth',1,'XGrid','on','YGrid','on');
    axis(ax,[0 0.005 lo_fenbu/1e3 hi_fenbu/1e3]);
end

%% ========================= Figure 2: baseline JS with band =========================
figure('Color','w','Position',[100 100 720 450],'WindowStyle','normal');
tl_js = tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
ax4 = nexttile(tl_js,1);
h_js_band = fill(ax4, [t_days fliplr(t_days)], [JS_low fliplr(JS_high)], ...
    [0.32 0.62 0.86], 'EdgeColor','none', 'FaceAlpha',0.58); hold(ax4,'on');
h_js_line = plot(ax4, t_days, JS_sandu, 'Color', [0, 82, 155]/255, 'LineWidth', 1); hold(ax4,'off');
xlabel('Time (days)');
ylabel('JS divergence');
ylim([0 0.05]);
legend(ax4, [h_js_line, h_js_band], {'Baseline JS divergence','95% uncertainty band'}, 'Location','northeast');
set(ax4,'FontSize',15,'LineWidth',1,'XGrid','off','YGrid','on');

%% ========================= console output =========================
fprintf('Mean baseline JS divergence = %.8f\n', mean(JS_sandu));
fprintf('Max baseline JS divergence  = %.8f\n', max(JS_sandu));
fprintf('Mean JS-band half-width     = %.8e\n', mean(1.96 * JS_se));
fprintf('Max  JS-band half-width     = %.8e\n', max(1.96 * JS_se));
fprintf('Non-finite JS jackknife values = %d\n', sum(~isfinite(JS_delete(:))));
fprintf('Complex JS jackknife flag      = %d\n', ~isreal(JS_delete));

%% ========================= export =========================
exportgraphics(tl_dist, fullfile(this_dir,'fig_closed_distribution_band.png'), ...
    'Resolution',600,'BackgroundColor','white');
exportgraphics(tl_js, fullfile(this_dir,'fig_closed_js_band.png'), ...
    'Resolution',600,'BackgroundColor','white');

%% ========================= local functions =========================
function JS_sandu = compute_js_from_distribution(xiao_cankao, xiao_tuiyan, Px_cankao, Px_tuiyan)
ns = size(xiao_cankao, 3) - 1;
JS_sandu = zeros(1, ns+1);
for i = 1:ns+1
    P_KL = xiao_cankao(1,:,i);
    Q_KL = xiao_tuiyan(1,:,i);
    for j = 1:length(P_KL)
        if P_KL(j) ~= 0
            JS_sandu(i) = JS_sandu(i) + 0.5 * P_KL(j) * log(P_KL(j) * 2 / (P_KL(j) + Q_KL(j)));
        end
        if Q_KL(j) ~= 0
            JS_sandu(i) = JS_sandu(i) + 0.5 * Q_KL(j) * log(Q_KL(j) * 2 / (P_KL(j) + Q_KL(j)));
        end
    end
    if Px_cankao(i,1) ~= 0
        JS_sandu(i) = JS_sandu(i) + 0.5 * Px_cankao(i,1) * log(Px_cankao(i,1) * 2 / (Px_cankao(i,1) + Px_tuiyan(i,1)));
    end
    if Px_cankao(i,2) ~= 0
        JS_sandu(i) = JS_sandu(i) + 0.5 * Px_cankao(i,2) * log(Px_cankao(i,2) * 2 / (Px_cankao(i,2) + Px_tuiyan(i,2)));
    end
    if Px_tuiyan(i,1) ~= 0
        JS_sandu(i) = JS_sandu(i) + 0.5 * Px_tuiyan(i,1) * log(Px_tuiyan(i,1) * 2 / (Px_cankao(i,1) + Px_tuiyan(i,1)));
    end
    if Px_tuiyan(i,2) ~= 0
        JS_sandu(i) = JS_sandu(i) + 0.5 * Px_tuiyan(i,2) * log(Px_tuiyan(i,2) * 2 / (Px_cankao(i,2) + Px_tuiyan(i,2)));
    end
end
end

function [group_id, ref_snap_delete, inf_snap_delete, JS_delete] = compute_jackknife_delete_group(...
    G, rng_seed, no_x, snap_idx, Ts, xiao_cankao, xiao_tuiyan, Px_cankao, Px_tuiyan, ...
    yxiao_cankao, hang_yxiao_cankao, miu_earth, tae_xiao_tuiyan, k5, lo_fenbu, hi_fenbu, dR_fenbu)

rng(rng_seed, 'twister');
perm = randperm(no_x);
group_id = zeros(no_x,1);
base_group = repelem(1:G, ceil(no_x/G));
base_group = base_group(1:no_x);
group_id(perm) = base_group(:);

n_bins = size(xiao_cankao, 2);
ns = size(xiao_cankao, 3) - 1;

ref_snap_delete = zeros(G, length(snap_idx), n_bins);
inf_snap_delete = zeros(G, length(snap_idx), n_bins);
JS_delete = zeros(G, ns+1);

full_ref_sum = squeeze(xiao_cankao(1,:,:)) * no_x;   % n_bins x (ns+1)
full_inf_sum = squeeze(xiao_tuiyan(1,:,:)) * no_x;   % n_bins x (ns+1)
full_ref_up = Px_cankao(:,1).' * no_x;
full_ref_low = Px_cankao(:,2).' * no_x;
full_inf_up = Px_tuiyan(:,1).' * no_x;
full_inf_low = Px_tuiyan(:,2).' * no_x;

for g = 1:G
    ids_grp = find(group_id == g);
    n_grp = numel(ids_grp);
    n_remain = no_x - n_grp;
    fprintf('Jackknife group %d/%d, removed fragments = %d\n', g, G, n_grp);

    % pointer for inferred fragment histories in the removed group
    ptr = ones(1, n_grp);
    snap_counter = 1;

    for it = 1:ns+1
        t_target = (it-1) * Ts;

        % removed-group reference contributions at the current epoch
        [ref_bins_removed, ref_up_removed, ref_low_removed] = ...
            one_epoch_removed_reference(it, ids_grp, yxiao_cankao, hang_yxiao_cankao, miu_earth, lo_fenbu, hi_fenbu, dR_fenbu);

        % removed-group inferred contributions at the current epoch
        [inf_bins_removed, inf_up_removed, inf_low_removed, ptr] = ...
            one_epoch_removed_inferred(t_target, ids_grp, ptr, tae_xiao_tuiyan, k5, lo_fenbu, hi_fenbu, dR_fenbu);

        % delete-group distributions
        P_ref_del = (full_ref_sum(:,it).' - ref_bins_removed) / n_remain;
        P_inf_del = (full_inf_sum(:,it).' - inf_bins_removed) / n_remain;
        Px_ref_up_del = (full_ref_up(it) - ref_up_removed) / n_remain;
        Px_ref_low_del = (full_ref_low(it) - ref_low_removed) / n_remain;
        Px_inf_up_del = (full_inf_up(it) - inf_up_removed) / n_remain;
        Px_inf_low_del = (full_inf_low(it) - inf_low_removed) / n_remain;

        % The full-population sum and removed-group contribution are
        % evaluated independently. Their subtraction can leave tiny negative
        % round-off residuals (e.g. -1e-16), which would make log() complex
        % in the JS calculation. Clip those numerical residues and
        % renormalize the complete probability vectors before computing JS.
        [P_ref_del, Px_ref_low_del, Px_ref_up_del] = ...
            sanitize_probability_components(P_ref_del, Px_ref_low_del, Px_ref_up_del);
        [P_inf_del, Px_inf_low_del, Px_inf_up_del] = ...
            sanitize_probability_components(P_inf_del, Px_inf_low_del, Px_inf_up_del);

        % JS for this delete-group replicate
        JS_delete(g,it) = one_epoch_js(P_ref_del, P_inf_del, ...
            Px_ref_up_del, Px_ref_low_del, Px_inf_up_del, Px_inf_low_del);

        % save distributions at the snapshot epochs
        if snap_counter <= length(snap_idx) && it == snap_idx(snap_counter)
            ref_snap_delete(g, snap_counter, :) = P_ref_del;
            inf_snap_delete(g, snap_counter, :) = P_inf_del;
            snap_counter = snap_counter + 1;
        end
    end
end
end

function [bins_removed, P_up_removed, P_low_removed] = one_epoch_removed_reference(...
    eval_idx, ids_grp, yxiao_cankao, hang_yxiao_cankao, miu_earth, lo_fenbu, hi_fenbu, dR_fenbu)

n_grp = numel(ids_grp);
a = zeros(1, n_grp);
e = zeros(1, n_grp);
count = 0;
n_reentered_removed = 0;

for ii = 1:n_grp
    idx = ids_grp(ii);
    if eval_idx <= hang_yxiao_cankao(idx)
        r = squeeze(yxiao_cankao(eval_idx,1:3,idx)).';
        v = squeeze(yxiao_cankao(eval_idx,4:6,idx)).';
        [a_now, e_now, ~] = r0v0_genshu(r, v, miu_earth);
        if a_now > 0
            count = count + 1;
            a(count) = a_now;
            e(count) = e_now;
        else
            n_reentered_removed = n_reentered_removed + 1;
        end
    else
        % In the full-population reconstruction, a reentered fragment
        % contributes unit probability to the lower boundary component.
        n_reentered_removed = n_reentered_removed + 1;
    end
end

[bins_removed, ~, P_up_removed, P_low_removed] = kongjianfenbu_local(lo_fenbu, hi_fenbu, dR_fenbu, a(1:count), e(1:count));
P_low_removed = P_low_removed + n_reentered_removed;
end

function [bins_removed, P_up_removed, P_low_removed, ptr] = one_epoch_removed_inferred(...
    t_target, ids_grp, ptr, tae_xiao_tuiyan, k5, lo_fenbu, hi_fenbu, dR_fenbu)

n_grp = numel(ids_grp);
a = zeros(1, n_grp);
e = zeros(1, n_grp);
count = 0;
n_reentered_removed = 0;

for ii = 1:n_grp
    idx = ids_grp(ii);
    kj = k5(idx);
    while ptr(ii) < kj && tae_xiao_tuiyan(1, ptr(ii)+1, idx) <= t_target
        ptr(ii) = ptr(ii) + 1;
    end
    a_now = tae_xiao_tuiyan(2, ptr(ii), idx);
    if a_now > 0
        count = count + 1;
        a(count) = a_now;
        e(count) = tae_xiao_tuiyan(3, ptr(ii), idx);
    else
        % Same convention as the nominal inferred distribution: once the
        % fragment has reentered, it contributes unit probability below the
        % reconstruction domain.
        n_reentered_removed = n_reentered_removed + 1;
    end
end

[bins_removed, ~, P_up_removed, P_low_removed] = kongjianfenbu_local(lo_fenbu, hi_fenbu, dR_fenbu, a(1:count), e(1:count));
P_low_removed = P_low_removed + n_reentered_removed;
end

function JS_one = one_epoch_js(P_ref, P_inf, P_ref_up, P_ref_low, P_inf_up, P_inf_low)
% Robust JS-divergence evaluation. Tiny negative values caused by floating-
% point subtraction are clipped to zero before taking logarithms.
P = real([P_ref_low, P_ref, P_ref_up]);
Q = real([P_inf_low, P_inf, P_inf_up]);

P(~isfinite(P)) = 0;
Q(~isfinite(Q)) = 0;
P(P < 0) = 0;
Q(Q < 0) = 0;

M = 0.5 * (P + Q);
JS_one = 0;

idxP = (P > 0) & (M > 0);
idxQ = (Q > 0) & (M > 0);
JS_one = JS_one + 0.5 * sum(P(idxP) .* log(P(idxP) ./ M(idxP)));
JS_one = JS_one + 0.5 * sum(Q(idxQ) .* log(Q(idxQ) ./ M(idxQ)));

JS_one = real(JS_one);
if ~isfinite(JS_one)
    error('Non-finite JS divergence encountered after probability sanitization.');
end
end

function [P_bins, P_low, P_up] = sanitize_probability_components(P_bins, P_low, P_up)
% Internal bins plus lower/upper boundary components should form a complete
% probability distribution. Remove tiny round-off residuals and renormalize.
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

function [yspace, xspace, P_shang, P_xia] = kongjianfenbu_local(lo, hi, dR, a, e)
% Same logic as the original kongjianfenbu in zhuchengxu.m.
xspace = linspace(lo, hi-dR, (hi-lo)/dR) + dR/2;
yspace = zeros(1, (hi-lo)/dR);
P_shang = 0;
P_xia = 0;

for i = 1:length(a)
    for j = 1:(hi-lo)/dR
        r_down = lo + (j-1) * dR;
        r_up = lo + j * dR;
        if a(i) * (1 + e(i)) > r_down && a(i) * (1 - e(i)) < r_up
            jifen_down = asin(max((r_down - a(i)) / a(i) / e(i), -1));
            jifen_shang = asin(min((r_up - a(i)) / a(i) / e(i), 1));
            yspace(j) = yspace(j) + ((jifen_shang - jifen_down) - e(i) * (cos(jifen_shang) - cos(jifen_down))) / pi;
        end
    end
    if a(i) * (1 + e(i)) > hi
        jifen_down = asin(max((hi - a(i)) / a(i) / e(i), -1));
        jifen_shang = asin(1);
        P_shang = P_shang + ((jifen_shang - jifen_down) - e(i) * (cos(jifen_shang) - cos(jifen_down))) / pi;
    end
    if a(i) * (1 - e(i)) < lo
        jifen_down = asin(-1);
        jifen_shang = asin(min((lo - a(i)) / a(i) / e(i), 1));
        P_xia = P_xia + ((jifen_shang - jifen_down) - e(i) * (cos(jifen_shang) - cos(jifen_down))) / pi;
    end
end

xspace = xspace / 1e3;
end
