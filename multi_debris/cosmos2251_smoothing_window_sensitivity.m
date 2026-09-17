%% Cosmos 2251: smoothing-window sensitivity analysis
% Uses the saved population-level inference results only; no O2UDI rerun is needed.
% The calculation follows the original plot_SC.m convention:
%   - smooth ONLY the inferred internal radial-bin probabilities;
%   - leave the true distribution unchanged;
%   - leave the two boundary probability components unchanged;
%   - evaluate Jensen-Shannon divergence at every daily epoch.

clear; clc;

%% ========================= user settings =========================
data_file = 'true_10_year.mat';

% Five smoothing-window settings (radial bins).
% With the baseline 1-km radial bin width, these correspond to
% approximately 5, 10, 15, 20, and 30 km.
window_bins = [5, 10, 15, 20, 30];
baseline_window = 10;

% Save the statistical summary as CSV.
save_csv = true;
out_csv = 'cosmos2251_smoothing_window_sensitivity.csv';

%% ========================= load saved results =========================
if ~exist(data_file, 'file')
    error('Data file not found: %s', data_file);
end

S = load(data_file);

required_vars = {'latter_real','latter_tuiyan','Px_real','Px_tuiyan'};
for k = 1:numel(required_vars)
    if ~isfield(S, required_vars{k})
        error('Required variable "%s" is missing from %s.', required_vars{k}, data_file);
    end
end

latter_real   = S.latter_real;
latter_tuiyan = S.latter_tuiyan;
Px_real       = S.Px_real;
Px_tuiyan     = S.Px_tuiyan;

DAY = size(latter_real, 3);
if size(latter_tuiyan, 3) ~= DAY
    error('latter_real and latter_tuiyan have different numbers of epochs.');
end
if size(Px_real,1) ~= DAY || size(Px_tuiyan,1) ~= DAY
    error('Boundary-probability arrays are inconsistent with the number of epochs.');
end

%% ========================= original unsmoothed JS =========================
% Prefer the saved original JS series when available, because this is the
% exact baseline already used in the manuscript. Otherwise recompute it.
if isfield(S, 'JS_sandu') && numel(S.JS_sandu) == DAY
    JS_original = reshape(S.JS_sandu, 1, []);
else
    JS_original = compute_js_series( ...
        latter_real, latter_tuiyan, Px_real, Px_tuiyan, 1);
end

mean_original = mean(JS_original);
max_original  = max(JS_original);

%% ========================= smoothing-window sensitivity =========================
n_case = numel(window_bins);
JS_all = zeros(n_case, DAY);
mean_JS = zeros(n_case,1);
max_JS  = zeros(n_case,1);
mean_reduction = zeros(n_case,1);

for icase = 1:n_case
    w = window_bins(icase);

    JS_now = compute_js_series( ...
        latter_real, latter_tuiyan, Px_real, Px_tuiyan, w);

    JS_all(icase,:) = JS_now;
    mean_JS(icase) = mean(JS_now);
    max_JS(icase) = max(JS_now);
    mean_reduction(icase) = mean_original - mean_JS(icase);
end

%% ========================= comparison with 10-bin baseline =========================
idx_base = find(window_bins == baseline_window, 1);
if isempty(idx_base)
    error('baseline_window = %d is not included in window_bins.', baseline_window);
end

mean_base = mean_JS(idx_base);
max_base  = max_JS(idx_base);

delta_mean_vs_10 = mean_JS - mean_base;
delta_max_vs_10  = max_JS - max_base;

% Physical smoothing width, if the saved result contains the radial-bin width.
if isfield(S, 'dR_fenbu')
    % Original code uses dR_fenbu in metres.
    window_km = window_bins(:) * S.dR_fenbu / 1e3;
else
    window_km = nan(n_case,1);
end

%% ========================= output table =========================
Result = table( ...
    window_bins(:), window_km, mean_JS, max_JS, mean_reduction, ...
    delta_mean_vs_10, delta_max_vs_10, ...
    'VariableNames', { ...
    'Window_bins', 'Window_km', 'Mean_JS', 'Max_JS', ...
    'Mean_reduction_vs_original', 'Delta_mean_vs_10bin', 'Delta_max_vs_10bin'});

fprintf('\n===============================================================\n');
fprintf('Cosmos 2251 smoothing-window sensitivity\n');
fprintf('===============================================================\n');
fprintf('Original unsmoothed: mean JS = %.8f, max JS = %.8f\n\n', ...
    mean_original, max_original);
disp(Result);

fprintf('10-bin baseline:      mean JS = %.8f, max JS = %.8f\n', ...
    mean_base, max_base);
fprintf('Sensitivity range:    mean JS = %.8f -- %.8f\n', ...
    min(mean_JS), max(mean_JS));
fprintf('Maximum |change| from 10-bin baseline in mean JS = %.8f\n', ...
    max(abs(delta_mean_vs_10)));
fprintf('Maximum |change| from 10-bin baseline in max JS  = %.8f\n', ...
    max(abs(delta_max_vs_10)));

if save_csv
    writetable(Result, out_csv);
    fprintf('\nSaved summary to: %s\n', out_csv);
end

% Also save the complete JS histories in case they are needed later.
save('cosmos2251_smoothing_window_sensitivity.mat', ...
    'window_bins', 'window_km', 'JS_all', 'JS_original', 'Result');

%% ========================= local function =========================
function JS = compute_js_series(latter_real, latter_tuiyan, Px_real, Px_tuiyan, window_size)
% Compute the Jensen-Shannon divergence using the same probability
% components as the original Cosmos 2251 analysis.
%
% Internal bins:
%   P = true distribution (unchanged)
%   Q = inferred distribution after movmean smoothing
%
% Boundary components:
%   Px_real and Px_tuiyan are retained unchanged, consistent with the
%   original plot_SC.m calculation.

DAY = size(latter_real, 3);
JS = zeros(1, DAY);

for i = 1:DAY
    P_bins = reshape(latter_real(1,:,i), 1, []);
    Q_bins = reshape(latter_tuiyan(1,:,i), 1, []);

    if window_size > 1
        Q_bins = movmean(Q_bins, window_size);
    end

    % Complete discrete distributions: upper boundary + internal bins + lower boundary.
    % The order of the two boundary components does not affect JS divergence as
    % long as P and Q use the same order.
    P = [Px_real(i,1), P_bins, Px_real(i,2)];
    Q = [Px_tuiyan(i,1), Q_bins, Px_tuiyan(i,2)];

    % Remove tiny round-off negatives only; do not otherwise alter the saved data.
    P(P < 0 & P > -1e-14) = 0;
    Q(Q < 0 & Q > -1e-14) = 0;

    if any(P < 0) || any(Q < 0)
        error('Negative probability encountered at epoch %d.', i);
    end

    M = 0.5 * (P + Q);

    idxP = P > 0;
    idxQ = Q > 0;

    JS(i) = 0.5 * sum(P(idxP) .* log(P(idxP) ./ M(idxP))) + ...
            0.5 * sum(Q(idxQ) .* log(Q(idxQ) ./ M(idxQ)));
end
end
