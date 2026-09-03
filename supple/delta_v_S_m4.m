% Combined distributions of ejection-velocity increment and area-to-mass ratio
% for Lc = 0.02 m, 0.05 m, 0.15 m, and 0.30 m.

clear;
format long

set(groot, 'defaultAxesFontName', 'Times New Roman');

%% (a) Ejection-velocity increment distributions
x_dv = logspace(-2, 2.6, 1000);

% Lc = 0.02 m
Lc_dv_002 = 0.02;
logLc_dv_002 = log10(Lc_dv_002);
mu_dv_002 = 0.2 * (-0.3 - 1.4 * (logLc_dv_002 + 1.75)) + 1.85;
sigma_dv_002 = sqrt((0.2 + 0.1333 * (logLc_dv_002 + 3.5))^2 * 0.04 + 0.16);
y_dv_002 = 1 ./ x_dv / log(10) / sqrt(2*pi) / sigma_dv_002 .* ...
    exp(-(log10(x_dv) - mu_dv_002).^2 / (2 * sigma_dv_002^2));

% Lc = 0.05 m
Lc_dv_005 = 0.05;
logLc_dv_005 = log10(Lc_dv_005);
mu_dv_005 = -0.28 * logLc_dv_005 + 1.3;
sigma_dv_005 = sqrt(0.16 + 0.04 * (0.1333 * logLc_dv_005 + 0.66655)^2);
y_dv_005 = 1 ./ x_dv / log(10) / sqrt(2*pi) / sigma_dv_005 .* ...
    exp(-(log10(x_dv) - mu_dv_005).^2 / (2 * sigma_dv_005^2));

% Lc = 0.15 m
Lc_dv_015 = 0.15;
logLc_dv_015 = log10(Lc_dv_015);
alpha_dv_015 = 0.3 + 0.4 * (logLc_dv_015 + 1.2);
mu_dv_015_1 = 1.85 + 0.2 * (-0.6 - 0.318 * (logLc_dv_015 + 1.1));
sigma_dv_015_1 = sqrt(0.16 + 0.04 * (0.1 + 0.2 * (logLc_dv_015 + 1.3))^2);
mu_dv_015_2 = 1.85 + 0.2 * (-1.2);
sigma_dv_015_2 = sqrt(0.16 + 0.04 * 0.5^2);
y_dv_015 = alpha_dv_015 ./ x_dv / log(10) / sqrt(2*pi) / sigma_dv_015_1 .* ...
    exp(-(log10(x_dv) - mu_dv_015_1).^2 / (2 * sigma_dv_015_1^2)) + ...
    (1 - alpha_dv_015) ./ x_dv / log(10) / sqrt(2*pi) / sigma_dv_015_2 .* ...
    exp(-(log10(x_dv) - mu_dv_015_2).^2 / (2 * sigma_dv_015_2^2));

% Lc = 0.30 m
Lc_dv_030 = 0.30;
logLc_dv_030 = log10(Lc_dv_030);
alpha_dv_030 = 0.3 + 0.4 * (logLc_dv_030 + 1.2);
mu_dv_030_1 = (-0.6 - 0.318 * (logLc_dv_030 + 1.1)) * 0.2 + 1.85;
sigma_dv_030_1 = sqrt((0.1 + 0.2 * (logLc_dv_030 + 1.3))^2 + 0.16);
mu_dv_030_2 = (-1.2 - 1.333 * (logLc_dv_030 + 0.7)) * 0.2 + 1.85;
sigma_dv_030_2 = sqrt(0.5^2 * 0.04 + 0.16);
y_dv_030 = alpha_dv_030 ./ x_dv / log(10) / sqrt(2*pi) / sigma_dv_030_1 .* ...
    exp(-(log10(x_dv) - mu_dv_030_1).^2 / (2 * sigma_dv_030_1^2)) + ...
    (1 - alpha_dv_030) ./ x_dv / log(10) / sqrt(2*pi) / sigma_dv_030_2 .* ...
    exp(-(log10(x_dv) - mu_dv_030_2).^2 / (2 * sigma_dv_030_2^2));

%% (b) Area-to-mass-ratio distributions
x_am = logspace(-5, 0, 1000);

% Lc = 0.02 m
Lc_am_002 = 0.02;
logLc_am_002 = log10(Lc_am_002);
mu_am_002 = -0.3 - 1.4 * (logLc_am_002 + 1.75);
sigma_am_002 = 0.2 + 0.1333 * (logLc_am_002 + 3.5);
y_am_002 = 1 ./ x_am / log(10) / sqrt(2*pi) / sigma_am_002 .* ...
    exp(-(log10(x_am) - mu_am_002).^2 / (2 * sigma_am_002^2));

% Lc = 0.05 m
Lc_am_005 = 0.05;
logLc_am_005 = log10(Lc_am_005);
mu_am_005 = -0.3 - 1.4 * (logLc_am_005 + 1.75);
sigma_am_005 = 0.2 + 0.1333 * (logLc_am_005 + 3.5);
y_am_005 = 1 ./ x_am / log(10) / sqrt(2*pi) / sigma_am_005 .* ...
    exp(-(log10(x_am) - mu_am_005).^2 / (2 * sigma_am_005^2));

% Lc = 0.15 m
Lc_am_015 = 0.15;
logLc_am_015 = log10(Lc_am_015);
alpha_am_015 = 0.3 + 0.4 * (logLc_am_015 + 1.2);
mu_am_015_1 = -0.6 - 0.318 * (logLc_am_015 + 1.1);
sigma_am_015_1 = 0.1 + 0.2 * (logLc_am_015 + 1.3);
mu_am_015_2 = -1.2;
sigma_am_015_2 = 0.5;
y_am_015 = alpha_am_015 ./ x_am / log(10) / sqrt(2*pi) / sigma_am_015_1 .* ...
    exp(-(log10(x_am) - mu_am_015_1).^2 / (2 * sigma_am_015_1^2)) + ...
    (1 - alpha_am_015) ./ x_am / log(10) / sqrt(2*pi) / sigma_am_015_2 .* ...
    exp(-(log10(x_am) - mu_am_015_2).^2 / (2 * sigma_am_015_2^2));

% Lc = 0.30 m
Lc_am_030 = 0.30;
logLc_am_030 = log10(Lc_am_030);
alpha_am_030 = 0.3 + 0.4 * (logLc_am_030 + 1.2);
mu_am_030_1 = -0.6 - 0.318 * (logLc_am_030 + 1.1);
sigma_am_030_1 = 0.1 + 0.2 * (logLc_am_030 + 1.3);
mu_am_030_2 = -1.2 - 1.333 * (logLc_am_030 + 0.7);
sigma_am_030_2 = 0.5;
y_am_030 = alpha_am_030 ./ x_am / log(10) / sqrt(2*pi) / sigma_am_030_1 .* ...
    exp(-(log10(x_am) - mu_am_030_1).^2 / (2 * sigma_am_030_1^2)) + ...
    (1 - alpha_am_030) ./ x_am / log(10) / sqrt(2*pi) / sigma_am_030_2 .* ...
    exp(-(log10(x_am) - mu_am_030_2).^2 / (2 * sigma_am_030_2^2));

%% Combined figure
f1 = figure('Color', 'w', ...
    'Position', [50 50 1200 1500], ...
    'WindowStyle', 'normal');
t = tiledlayout(f1, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(t, 1);
plot(ax1, x_dv, y_dv_002, ...
          x_dv, y_dv_005, ...
          x_dv, y_dv_015, ...
          x_dv, y_dv_030, 'LineWidth', 2);
xlabel(ax1, 'Magnitude of ejection-velocity increment (m/s)');
ylabel(ax1, 'Probability density');
title(ax1, '(a) Ejection-velocity increment distribution', ...
    'FontSize', 18, 'FontWeight', 'normal');
legend(ax1, '$L_c = 0.02\;\mathrm{m}$', ...
            '$L_c = 0.05\;\mathrm{m}$', ...
            '$L_c = 0.15\;\mathrm{m}$', ...
            '$L_c = 0.30\;\mathrm{m}$', ...
            'Interpreter', 'latex', 'Location', 'northeast');
xlim(ax1, [0 400]);
ylim(ax1, [0 0.016]);
set(ax1, 'FontSize', 20, 'LineWidth', 1);
grid(ax1, 'on');

ax2 = nexttile(t, 2);
plot(ax2, x_am, y_am_002, ...
          x_am, y_am_005, ...
          x_am, y_am_015, ...
          x_am, y_am_030, 'LineWidth', 2);
xlabel(ax2, 'Area-to-mass ratio ($\mathrm{m^2/kg}$)', 'Interpreter', 'latex');
ylabel(ax2, 'Probability density');
title(ax2, '(b) Area-to-mass ratio distribution', ...
    'FontSize', 18, 'FontWeight', 'normal');
legend(ax2, '$L_c = 0.02\;\mathrm{m}$', ...
            '$L_c = 0.05\;\mathrm{m}$', ...
            '$L_c = 0.15\;\mathrm{m}$', ...
            '$L_c = 0.30\;\mathrm{m}$', ...
            'Interpreter', 'latex', 'Location', 'northeast');
xlim(ax2, [0 1]);
ylim(ax2, [0 8]);
set(ax2, 'FontSize', 20, 'LineWidth', 1);
grid(ax2, 'on');

drawnow;

outfile = fullfile(pwd, 'fig_delta_v_S_m_combined.png');
exportgraphics(f1, outfile, 'Resolution', 600, 'BackgroundColor', 'white');
