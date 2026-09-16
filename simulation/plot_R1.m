%% Matching-metric ablation: JS divergence comparison
% Load the saved R1 results and plot the JS-divergence histories for
% baseline, +i, +Omega, and +i+Omega.
% Plotting style follows plot_SC_4.m.

clear;
clc;

%% ==== data file ====
zhongzi = 2;
ns_file = 21000;
data_dir = fullfile(pwd,'data');
filename = fullfile(data_dir,sprintf('zhongzi=%d_ns=%d_R1.mat',zhongzi,ns_file));

% Fallback: allow the .mat file to be placed directly in the current folder.
if ~isfile(filename)
    filename = fullfile(pwd,sprintf('zhongzi=%d_ns=%d_R1.mat',zhongzi,ns_file));
end

if ~isfile(filename)
    error('Cannot find R1 result file: %s',filename);
end

S = load(filename);

JS_sandu_R1   = S.JS_sandu_R1;
JS_pingjun_R1 = S.JS_pingjun_R1;
JS_zuida_R1   = S.JS_zuida_R1;
Ts             = S.Ts;
ns             = S.ns;

% Compatibility in case the array orientation is reversed.
if size(JS_sandu_R1,1) ~= 4 && size(JS_sandu_R1,2) == 4
    JS_sandu_R1 = JS_sandu_R1.';
end

mode_name = {'Baseline','+$i$','+$\Omega$','+$i+\Omega$'};

%% ==== output numerical results ====
fprintf('\nMatching-metric ablation results:\n');
fprintf('------------------------------------------------------\n');
fprintf('%-14s  %-14s  %-14s\n','Case','Mean JS','Max JS');
fprintf('------------------------------------------------------\n');
fprintf('%-14s  %.8f      %.8f\n','Baseline',JS_pingjun_R1(1),JS_zuida_R1(1));
fprintf('%-14s  %.8f      %.8f\n','+i',JS_pingjun_R1(2),JS_zuida_R1(2));
fprintf('%-14s  %.8f      %.8f\n','+Omega',JS_pingjun_R1(3),JS_zuida_R1(3));
fprintf('%-14s  %.8f      %.8f\n','+i+Omega',JS_pingjun_R1(4),JS_zuida_R1(4));
fprintf('------------------------------------------------------\n\n');

%% ==== JS divergence histories ====
set(groot,'defaultAxesFontName','Times New Roman');

figure('Color','w','Position',[100 100 720 450],'WindowStyle','normal');
tl_js = tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
ax = nexttile(tl_js,1);

% Same blue as plot_SC_4.m for the baseline; the other colors distinguish
% the three augmented matching metrics.
c1 = [0, 82, 155]/255;
c2 = [213, 94, 0]/255;
c3 = [0, 135, 95]/255;
c4 = [204, 121, 167]/255;

t_days = (0:ns)*Ts/86400;

plot(t_days,JS_sandu_R1(1,:),'Color',c1,'LineWidth',1); hold on;
plot(t_days,JS_sandu_R1(2,:),'Color',c2,'LineWidth',1);
plot(t_days,JS_sandu_R1(3,:),'Color',c3,'LineWidth',1);
plot(t_days,JS_sandu_R1(4,:),'Color',c4,'LineWidth',1); hold off;

xlabel('Time (days)');
ylabel('JS divergence');
legend(mode_name,'Interpreter','latex','Location','northeast','FontSize',12);

% Keep the axis treatment consistent with plot_SC_4.m.
% In particular, do not force xlim to the exact final epoch; MATLAB chooses
% clean tick-aligned limits automatically, as in the original JS figure.
ylim([0 0.05]);
set(ax,'FontSize',15,'LineWidth',1,'XGrid','off','YGrid','on');

drawnow;

%% ==== export ====
outfile_js = fullfile(pwd,'fig_R1_matching_metric_js.png');
exportgraphics(tl_js,outfile_js,'Resolution',600,'BackgroundColor','white');

savefig(gcf,fullfile(pwd,'fig_R1_matching_metric_js.fig'));

fprintf('Figure saved to:\n%s\n',outfile_js);
