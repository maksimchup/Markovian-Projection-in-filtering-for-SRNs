%% This script plots solution for bistable network 
% (figure 3 in the paper)
% you need to run 'main' first to get *.mat files with simulation results
%% set plot params

time_hist_1 = 2; % time for hist1
time_hist_2 = 4.5; % time for hist2

LineWidth = 4;

set(0,'DefaultLineLineWidth', LineWidth)
set(0, 'defaultAxesFontSize', 27)
set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on','DefaultAxesZGrid','on');
set(groot, 'defaultTextInterpreter', 'latex'); 
set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(groot, 'defaultLegendInterpreter', 'latex'); 
set(groot, 'defaultTextFontWeight', 'bold'); 

WindowStyle = 'docked';

%%
figure('WindowStyle', WindowStyle, 'Units', 'Inches', 'Position', [0, 0, 18, 12]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%% Observations vs time

load Data/results/bistable_network_sol.mat

subplot1 = nexttile;
hold on;
stairs(t_obs, Y_obs(1, :), 'LineWidth', LineWidth, 'DisplayName', '\textbf{Protein 1}')
stairs(t_obs, Y_obs(2, :), 'LineWidth', LineWidth, 'DisplayName', '\textbf{Protein 2}')
legend('Location', 'northwest')
xlabel('\textbf{Time}');
ylabel('\textbf{Copy number}');
xlim([0 5])
%text(0.05, 0.9, '\textbf{(a)}', 'Units', 'normalized'); 
title('\textbf{(a)}')

%% Mean value vs time
 
mean_pf = sum(squeeze(X(6, :, :)) .* w, 1);
mean_ffsp = (0:max_mRNA) * squeeze(sum(pi_ffsp, 1:5));
mean_mp = (0:max_mRNA) * pi_mp;
mean_fmp = (0:max_mRNA) * pi_fmp;

subplot2 = nexttile; 
hold on;
stairs(t_true, Z_true(estimated_ind, :), 'black', 'LineWidth', LineWidth, 'DisplayName', '\textbf{Hidden trajectory}');
plot(t, mean_ffsp, 'b', 'DisplayName', '\textbf{FFSP}'); 
plot(t, mean_mp, 'r', 'DisplayName', '\textbf{Standard MP}'); 
plot(t, mean_fmp, 'g', 'DisplayName', '\textbf{Filtered MP}'); 
L = legend('Location', 'southeast', 'NumColumns', 1);
L.AutoUpdate = 'off'; 
xlabel('\textbf{Time}');
ylabel('\textbf{Copy number}');
xline(time_hist_1, '--', 'LineWidth', LineWidth);
xline(time_hist_2, '--', 'LineWidth', LineWidth);
%text(0.05, 0.9, '\textbf{(b)}', 'Units', 'normalized'); 
xlim([0 5])
title('\textbf{(b)}')



%% Hist 1

it = find(t >= time_hist_1, 1);

subplot3 = nexttile;
hold on;
plot(0:max_mRNA, squeeze(sum(pi_ffsp(:, :, :, :, :, :, it), [1 2 3 4 5])), '-ob', 'DisplayName', '\textbf{FFSP}')
plot(0:max_mRNA, pi_mp(:, it), '-or', 'DisplayName', '\textbf{Standard MP}')
plot(0:max_mRNA, pi_fmp(:, it), '-og', 'DisplayName', '\textbf{Filtered MP}')
legend('Location', 'northeast', 'NumColumns', 1);
xlabel('\textbf{Copy number}');
ylabel('\textbf{Probability}');
xlim([0 10]);
%text(0.05, 0.9, '\textbf{(c)}', 'Units', 'normalized'); 
title('\textbf{(c)}')

%% Hist 2
it = find(t >= time_hist_2, 1);

subplot4 = nexttile;
hold on;
plot(0:max_mRNA, squeeze(sum(pi_ffsp(:, :, :, :, :, :, it), [1 2 3 4 5])), '-ob', 'DisplayName', '\textbf{FFSP}')
plot(0:max_mRNA, pi_mp(:, it), '-or', 'DisplayName', '\textbf{Standard MP}')
plot(0:max_mRNA, pi_fmp(:, it), '-og', 'DisplayName', '\textbf{Filtered MP}')
legend('Location', 'northeast', 'NumColumns', 1);
xlabel('\textbf{Copy number}');
ylabel('\textbf{Probability}');
xlim([0 10]);
%text(0.05, 0.9, '\textbf{(d)}', 'Units', 'normalized'); 
title('\textbf{(d)}')


%% Arrows to hists

pos2 = get(subplot2, 'Position'); 
pos3 = get(subplot3, 'Position'); 
pos4 = get(subplot4, 'Position'); 

% to hist 1
x_start_1 = pos2(1) + pos2(3) * (time_hist_1 - min(t)) / (max(t) - min(t));
y_start_1 = pos2(2); 
x_end_1 = pos3(1) + pos3(3) / 1.5; 
y_end_1 = pos3(2) + pos3(4) + 0.02; 
annotation('arrow', [x_start_1, x_end_1], [y_start_1, y_end_1], 'Color', 'black', 'LineWidth', 3);

% to hist 2
x_start_2 = pos2(1) + pos2(3) * (time_hist_2 - min(t)) / (max(t) - min(t)); 
y_start_2 = pos2(2); 
x_end_2 = pos4(1) + pos4(3) / 1.5; 
y_end_2 = pos4(2) + pos4(4) + 0.02; 
annotation('arrow', [x_start_2, x_end_2], [y_start_2, y_end_2], 'Color', 'black', 'LineWidth', 3);

saveas(gcf, 'bistable_network_sol.png');