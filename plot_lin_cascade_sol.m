%% This script plots solution for linear cascade model
% (figures 5-8 in the paper)
% you need to run 'main' first to get *.mat files with simulation results

%% set plot params
LineWidth = 4;

set(0,'DefaultLineLineWidth', LineWidth)
set(0, 'defaultAxesFontSize', 27)
set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on', ...
    'DefaultAxesZGrid','on');
set(groot, 'defaultTextInterpreter', 'latex'); 
set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(groot, 'defaultLegendInterpreter', 'latex'); 
set(groot, 'defaultTextFontWeight', 'bold'); 


WindowStyle = 'normal';
%% observations + estimated mean (Figure 5)
figure('WindowStyle', WindowStyle, 'Units', 'Inches', ...
    'Position', [0, 0, 18, 6]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%Observations
load Data/results/lin_cascade_5_sol.mat

subplot1 = nexttile;
hold on;
stairs(t_obs, Y_obs(1, :), 'LineWidth', LineWidth, ...
    'DisplayName', '\textbf{Observed $Z_d$}')
legend('Location', 'northeast')
xlabel('\textbf{Time}');
ylabel('\textbf{Copy number}');
xlim([0 5])
%text(0.05, 0.9, '\textbf{(a)}', 'Units', 'normalized'); 
title('\textbf{(a)}')



%Mean value 
mean_pf = sum(squeeze(X(1, :, :)) .* w, 1);
mean_ffsp = (0:max_S) * squeeze(sum(pi_ffsp, 2:(d-1)));
mean_mp = (0:max_S) * pi_mp;
mean_fmp = (0:max_S) * pi_fmp;

subplot2 = nexttile; 
hold on;
stairs(t_true, Z_true(1, :), 'black', 'LineWidth', LineWidth, ...
    'DisplayName', '\textbf{Hidden trajectory}');
plot(t, mean_ffsp, 'b', 'DisplayName', '\textbf{FFSP}'); 
plot(t, mean_mp, 'r', 'DisplayName', '\textbf{Unconditional MP filter}'); 
plot(t, mean_fmp, 'g', 'DisplayName', '\textbf{Conditional MP filter}'); 
L = legend('Location', 'northeast', 'NumColumns', 1);
L.AutoUpdate = 'off'; 
xlabel('\textbf{Time}');
ylabel('\textbf{Copy number}');
%xline(time_hist_1, '--', 'LineWidth', LineWidth);
%xline(time_hist_2, '--', 'LineWidth', LineWidth);
%text(0.05, 0.9, '\textbf{(b)}', 'Units', 'normalized'); 
xlim([0 5])
title('\textbf{(b)}')

% save
exportgraphics(gcf, 'linear_cascade_sol.pdf', 'ContentType', 'vector')
saveas(gcf, 'linear_cascade_sol.png');



%% PMF in log scale (Figure 6)
figure('WindowStyle', WindowStyle, 'Units', 'Inches', ...
    'Position', [0, 0, 10, 8]);
hold on;
plot(0:10, squeeze(sum(pi_ffsp(:, :, :, :, end), 2:4)), '-ob', ...
    'MarkerSize', 14, 'MarkerFaceColor', 'b', ...
    'DisplayName', '\textbf{FFSP}')

% PF 
pi_pf = zeros(11,1);
for i = 1:M
    ind = 1+X(1,i,end);
    pi_pf(ind) = pi_pf(ind) + w(i,end);
end

plot(0:10, pi_pf, ...
    '-om', 'DisplayName','\textbf{Particle Filter}') 


plot(0:10, pi_mp(:, it), '-or', 'DisplayName', '\textbf{Unconditional MP filter}')
plot(0:10, pi_fmp(:, it), '-og', 'DisplayName', '\textbf{Conditional MP filter}')
legend('Location', 'southwest', 'NumColumns', 1);
xlabel('\textbf{Copy number}');
ylabel('\textbf{Probability}');
xlim([0 10]);
yscale('log')
%text(5, 0.2, '\textbf{d = 5}', 'FontSize', 22); 

exportgraphics(gcf, 'linear_cascade_log_distr.pdf', 'ContentType', 'vector')
saveas(gcf, 'linear_cascade_log_distr.png');


%% Convergence of error in estimating tail (Figure 7)

figure('WindowStyle', WindowStyle, 'Units', 'Inches', ... 
    'Position', [0, 0, 10, 8]);
tiledlayout('TileSpacing', 'compact', 'Padding', 'compact');
hold on


load Data/results/lin_cascade_5_tail_estimation.mat


% PF 
mean_error_pf = mean(abs(p_pf - p_ffsp), 2)/p_ffsp;
ci_error_pf = 1.96./sqrt(M_rep).*std(abs(p_pf - p_ffsp)')/p_ffsp;
errorbar(M_arr, mean_error_pf, ci_error_pf, ...
    '-Om', 'LineWidth', LineWidth, 'DisplayName','\textbf{Particle Filter}') 

% FMP
mean_error_fmp = mean(abs(p_fmp - p_ffsp), 2)/p_ffsp;
ci_error_fmp = 1.96./sqrt(M_rep).*std(abs(p_fmp - p_ffsp)')/p_ffsp;
errorbar(M_arr, mean_error_fmp, ci_error_fmp, ...
    '-Og', 'LineWidth', LineWidth, 'DisplayName','\textbf{Conditional MP filter}')  

% ref line
plot(M_arr, 25./sqrt(M_arr), ...
    '--', 'LineWidth', 2, 'Color', 'black', ...
    'DisplayName', '$\propto M^{-1/2}$') 


ylim([3e-2 10])
xlim([512 65536])
xscale('log')
yscale('log')
xlabel('\textbf{Sample size, M}')
ylabel('\textbf{Relative Error}')
legend()



exportgraphics(gcf, 'linear_cascade_errors.pdf', 'ContentType', 'vector')
saveas(gcf, 'linear_cascade_errors.png');




%% CPU times (Figure 8)
set(0, 'defaultAxesFontSize', 18)

figure('WindowStyle', WindowStyle, 'Units', 'Inches', ...
    'Position', [0, 0, 8, 6]);
hold on;

% this cpu times obtained on Intel Xeon 6230R (sequential implantation)
d = 4:8;
mp_time = [4.1, 4.2, 4.1, 4.1, 4.];
fmp_time = [1.8, 1.7, 1.9, 2.0, 1.9];
ffsp_time = [3, 42, 366, 3915, 38166];

X = d;
Y = [ffsp_time; mp_time; fmp_time]';
b = bar(X, Y, 'grouped', 'BarWidth', 0.9);
b(1).FaceColor = 'b';
b(2).FaceColor = 'r';
b(3).FaceColor = 'g';
xlabel('\textbf{Dimensionality, d}');
ylabel('\textbf{CPU Time (s)}');
legend('\textbf{FFSP for full model}', '\textbf{Unconditional MP filter}', ...
    '\textbf{Conditional MP filter}', 'Location', 'northwest');
yscale('log')


exportgraphics(gcf, 'linear_cascade_cpu_times.pdf', 'ContentType', 'vector')
saveas(gcf, 'linear_cascade_cpu_times.png');



