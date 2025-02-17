%% This script plots solution for linear cascade model
% (figures 5-7 in the paper)
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
plot(t, mean_mp, 'r', 'DisplayName', '\textbf{Standard MP}'); 
plot(t, mean_fmp, 'g', 'DisplayName', '\textbf{Filtered MP}'); 
L = legend('Location', 'northeast', 'NumColumns', 1);
L.AutoUpdate = 'off'; 
xlabel('\textbf{Time}');
ylabel('\textbf{Copy number}');
%xline(time_hist_1, '--', 'LineWidth', LineWidth);
%xline(time_hist_2, '--', 'LineWidth', LineWidth);
%text(0.05, 0.9, '\textbf{(b)}', 'Units', 'normalized'); 
xlim([0 5])
title('\textbf{(b)}')

saveas(gcf, 'lin_cascade_sol.png');


%% Convergence of error in estimating tail (Figure 7)
set(0, 'defaultAxesFontSize', 27)

figure('WindowStyle', WindowStyle, 'Units', 'Inches', ... 
    'Position', [0, 0, 18, 8]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');



%% errors d=3

load Data/results/lin_cascade_3_tail_estimation.mat


subplot1 = nexttile;
hold on

% PF 
mean_error_pf = mean(abs(p_pf - p_ffsp), 2)/p_ffsp;
ci_error_pf = 1.96./sqrt(M_rep).*std(abs(p_pf - p_ffsp)')/p_ffsp;
errorbar(M_arr, mean_error_pf, ci_error_pf, ...
    '-Om', 'LineWidth', LineWidth, 'DisplayName','\textbf{Particle Filter}') 

% FMP
mean_error_fmp = mean(abs(p_fmp - p_ffsp), 2)/p_ffsp;
ci_error_fmp = 1.96./sqrt(M_rep).*std(abs(p_fmp - p_ffsp)')/p_ffsp;
errorbar(M_arr, mean_error_fmp, ci_error_fmp, ...
    '-Og', 'LineWidth', LineWidth, 'DisplayName','\textbf{Filtered MP}')  

% MP
mean_error_mp = mean(abs(p_mp - p_ffsp), 2)/p_ffsp;
ci_error_mp = 1.96./sqrt(M_rep).*std(abs(p_mp - p_ffsp)')/p_ffsp;

ind = find(ci_error_mp' >  mean_error_mp); %plot CI only if CI < mean
ci_error_mp(ind) = NaN;
errorbar(M_arr, mean_error_mp, ci_error_mp, ...
     '-Or', 'LineWidth', LineWidth, 'DisplayName','\textbf{Standard MP}')  

% ref line
plot(M_arr, 100./sqrt(M_arr), ...
    '--', 'LineWidth', 2, 'Color', 'black', ...
    'DisplayName', '$\propto M^{-0.5}$') 


ylim([1e-2 1e2])
xlim([M_arr(1) M_arr(end)])
title('\textbf{(b)}')
text(M_arr(floor(length(M_arr)/2)), 3e-2, '\textbf{d = 5}');
xscale('log')
yscale('log')
xlabel('\textbf{Sample size, M}')
legend()

%% errors d=5 

load Data/results/lin_cascade_5_tail_estimation.mat

subplot2 = nexttile;
hold on

% PF 
mean_error_pf = mean(abs(p_pf - p_ffsp), 2)/p_ffsp;
ci_error_pf = 1.96./sqrt(M_rep).*std(abs(p_pf - p_ffsp)')/p_ffsp;
errorbar(M_arr, mean_error_pf, ci_error_pf, ...
    '-Om', 'LineWidth', LineWidth, 'DisplayName','\textbf{Particle Filter}') 

% FMP
mean_error_fmp = mean(abs(p_fmp - p_ffsp), 2)/p_ffsp;
ci_error_fmp = 1.96./sqrt(M_rep).*std(abs(p_fmp - p_ffsp)')/p_ffsp;
errorbar(M_arr, mean_error_fmp, ci_error_fmp, ...
    '-Og', 'LineWidth', LineWidth, 'DisplayName','\textbf{Filtered MP}')  

% MP
mean_error_mp = mean(abs(p_mp - p_ffsp), 2)/p_ffsp;
ci_error_mp = 1.96./sqrt(M_rep).*std(abs(p_mp - p_ffsp)')/p_ffsp;

ind = find(ci_error_mp' >  mean_error_mp); %plot CI only if CI < mean
ci_error_mp(ind) = NaN;
errorbar(M_arr, mean_error_mp, ci_error_mp, ...
     '-Or', 'LineWidth', LineWidth, 'DisplayName','\textbf{Standard MP}')  

% ref line
plot(M_arr, 100./sqrt(M_arr), ...
    '--', 'LineWidth', 2, 'Color', 'black', ...
    'DisplayName', '$\propto M^{-0.5}$') 


ylim([1e-2 1e2])
xlim([M_arr(1) M_arr(end)])
title('\textbf{(b)}')
text(M_arr(floor(length(M_arr)/2)), 3e-2, '\textbf{d = 5}');
xscale('log')
yscale('log')
xlabel('\textbf{Sample size, M}')
legend()



saveas(gcf, 'linear_cascade_errors.png');


%% CPU times (Figure 6)
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
legend('\textbf{FFSP for full model}', '\textbf{Standard MP}', ...
    '\textbf{Filtered MP}', 'Location', 'northwest');
yscale('log')


saveas(gcf, 'cpu_times.png');



