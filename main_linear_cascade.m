clear;
rng(0);
%% load original model
d = 5; % # of species
load(['Data/models/linear_cascade_' num2str(d) '.mat'])



M = 500;

%% generate observations from original model
observed_ind = d;
estimated_ind = 1;
T = 5;
[Y_obs, t_obs, Z_true, t_true] = generate_observations(model,Z0,T,observed_ind);

%% MP
proj_ind = [estimated_ind; observed_ind];

max_S = 10;

max_Z = max_S*ones(2,1);
dt = 1e-2;


disp('running MP...')

tic
model_mp = markovian_projection_extrap(model, Z0, T, proj_ind, ...
    max_Z, M , dt);
toc

filter_ffsp_mp = FFSP(model_mp, 2, max_S);
filter_ffsp_mp.dt = 1e-2;

tic 
disp('running FFSP for MP model...')
pi_mp = filter_ffsp_mp.fit(Y_obs, t_obs, Z0(proj_ind));
toc


%% Particle filter

filter_pf = PF(model, d, M);
filter_pf.dt = 1e-2;

disp('running PF:')
tic
[X, w, t] = filter_pf.fit(Y_obs, t_obs, Z0);
toc

%% estimate a_tilde 

disp('estimating a_tile for FMP ...')
tic

nt = length(t);

a_tilde = zeros(model.r, max_S+1, nt);

counts = zeros(max_S+1, nt);

for i = 1:M
    for it = 1:nt
        Xi = X(:, i, it);

        if Xi(estimated_ind) > max_S 
            continue
        end

        a_tilde(:, Xi(estimated_ind)+1, it) = ...
        a_tilde(:, Xi(estimated_ind)+1, it) + w(i, it)*model.a(Xi, t(it));

        counts(Xi(estimated_ind)+1, it) = counts(Xi(estimated_ind)+1, it) + w(i, it);
    end
end

for ir = 1:model.r
    a_tilde(ir,:,:) = squeeze(a_tilde(ir,:,:)) ./ counts;
end


a_tilde = extrapolate_a_bar_1d(a_tilde, T, 'nearest');

a = @(x, t_) a_tilde(:, x(1)+1, find(t >= t_, 1));
model_fmp = SRN(model.V(proj_ind, :), a);
model_fmp.is_homogeneous = false;

toc

%% FMP

filter_ffsp_fmp = FFSP(model_fmp, 2, max_S);
filter_ffsp_fmp.dt = 1e-2;

disp('running FFSP for FMP model...')
tic 
pi_fmp = filter_ffsp_fmp.fit(Y_obs, t_obs, Z0(proj_ind));
toc



%% FFSP 

filter_ffsp_fmp= FFSP(model, d, max_S*ones(d-1, 1));
filter_ffsp_fmp.dt  = 1e-2;
filter_ffsp_fmp.show_progress = true;

disp('running FFSP for full model...')
tic 
pi_ffsp  = filter_ffsp_fmp.fit(Y_obs, t_obs, Z0);
toc



%% save 
filename = ['Data/results/lin_cascade_' num2str(d) '_sol.mat'];
save(filename)
disp(['Results saved to ' filename]);
disp(' ')




%% TAIL ESTIMATION

disp(' ')
disp('TAIL ESTIMATION')


d = 3;
tail_estimation_linear_cascade;

d = 5;
tail_estimation_linear_cascade;


%% plot
plot_lin_cascade_sol