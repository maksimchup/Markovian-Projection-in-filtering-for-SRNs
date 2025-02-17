clear;
rng(0);
%% load original model
load Data/models/bistable_gene_network.mat

M = 1e3;

%% generate observations from original modeel
T = 5;
observed_ind = (7:8)';
estimated_ind = 6;

[Y_obs,t_obs,Z_true,t_true]=generate_observations(model,Z0,T,observed_ind);



%% MP
proj_ind = [estimated_ind; observed_ind];


max_mRNA = 30;
max_protein1 = max(Y_obs(1, :));
max_protein2 = max(Y_obs(2, :));

max_Z = [max_mRNA; max_protein1; max_protein2];
dt = 1e-2;

disp('running MP...')
tic
model_mp = markovian_projection_extrap(model, Z0, T, proj_ind, max_Z,M,dt);
toc



filter_ffsp_mp = FFSP(model_mp, [2;3], max_mRNA);
filter_ffsp_mp.dt = 1e-2;

tic 
disp('running FFSP for MP model...')
pi_mp = filter_ffsp_mp.fit(Y_obs, t_obs, Z0(proj_ind));
toc


%% Particle filter

filter_pf = PF(model, [7;8], M);
filter_pf.dt = 1e-2;

disp('running PF:')
tic
[X, w, t] = filter_pf.fit(Y_obs, t_obs, Z0);
toc

%% estimate a_tilde 


tic
disp('estimating a_tile for FMP ...')

nt = length(t);

a_tilde = zeros(model.r, max_mRNA+1, nt);

counts = zeros(max_mRNA+1, nt);

for i = 1:M
    for it = 1:nt
        Xi = X(:, i, it);

        if Xi(6) > max_mRNA 
            continue
        end

        a_tilde(:, Xi(6)+1, it) = ...
        a_tilde(:, Xi(6)+1, it) + w(i, it)*model.a(Xi, t(it));

        counts(Xi(6)+1, it) = counts(Xi(6)+1, it) + w(i, it);
    end
end

for ir = 1:model.r
    a_tilde(ir,:,:) = squeeze(a_tilde(ir,:,:)) ./ counts;
end


a_tilde = extrapolate_a_bar_1d(a_tilde, T, 'linear');

a = @(x, t_) a_tilde(:, x(1)+1, find(t >= t_, 1));
model_fmp = SRN(model.V(proj_ind, :), a);
model_fmp.is_homogeneous = false;

toc

%% FMP


filter_ffsp_fmp = FFSP(model_fmp, [2;3], max_mRNA);
filter_ffsp_fmp.dt = 1e-2;

disp('running FFSP for FMP model...')

tic 
pi_fmp = filter_ffsp_fmp.fit(Y_obs, t_obs, Z0(proj_ind));
toc



%% FFSP 

filter_ffsp_full= FFSP(model, [7;8], [1; 1; 1; 1; max_mRNA; max_mRNA]);
filter_ffsp_full.dt  = 1e-2;
filter_ffsp_full.show_progress = true;

disp('running FFSP for full model...')
tic 
pi_ffsp  = filter_ffsp_full.fit(Y_obs, t_obs, Z0);
toc



%% save
save('Data/results/bistable_network_sol.mat', '-v7.3');

%% plot
plot_bistable_network_sol;
