
rng(0)
%%
if isempty(d)
    d = 5; % # of species
end
disp(['d = ' num2str(d)]);

r = 8; % we estimate p := P(Z_1(T) >= r | Z_n(s), s <= T)

M_rep = 30; % # of repetitions for each sample size


load(['Data/models/linear_cascade_' num2str(d) '.mat'])

%% generate observations from original model
observed_ind = d;
estimated_ind = 1;
T = 5;
[Y_obs, t_obs, Z, t_true] = generate_observations(model,Z0,T,observed_ind);


%% FFSP

T = 5;
max_X = 10*ones(model.d-1, 1);

filter = FFSP(model, observed_ind, max_X);
filter.show_progress = true;
[pi, t] = filter.fit(Y_obs, t_obs, Z0);


pi1 = squeeze(sum(pi, 2:numel(max_X)));
p_ffsp = sum(pi1(r+1:end, end));




%% PF,FMP,MP with different sample size
figure; 

M_arr = 2.^(10:18);

proj_ind = [estimated_ind; observed_ind];

p_mp = nan(numel(M_arr), M_rep);
p_pf = nan(numel(M_arr), M_rep);
p_fmp = nan(numel(M_arr), M_rep);

filename = ['Data/results/lin_cascade_' num2str(d) '_tail_estimation.mat'];

for idx = 1:length(M_arr)
    M = M_arr(idx);

    for irep = 1:M_rep
        disp(['M = ' num2str(M), 'rep = ', num2str(irep)])
        %% MP
        dt = 5e-2;
        max_Z = [15; 15];

        tic
        model_mp = markovian_projection_extrap(model, Z0, t_obs(end), ...
            proj_ind, max_Z, M , dt);


        filter_mp = FFSP(model_mp, 2, max_X(1));
        [pi_mp, t_mp] = filter_mp.fit(Y_obs, t_obs, Z0([1 d]));

        toc
        p_mp(idx, irep) = sum(pi_mp(r+1:end, end));

        %% PF
        filter_pf = PF(model, observed_ind, M);
        filter_pf.dt = 5e-2;

        tic
        [X, w, t_pf] = filter_pf.fit(Y_obs, t_obs, Z0);

        %estimate a_tilde
        max_X1 = max_X(1);
        nt = length(t_pf);
        a_tilde = zeros(model.r, max_X1+1, nt);
        counts = zeros(max_X1+1, nt);

        for i = 1:M
            for it = 1:nt
                Xi = X(:, i, it);

                if Xi(estimated_ind) > max_X1
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


        a_tilde = extrapolate_a_bar_1d(a_tilde, T, 'linear');
        %a_tilde(isnan(a_tilde)) = 10;

        a = @(x, t_) a_tilde(:, x(1)+1, find(t_pf >= t_, 1));
        model_fmp = SRN(model.V([1 d], :), a);
        model_fmp.is_homogeneous = false;

        toc
        %% FFSP for FMP model

        tic
        filter_fmp = FFSP(model_fmp, 2, max_X(1));
        [pi_fmp, t_fmp] = filter_fmp.fit(Y_obs, t_obs, Z0([1 d]));
        toc

        p_fmp(idx, irep) = sum(pi_fmp(r+1:end, end));

        %% PF estimation of QoI

        p_pf(idx, irep) = sum(w(X(1, :, end) >= r, end));
        clear X w
        
        


    end
    p_mp(idx, isnan(p_mp(idx, :))) = 0;
    p_fmp(idx, isnan(p_mp(idx, :))) = 0;

    %
    disp(['abs error MP = ' num2str(abs(mean(p_mp(idx,:)) - p_ffsp))])
    disp(['abs error PF = ' num2str(abs(mean(p_pf(idx,:)) - p_ffsp))])
    disp(['abs error FMP= ' num2str(abs(mean(p_fmp(idx,:)) - p_ffsp))])

    
    save(filename)
    disp(['Intermediate results saved to ' filename]);
end

save(filename)
disp(['All results saved to ' filename]);

