
rng(0)

%%
d = 5;
disp(['d = ' num2str(d)]);

r = 8; % we estimate p := P(Z_1(T) >= r | Z_n(s), s <= T)

M_rep = 100; % # of repetitions for each sample size


load(['Data/models/linear_cascade_' num2str(d) '.mat'])

%% generate observations from original model
observed_ind = d;
estimated_ind = 1;
T = 5;
[Y_obs, t_obs, Z, t_true] = generate_observations(model,Z0,T,observed_ind);


%% FFSP ('true' solution)

T = 5;
max_X = 20*ones(model.d-1, 1);

filter = FFSP(model, observed_ind, max_X);
filter.show_progress = true;
[pi, t] = filter.fit(Y_obs, t_obs, Z0);

pi1 = squeeze(sum(pi, 2:numel(max_X)));
p_ffsp = sum(pi1(r+1:end, end));




%% PF,FMP,MP with different sample size

M_arr = 2.^(9:16); 
dt = 1e-2;

proj_ind = [estimated_ind; observed_ind];

p_mp = nan(numel(M_arr), M_rep);
p_pf = nan(numel(M_arr), M_rep);
p_fmp = nan(numel(M_arr), M_rep);

filename = ['Data/results/lin_cascade_' num2str(d) '_tail_estimation.mat'];

for idx = 1:length(M_arr)
    M = M_arr(idx);

    parfor irep = 1:M_rep
        disp(['M = ' num2str(M), 'rep = ', num2str(irep)])        
        %% run PF
        filter_pf = PF(model, observed_ind, M);
        filter_pf.dt = dt;
        
        tic
        [X, w, t_pf] = filter_pf.fit(Y_obs, t_obs, Z0);

        % PF estimation of QoI
        p_pf(idx, irep) = sum(w(X(1, :, end) >= r, end));
        toc

        %% estimate a_tilde
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
                    a_tilde(:, Xi(estimated_ind)+1, it) + w(i, it)*model.a(Xi, t_pf(it));

                counts(Xi(estimated_ind)+1, it) = counts(Xi(estimated_ind)+1, it) + w(i, it);
            end
        end

        for ir = 1:model.r
            a_tilde(ir,:,:) = squeeze(a_tilde(ir,:,:)) ./ counts;
        end


        a_tilde = extrapolate_a_bar_1d(a_tilde, T, 'linear');

        a = @(x, t_) a_tilde(:, x(1)+1, find(t_pf >= t_, 1));
        model_fmp = SRN(model.V([1 d], :), a);
        model_fmp.is_homogeneous = false;
        %% FFSP for FMP model

        tic
        filter_fmp = FFSP(model_fmp, 2, max_X(1));
        filter_fmp.dt = dt;
        [pi_fmp, t_fmp] = filter_fmp.fit(Y_obs, t_obs, Z0([1 d]));
        toc

        p_fmp(idx, irep) = sum(pi_fmp(r+1:end, end));
        
    end


    %
    disp(['abs error PF = ' num2str(mean(abs(p_pf(idx,:) - p_ffsp)))])
    disp(['abs error FMP= ' num2str(mean(abs(p_fmp(idx,:) - p_ffsp)))])

    
    save(filename, '-v7.3')
    disp(['Intermediate results saved to ' filename]);
end

save(filename, '-v7.3')
disp(['All results saved to ' filename]);

