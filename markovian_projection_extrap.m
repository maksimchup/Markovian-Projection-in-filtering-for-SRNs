function model_proj = markovian_projection_extrap(model, Z0, T, ...
    proj_ind, max_Z, M, dt)
%Markovian projection for SRN
% new propensities are estimated with MC and linear space extrapolation
% -------------------------------------------------------------------------
%INPUT
% model         : Stochastic Reaction Network         | object of class SRN
% Z0            : initial state at time 0             | (d,1) array
% T             : final time
% proj_ind      : indices of projected species        | (~,d) array
% max_Z         : 
% M             : # of samples 
% dt            : time step for array
% -------------------------------------------------------------------------
%OUTPUT
% model_reduced : projected SRN                       | object of class SRN
% -------------------------------------------------------------------------

%% Project stoichiometric matrix

V_reduced = model.V(proj_ind, :);
% reactions with zero stoichiometric vector
reactions_ind = any(V_reduced ~= 0, 1);
V_reduced = V_reduced(:, reactions_ind);

%% 

d_reduced = length(proj_ind);
r_reduced = size(V_reduced, 2); % # of reactions in reduced model
nt = ceil(T/dt);                % # of time nodes

a_bar = zeros([r_reduced, max_Z'+1, nt]);
counts = zeros([1, max_Z'+1, nt]);

%% sample with SSA and estimate propensities a_bar
for im = 1:M
    [t_path, Z_path] = model.sample_exact_path(Z0, T);

    for it = 1:nt % save values on on this time grid
        t = min(dt*(it-1), T);        
        Z = Z_path(:, find(t_path >= t, 1));
        Z_reduced = Z(proj_ind);

        if any(Z_reduced > max_Z)
            continue
        end

        a = model.a(Z, t); 
        a_reduced = a(reactions_ind);
        
        idx_Z = num2cell(Z_reduced+1);
        
        a_bar(:, idx_Z{:}, it) = a_bar(:, idx_Z{:}, it) + a_reduced;
        

        counts(1, idx_Z{:}, it) = counts(1, idx_Z{:}, it) + 1;
    end
end

a_bar = a_bar ./ counts;

%% extrapolate propensities a_bar

if length(proj_ind) == 2
    a_bar = extrapolate_a_bar_2d(a_bar, T, 'linear');
elseif lentgh(proj_ind) == 3
    a_bar = extrapolate_a_bar_3d(a_bar, T, 'linear');
else
    warning(['Extrapolation for a_bar for ' num2str(length(proj_ind)) ...
        ' projected species is not implimented! Setting unavailable ' ...
        'values to zero...'])
end

a_bar(isnan(a_bar)) = 0;

a_reduced = @(z,t) get_a(a_bar, z, t, dt);

%% Create reduced SRN

model_proj = SRN(V_reduced, a_reduced);

end




%%
function a = get_a(a_bar, Z_reduced, t, dt)
idx = num2cell(Z_reduced+1); 
a = a_bar(:, idx{:}, floor(t/dt)+1);
end