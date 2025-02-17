function [Z, w, t] = fit(filter, Y_obs, t_obs, Z0)
%Solves filtering problem for SRN with exact observations
%
% -------------------------------------------------------------------------
%INPUT
% filter      : object of class PF
% Y_obs       : observed states at t_obs (including Y(0))   | (~,1) array
% t_obs       : jump time points (including 0)              | (~,1) array 
% Z0          : initial state                               | (d,1) array
% -------------------------------------------------------------------------
%OUTPUT
% Z           : states of each particle
% w           : weights of particles
% t           : time points (grid with step 'filter.dt' between jumps)
% -------------------------------------------------------------------------

nt = length(t_obs);
Z = [];
w = [];
t = [];


filter.Z = repmat(Z0, 1, filter.M);
filter.w = ones(filter.M, 1);

for i = 1:nt-1 % for each interval (t_{k}, t_{k+1}]
    
    [filter, Z_new, w_new, t_new] = filter.evolve(t_obs(i), t_obs(i+1));
    
    filter = filter.jump(t_obs(i+1), Y_obs(:, i+1)-Y_obs(:, i));
    filter = filter.resample();

    % save results
    Z = cat(3, Z, Z_new);
    w = [w, w_new];
    t = [t, t_new];
end

% normalize weights
w = w ./ sum(w, 1);

end