function [Y_obs, t_obs, Z, t] = generate_observations(model, Z0, ...
                                                           T, observed_ind)
% Samples one trajectory with SSA, returns observed part & full trajectory
% -------------------------------------------------------------------------
%INPUT
% model         : Stochastic Reaction Network         | object of class SRN
% Z0            : initial state at time 0             | (d,1) array
% T             : final time
% observed_ind  : index/indices of the observed species
% -------------------------------------------------------------------------
%OUTPUT
% Y_obs         : observed part of states at observed jump points
% t_obs         : jump times in observed part
% Z             : all states at all jump points
% t             : all jump times 
% -------------------------------------------------------------------------


[t, Z] = model.sample_exact_path(Z0, T);

Y_obs = Z(observed_ind, :); 
t_obs = t;

% remove points with no jumps in 'Y'
i = 1;
while i < length(t_obs)
    if all(Y_obs(:, i+1) == Y_obs(:, i))
        Y_obs(:, i+1) = [];
        t_obs(i+1) = [];
    else
        i = i + 1;
    end
end

end