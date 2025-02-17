function [pi, t] = fit(filter, Y_obs, t_obs, Z0)
%Solves filtering problem for SRN with exact observations
%
% -------------------------------------------------------------------------
%INPUT
% filter      : object of class FFSP
% Y_obs       : observed process (including Y(0))   | (~,1) array
% t_obs       : jump time points (including 0)      | (~,1) array 
% Z0          : initial state                       | (d,1) array
% -------------------------------------------------------------------------
%OUTPUT
% pi          : conditional PMF for hidden species  | dim(X)+1 array
% t           : time points for pi
% -------------------------------------------------------------------------


nt = length(t_obs);
pi_t = zeros(filter.N, 1);

t = [];
pi = [];

% initial distribution 
Z0(filter.observed_ind) = [];
ind = filter.state2ind(Z0);
pi_t(ind, 1) = 1;

% solve filtering problem
for k = 1:nt-1 
    if filter.show_progress
        disp(['FFSP, model.d = ' num2str(filter.model.d) ...
            ', filter.N = ' num2str(filter.N) ...
            ', Progress: ' num2str((t_obs(k)/t_obs(end)*100), 3) '%' ]);	
    end

    % between jumps (t_k, t_{k+1}):
    [rho, c, t_span] = filter.evolve(pi_t, t_obs(k), t_obs(k+1), ...
                                     Y_obs(:, k));
    pi = [pi, rho ./ c];
    t  = [t, t_span];

    % jump at t_{k+1}:
    dY = Y_obs(:, k+1)-Y_obs(:, k);
    [rho_t, c_t] = filter.jump(rho(:,end), c(end), t_obs(k+1), ...
                               Y_obs(:, k), dY);

    % normalization
    pi_t  = rho_t ./ c_t;
    
end

pi = reshape(pi, [filter.max_X+1; length(t)]');

end




























