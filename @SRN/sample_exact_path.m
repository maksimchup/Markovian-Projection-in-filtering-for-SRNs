function [t_path, Z_path] = sample_exact_path(model, Z0, T)
%Modified next reaction method for time homogeneous(!) SRN on [0, T], 
% returns (time, state) at jump times
%
% -------------------------------------------------------------------------
%INPUT
% model       : object of class SRN
% Z0          : initial state                   | (d,1) array
% T           : final time
% -------------------------------------------------------------------------
%OUTPUT
% t_path      : jump times (including 0 and T)  | (~,1) array
% Z_path      : stats at jump times             | (d,~) array 
% -------------------------------------------------------------------------

if ~(model.is_homogeneous)
    warning('Exact sampling is not implemented for non-homogeneous SRNs')
    % see part 5 in [Anderson2007], https://doi.org/10.1063/1.2799998
end

% preallocate memory for Nt steps
Nt = 500;
t_path = zeros(1, Nt);
Z_path = zeros(model.d, Nt);


Tau = zeros(model.r, 1);
P = log(1./rand(model.r, 1));

Z = Z0;
t = 0;
i = 1;
while t < T
    t_path(i) = t;
    Z_path(:, i) = Z;

    a = model.a(Z, t);
    dt = (P - Tau) ./ a;
    [delta, mu] = min(dt);  % find next reaction
    
    Tau = Tau + a.*delta;   % update internal times
    P(mu) = P(mu) + log(1/rand());
    
    t = t + delta;          % update time
    Z = Z + model.V(:, mu); % update state
    i = i + 1;
end

% save state at final time T
t_path(i) = T;
Z_path(:, i) = Z_path(:, i-1);

% delete excess preallocated memory
t_path = t_path(1:i);
Z_path = Z_path(:, 1:i);

end
