function [rho, c, t_span] = evolve(filter, pi, t, T, Yt)
%Solves ODE system for cond PMF between jumps
% 
% -------------------------------------------------------------------------
%INPUT
% pi        : initial distribution                          | (~,1) array
% t         : initial time 
% T         : final time
% Yt        : observed value on [t, T]                      | int
% -------------------------------------------------------------------------
%OUTPUT
% rho       : unnormalized cond distribution in time(t, T)  | (~,N) array
% c         : normaliztion factor in [t, T]                 | (1,~) array
% t_span    : time grid on [t, T]                           | (1,~) array
% -------------------------------------------------------------------------

t = t+1e-15;

t_span = t:filter.dt:T;
if t_span(end) ~= T
    t_span = [t_span, T];
end

% Solve ode system for 'rho'
if filter.model.is_homogeneous
    A = filter.A(0, Yt);
    odefun = @(t, rho) A*rho;
    %opts = odeset('Jacobian', A);
else
    odefun = @(t, rho) filter.A(t, Yt)*rho;
    %opts = odeset('Jacobian', @(t, rho) filter.A(t, Yt));
end

%[t, rho] = ode15s(odefun, t_span, pi_0, opts); % for stiff systems
[t_span_ode, rho] = ode45(odefun, t_span, pi);
rho = rho';

% Compute normaliztion factor 'c'
if filter.use_conserv
    % FFSP 1 algorithm
    c = sum(rho, 1); 
else
    % FFSP 2 algorithm
    c = evolve_c(filter, rho, t_span_ode, pi, Yt); 
end



% delete time points not in time grid 't_span'
ind = ismember(t_span_ode, t_span);
rho = rho(:, ind);
c = c(ind);


end







function c = evolve_c(filter, rho, t_span_ode, pi, Yt)
%compute upper bound for normalization factor according to FFSP 2 alg
I = zeros(1, size(rho, 2));
for it = 1:length(t_span_ode)-1
    t = t_span_ode(it);
    Xi = zeros(filter.model.d-numel(filter.observed_ind), 1);
    f = 0;
    for i = 1:filter.N
        % full state vector
        Zi = Xi;
        for io = 1:numel(Yt)
            Zi = insert(Zi, filter.observed_ind(io), Yt(io));
        end
            
        ai = filter.model.a(Zi, t);
        for j = filter.unobs_reac_ind
            Zj = Zi + filter.model.V(:,j);
            Xj = remove(Zj, filter.observed_ind);
            if ai(j) > 0 && filter.state2ind(Xj) == 0
                f = f + rho(i, it)*ai(j);
            end
        end

        % state for the next row
        Xi = filter.next_state(Xi);
    end
    I(it+1) = I(it) + (t_span_ode(it+1)-t_span_ode(it))*f;
end
c = 1 + (sum(rho, 1) - sum(pi)) + I;

end



function new_arr = insert(arr, i, value)
%insert value in array (size (~,1))
new_arr = [arr(1:i-1); value; arr(i:end)];
end

function new_arr = remove(arr, i)
%remove element from array (size (~,1))
new_arr = arr;
new_arr(i) = [];
end
