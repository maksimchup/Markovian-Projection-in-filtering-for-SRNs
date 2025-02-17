function [rho_t, c_t] = jump(filter, rho_t_minus, c_t_minus, t, Y, dY)
%Adjusts the distribution and norm. factor at the jump point 
% 
% -------------------------------------------------------------------------
%INPUT
% rho_t_minus   : unnormalized distribution before jump   | (~,1) array
% c_t_minus     : normalization factor before jump 
% t             : time point of jump
% Y             : observed value before the jump
% dY            : change in the observed agent
% -------------------------------------------------------------------------
%OUTPUT
% rho_t         : corrected unnormalized distribution     | (~,1) array
% c_t           : normaliztion factor for time t          | (~,1) array
% -------------------------------------------------------------------------


% split stoichiometric matrix 
V_y = filter.model.V(filter.observed_ind, :);
V_x = filter.model.V;
V_x(filter.observed_ind, :) = [];


% suitable reactions for the jump
if isscalar(dY)
    reactions_ind = find(dY == V_y);
else
    reactions_ind = find(all(dY == V_y));
end

if isempty(reactions_ind)
    error(['Incorrect observations. Each increment of Y must be ' ...
        'caused by exactly one reaction'])
end

rho_t = zeros(filter.N, 1);
c_t = 0;

Xi = zeros(filter.model.d-numel(filter.observed_ind), 1);
for i = 1:filter.N
    Zi = insert(Xi, filter.observed_ind, Y);
    ai = filter.model.a(Zi, t);
    
    for ir = reactions_ind
        c_t = c_t + ai(ir)*rho_t_minus(i);

        Xj = Xi - V_x(:, ir);
        Zj = insert(Xj, filter.observed_ind, Y);
        j = filter.state2ind(Xj);
        
        aj = filter.model.a(Zj, t);        
        if j ~= 0    
            rho_t(i) = rho_t(i) + aj(ir)*rho_t_minus(j);
        end

    end
    
    % state for the next row
    Xi = filter.next_state(Xi);
end

if ~(filter.use_conserv)
    % we assume that the max propensity is reached in the largest state 
    c_t = c_t + sum(ai(reactions_ind)) * (c_t_minus - sum(rho_t_minus));
end


if sum(rho_t) == 0
    disp('rho = 0 ')
end

end



function new_arr = insert(arr, i, value)
%insert value in array (size (~,1))
if i > length(arr)  % If inserting at the end
    new_arr = [arr; value];
else
    new_arr = [arr(1:i-1); value; arr(i:end)];
end
end
