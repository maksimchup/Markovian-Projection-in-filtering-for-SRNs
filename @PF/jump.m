function filter = jump(filter, t, dY)
%Updates particle states at jump time 
%
% -------------------------------------------------------------------------
%INPUT
% filter      : object of class PF
% t           : time of jump
% dY          : change in the observed agent
% -------------------------------------------------------------------------
%OUTPUT
% filter      : PF with updated states & weights
% -------------------------------------------------------------------------


% split stoichiometric matrix 
V_x = filter.model.V(filter.observed_ind, :);

% suitable reactions for the jump
if isscalar(dY)
    jump_reac = find(V_x == dY);
else
    jump_reac = find(all(V_x == dY));
end

if isempty(jump_reac)
    error(['Incorrect observations. Each increment of Y must be ' ...
        'caused by exactly one reaction'])
end


% update states 
for i = 1:filter.M
    a = filter.model.a(filter.Z(:, i), t);

    j = jump_reac(sample_reaction(length(jump_reac)));

    filter.Z(:, i) = filter.Z(:, i) + filter.model.V(:, j);
    filter.w(i) = filter.w(i) * a(j);
end

end




function j = sample_reaction(L)
%sample random reaction index 

j = 1 + floor(L*rand());

end
