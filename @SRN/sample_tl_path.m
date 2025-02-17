function Z_path = sample_tl_path(model, Z0, dt, T)
%Tau-leap method for SRN sampling, returns states at all time steps
%
% -------------------------------------------------------------------------
%INPUT
% model       : object of class SRN
% Z0          : initial state            | (d,1) array
% dt          : time step
% T           : final time
% -------------------------------------------------------------------------
%OUTPUT
% Z_path      : states at each time step | (d,T/dt) array 
% -------------------------------------------------------------------------


Z_path = zeros(model.d, 1+floor(T/dt));
Z_path(:, 1) = Z0;

i = 1;
for t = 0:dt:(T-dt)
    Z_path(:, i+1) = model.tl_step(Z_path(:, i), t, dt); 
    i = i + 1;
end


if t < (T-dt)
    Z_path(:, end) = model.tl_step(Z_path(:, end-1), t, T-t);
end

end