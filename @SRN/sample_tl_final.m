function Z = sample_tl_final(model, Z0, dt, T)
%Tau-leap method for sampling SRN on [0, T], returns only final state
%
% -------------------------------------------------------------------------
%INPUT
% model       : object of class SRN
% Z0          : initial state        | (d,1) array
% dt          : time step
% T           : final time
% -------------------------------------------------------------------------
%OUTPUT
% Z           : state at timet T     | (d,1) array 
% -------------------------------------------------------------------------


Z = Z0;

for t = 0:dt:(T-dt)
    Z = model.tl_step(Z, t, dt); 
end

if t < (T-dt)
    Z = model.tl_step(Z, t, T-t);
end

end