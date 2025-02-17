function Z = tl_step(model, Z0, t, dt)
%One step of Tau-leap method for sampling SRN
%
% -------------------------------------------------------------------------
%INPUT
% model       : object of class SRN
% Z0          : initial state at time    | (d,1) array
% t           : current time
% dt          : time step
% -------------------------------------------------------------------------
%OUTPUT
% Z           : final state at time t+dt | (d,1) array
% -------------------------------------------------------------------------


a0 = max(model.a(Z0, t), 0); 

dZ = model.V * poissrnd(dt * a0);

Z = Z0 + dZ;
Z = max(Z, 0);

end