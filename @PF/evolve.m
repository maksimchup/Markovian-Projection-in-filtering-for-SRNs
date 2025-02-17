function [filter, Z, w, t_span] = evolve(filter, t, T)
%Simulates particle paths between jumps 
% 
% -------------------------------------------------------------------------
%INPUT
% filter      : object of class PF
% t           : initial time
% T           : final time
% -------------------------------------------------------------------------
%OUTPUT
% filter      : PF with updated states & weights
% Z           : states of each particle                   | (d,M,Nt) array
% w           : weights of particles                      | (M,Nt) array
% t_span      : time points (grid with step 'filter.dt')  | (1,Nt) array
% -------------------------------------------------------------------------


t_span = t:filter.dt:T;
if t_span(end) ~= T
    t_span = [t_span, T];
end

Z = zeros(filter.model.d, filter.M, length(t_span));
w = zeros(filter.M, length(t_span));

Z(:, :, 1) = filter.Z;
filter.w = zeros(filter.M, 1); 

for i = 1:filter.M 
    tau = t;
    it = 1;
    
    while tau < T
        % compute propensities
        a = filter.model.a(filter.Z(:,i), tau);
        a_obs = a(filter.observed_reac);
        a_unobs = a(filter.unobserved_reac);
        
        dtau = -log(rand()) / sum(a_unobs);
        dtau = min(dtau, T-tau);
      
        tau = tau + dtau;
        
        filter.w(i) = filter.w(i) - dtau*sum(a_obs);
        
        % 
        while t_span(it) <= tau 
            Z(:, i, it) = filter.Z(:, i);
            w(i, it) = filter.w(i) + (tau-t_span(it))*sum(a_obs);
            
            it = it + 1;
            if it > length(t_span)
                break
            end
        end

        if tau < T
            j = filter.unobserved_reac(sample_reaction(a_unobs));
            filter.Z(:, i) = filter.Z(:, i) + filter.model.V(:, j);
        end
    end

end

w = exp(w);                % all weights at times t_span
filter.w = exp(filter.w);  % weights at time T

end



function j = sample_reaction(a)
%sample random reaction index according to the propensities a
a = a / sum(a);

u = rand();
for j = 1:length(a)
    u = u - a(j);
    if u <= 0
        return
    end
end

end