classdef SRN
    % Basic implementation of Stochastic Reaction Network (SRN)
    % see examples in 'create_models.m'

    properties 
        d int32,     % # of species
        r int32,     % # of reactions
        V (:,:),     % stoichiometric matrix | (d,r) array 
        a,           % propensities          | lambda function of (z,t)
        is_homogeneous=false,  % true if 'a' are independent of time | bool
    end
    
    methods
        function SRN = SRN(V, a)
            SRN.V = V;
            SRN.a = a;
            [SRN.d, SRN.r] = size(V);
        end
    end
    
    methods (Access=public) 
        % sample with modified next reaction method
        [t_path, Z_path] = sample_exact_path(model, Z0, T)

        % sample with tau-leap
        Z = sample_tl_final(model, Z0, dt, T); 
        Z_path = sample_tl_path(model, Z0, dt, T); 
    end
    
    methods (Access=private)
        Z = tl_step(model, Z0, t, dt)
    end
end