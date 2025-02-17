classdef FFSP
    % Filtered Finite State Projection for SRN 
    % Solves filtering equation by truncating hidden state space:
    %       X^d  ->  [0, max_X(1)] x ... x [0, max_X(d)]
    %
    % https://doi.org/10.1101/2022.10.18.512737
    
    properties 
        model,            % object of class SRN
        observed_ind,     % index of observed agent        | int
        max_X,            % upper bound for hidden agents  | (d-1,1) array
        dt=1e-2,          % time step to save the solution
        N,                % size of the state space        | int
        unobs_reac_ind,   % unobserved reaction indices
        ind_offsets,      % index offstes (see 'state2ind.m')
        use_conserv=true, % use conservative normalization (FFSP 1 alg) 
        show_progress = false,
    end
    
    methods
        function filter = FFSP(model, observed_ind, max_X, varargin)
            filter.model = model;
            filter.observed_ind = observed_ind;
            filter.max_X = max_X;
            
            filter.N = prod(max_X+1);
            filter.ind_offsets = [1; cumprod(filter.max_X(1:end-1)+1)]';
            
            V_y = model.V(observed_ind, :);
            if isscalar(observed_ind)
                filter.unobs_reac_ind = find(V_y==0);
            else
                filter.unobs_reac_ind = find(all(V_y==0));
            end

            if ~isempty(varargin)
                filter.dt = varargin{1};                
            end
        end
    end
    
    methods (Access=public) 
        % solve filtering problem
        [pi, t] = fit(filter, Y, t_obs, Z0); 
    end
    
    methods (Access=private)
        % solve ode in between observed jumps
        [rho, c, t] = evolve(filter, pi_0, t0, T, Yt); 

        % update distribution at observed jump
        [rho_t, c_t] = jump(filter, rho_t_minus, c_t_minus, t, Yt, dY); 
        
        % get matrix for ode system in between observed jumps
        A = A(filter, t, Yt);
        
        % get index of hiddent state in truncated space
        ind = state2ind(filter, X);

        % get next state in truncated space (lexicographic order)
        X = next_state(filter, X);
    end
end