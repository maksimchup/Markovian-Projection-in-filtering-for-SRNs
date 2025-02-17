classdef PF
    % Particle filter for SRN for exact and continuous observations
    % ref: https://doi.org/10.1063/5.0032539

    properties
        model,            % object of class SRN
        observed_ind,     % index/indices of the observed species
        M,                % # of particles
        dt=1e-2,          % time step to save solution
        Z,                % current states of particles      | (d,M) array
        w,                % weight of each particle
        observed_reac,    % observed reaction ind
        unobserved_reac,  % unobserved reaction ind
    end

    methods
        function filter = PF(model, observed_ind, M, varargin)
            filter.model = model;
            filter.observed_ind = observed_ind;
            filter.M = M;

            if ~isempty(varargin)
                filter.dt = varargin{1};
            end
            
            V_x = model.V(observed_ind, :);
            filter.observed_reac = find(any(V_x ~= 0, 1));
            filter.unobserved_reac = find(all(V_x == 0, 1));
        end
    end

    methods (Access=public)
        % solves filtering problem
        [X, w, t] = fit(filter, Y, t_obs, X0); 
    end
    
    methods (Access=private)
        % sample particles in between observed jumps 
        [filter, X, w, t] = evolve(filter, t0, T); 

        % update particles at observed jump 
        filter = jump(filter, t, dY); 

        % resample particles
        filter = resample(filter);
    end
end