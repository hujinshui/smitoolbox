classdef dpmm_solution < handle
    % The class to represent a solution of DPMM inference
    %
    
    % Created by Dahua Lin, on Sep 17, 2011
    %
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='private')   
        model;      % th associated non-parametric model
        alpha;      % the concentration parameter
        
        natoms;     % the number of atoms (K)        
        nobs;       % the number of observations (n)        
    end
    
    
    properties(GetAccess='private', SetAccess='private')                        
        atoms;      % the (growable) cell array of atoms
        counts;     % the (growable) vector of atom appearance counts 
        logliks0;   % the vector of log-likelihood w.r.t. base
        logliks;    % the (growable) matrix of log-likelihood values
        
        labels;     % the labels of observations (1 x n)
        groups;     % the cell array of grouped indices
    end
    
    
    methods
        
        function r = get_atoms(sol)
            r = sol.atoms(1:sol.natoms);
        end
    
        function r = get_atom_counts(sol)
            r = sol.counts(1:sol.natoms);
        end
        
        function r = get_logliks(sol)
            r = sol.logliks(1:sol.natoms, :);
        end
    
        function r = get_labels(sol)
            r = sol.labels(1:sol.nobs);
        end
        
        function r = get_capacity(sol)
            r = numel(sol.atoms);
        end       
    end
    
    
    %% Methods
    
    methods
        
        function sol = dpmm_solution(model, alpha, rK)
            % Creates an empty DPMM soltion object
            %
            %   sol = dpmm_solution(model, alpha);
            %
            %       creates an empty solution, based on the input
            %       non-parametric model.
            %
            %       alpha is the concentration parameter.
            %
            %   sol = dpmm_solution(model, alpha, rK)
            %
            %       Creates an empty solution such that using its initial
            %       capacity, it can host at least rK atoms.
            %
            
            % verify inputs
            
            if ~isa(model, 'nonparam_model')
                error('dpmm_solution:invalidarg', ...
                    'The input model should be of class nonparam_model.');
            end
            
            if ~(isfloat(alpha) && isscalar(alpha) && isreal(alpha) && alpha > 0)
                error('dpmm_solution:invalidarg', ...
                    'alpha should be a positive real number.');
            end
            
            if nargin < 3
                rK = 8;
            else
                if ~(isnumeric(rK) && isscalar(rK) && rK > 0)
                    error('dpmm_solution:invalidarg', ...
                        'rK should be a positive scalar.');
                end
                rK = 2^(ceil(log2(double(rK))));
            end
            
            % set fields
            
            n = model.nobs;
            
            sol.model = model;
            sol.alpha = double(alpha);
            sol.natoms = 0;
            sol.nobs = n;
            
            sol.atoms = cell(1, rK);
            sol.counts = zeros(1, rK);
            sol.logliks0 = model.evaluate_loglik();
            sol.logliks = zeros(rK, n);
            sol.labels = zeros(1, n);
        end
        
        
        function update_atoms(sol, k)
            % Updates the values of atoms
            %
            %   sol.update_atoms();
            %       Re-draws all atoms conditioned on their associated
            %       observations.
            %
            %   sol.update_atoms(k);
            %       Re-draws the atoms selected by k (an index vector)
            %       conditioned on their associated observations.
            %   
            
            K = sol.natoms;           
            g = sol.groups;
            mdl = sol.model;
            
            if nargin < 2
                for k = 1 : K
                    a = mdl.posterior_atom(g{k});
                    sol.atoms{k} = a;
                    sol.logliks(k, :) = mdl.evaluate_loglik(a);
                end
            else
                for j = 1 : numel(k)
                    kj = k(j);
                    if ~(kj >= 1 && kj <= K)
                        error('dpmm_solution:invalidarg', ...
                            'The atom index exceeds valid range.');
                    end
                    
                    a = mdl.posterior_atom(g{kj});
                    sol.atoms{kj} = a;
                    sol.logliks(kj, :) = mdl.evaluate_loglik(a); 
                end
            end
        end
               
                
        function update_labels(sol, inds)
            % Performs sequential update 
            %
            %   sol.update_labels();
            %       
            %       Updates the labels of all observations.
            %
            %   sol.update_labels(inds);
            %       
            %       Updates the labels associated with the observations 
            %       selected by inds.
            %
            %       Note that inds should not contain repeated indices
            %
            %
            %   This function will create new atoms.
            %
                                    
            if nargin < 2
                ni = sol.nobs;
                s = 1:ni;
            else
                s = inds;
            end
            
            if sol.natoms == 0
                sol.new_atom(s(1));
                s(1) = [];
            end
            
            first_iter = true;
            
            while ~isempty(s)
                
                % re-draw labels
                
                K = sol.natoms;
                E = bsxfun(@plus, log(sol.counts(1:K)).', sol.logliks(1:K, s));
                ev0 = log(sol.alpha) + sol.logliks0(1, s);
                
                z = dpmm_redraw_labels(E, ev0, rand(1, numel(s)));
                
                % update labels to the solution
                
                sol.labels(s) = z;
                if first_iter
                    sol.counts(1:K) = intcount(K, sol.labels);
                else
                    sol.counts(1:K) = sol.counts(1:K) + intcount(K, z);
                end
                
                % assert(isequal(get_atom_counts(sol), intcount(sol.natoms, sol.labels)));
                
                % create new atom if necessary
                
                s = s(z == 0);
                
                if ~isempty(s)
                    sol.new_atom(s(1));
                    s(1) = [];
                    K = K + 1;
                end
                
                first_iter = false;
                
                % assert(isequal(get_atom_counts(sol), intcount(sol.natoms, sol.labels)));
            end
            
            sol.groups = intgroup(K, sol.labels);
        end
        
    end
    
    
    methods(Access='private')
        
        function new_atom(sol, i)
            % Creates a new atom based on the i-th observation
            %
            % (pre-condition, the i-th obs has no label)
            %
            
            assert(sol.labels(i) == 0);
            
            K = sol.natoms;
            mdl = sol.model;
            
            a = mdl.posterior_atom(i);
            llik = mdl.evaluate_loglik(a);
            
            if K == get_capacity(sol)
                new_capa = 2 * K;
                sol.atoms{1, new_capa} = [];
                sol.counts(1, new_capa) = 0;
                sol.logliks(new_capa, end) = 0;
            end
            
            K = K + 1;
            sol.natoms = K;
            sol.atoms{K} = a;
            sol.counts(K) = 1;
            sol.logliks(K, :) = llik;
            
            sol.labels(i) = K;
        end
                
    end
    
    
end


