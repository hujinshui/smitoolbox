classdef dpmm < handle
    % The class to represent a DP mixture model
    %
    
    % Created by Dahua Lin, on Sep 17, 2011
    %
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='private')   
        np_model;   % the associated non-parametric model
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
        
        function mdl = dpmm(model, alpha, rK)
            % Creates an empty DPMM object
            %
            %   sol = dpmm(model, alpha);
            %
            %       creates an empty solution, based on the input
            %       non-parametric model.
            %
            %       alpha is the concentration parameter.
            %
            %   sol = dpmm(model, alpha, rK)
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
            
            mdl.np_model = model;
            mdl.alpha = double(alpha);
            mdl.natoms = 0;
            mdl.nobs = n;
            
            mdl.atoms = cell(1, rK);
            mdl.counts = zeros(1, rK);
            mdl.logliks0 = model.evaluate_loglik();
            mdl.logliks = zeros(rK, n);
            mdl.labels = zeros(1, n);
        end
        
        
        function update_atoms(mdl, k)
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
            
            K = mdl.natoms;           
            g = mdl.groups;
            np_mdl = mdl.np_model;
            
            if nargin < 2
                for k = 1 : K
                    a = np_mdl.posterior_atom(g{k});
                    mdl.atoms{k} = a;
                    mdl.logliks(k, :) = np_mdl.evaluate_loglik(a);
                end
            else
                for j = 1 : numel(k)
                    kj = k(j);
                    if ~(kj >= 1 && kj <= K)
                        error('dpmm_solution:invalidarg', ...
                            'The atom index exceeds valid range.');
                    end
                    
                    a = np_mdl.posterior_atom(g{kj});
                    mdl.atoms{kj} = a;
                    mdl.logliks(kj, :) = np_mdl.evaluate_loglik(a); 
                end
            end
        end
               
                
        function update_labels(mdl, inds)
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
                ni = mdl.nobs;
                s = 1:ni;
            else
                s = inds;
            end
            
            if mdl.natoms == 0
                mdl.new_atom(s(1));
                s(1) = [];
            end
            
            first_iter = true;
            
            while ~isempty(s)
                
                % re-draw labels
                
                K = mdl.natoms;
                E = bsxfun(@plus, log(mdl.counts(1:K)).', mdl.logliks(1:K, s));
                ev0 = log(mdl.alpha) + mdl.logliks0(1, s);
                
                z = dpmm_redraw_labels(E, ev0, rand(1, numel(s)));
                
                % update labels to the solution
                
                mdl.labels(s) = z;
                if first_iter
                    mdl.counts(1:K) = intcount(K, mdl.labels);
                else
                    mdl.counts(1:K) = mdl.counts(1:K) + intcount(K, z);
                end
                
                % assert(isequal(get_atom_counts(sol), intcount(sol.natoms, sol.labels)));
                
                % create new atom if necessary
                
                s = s(z == 0);
                
                if ~isempty(s)
                    mdl.new_atom(s(1));
                    s(1) = [];
                    K = K + 1;
                end
                
                first_iter = false;
                
                % assert(isequal(get_atom_counts(sol), intcount(sol.natoms, sol.labels)));
            end
            
            mdl.groups = intgroup(K, mdl.labels);
        end
        
    end
    
    
    methods(Access='private')
        
        function new_atom(mdl, i)
            % Creates a new atom based on the i-th observation
            %
            % (pre-condition, the i-th obs has no label)
            %
            
            assert(mdl.labels(i) == 0);
            
            K = mdl.natoms;
            np_mdl = mdl.np_model;
            
            a = np_mdl.posterior_atom(i);
            llik = np_mdl.evaluate_loglik(a);
            
            if K == get_capacity(mdl)
                new_capa = 2 * K;
                mdl.atoms{1, new_capa} = [];
                mdl.counts(1, new_capa) = 0;
                mdl.logliks(new_capa, end) = 0;
            end
            
            K = K + 1;
            mdl.natoms = K;
            mdl.atoms{K} = a;
            mdl.counts(K) = 1;
            mdl.logliks(K, :) = llik;
            
            mdl.labels(i) = K;
        end
                
    end
    
    
end


