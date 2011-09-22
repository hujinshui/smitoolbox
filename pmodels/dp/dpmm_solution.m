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
        logliks;    % the (growable) matrix of log-likelihood values
        
        labels;     % the labels of observations (1 x n)                
        groups;     % the observation index groups for different labels
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
            if isempty(sol.groups)
                sol.groups = intgroup(K, sol.labels);
            end            
            g = sol.groups;
            mdl = sol.model;
            
            if nargin < 2
                for k = 1 : K
                    a = mdl.create_atom(g{k});
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
                    
                    a = mdl.create_atom(g{kj});
                    sol.atoms{kj} = a;
                    sol.logliks(kj, :) = mdl.evaluate_loglik(a); 
                end
            end
        end
               
                
        function update_labels(sol, i)
            % Performs sequential update 
            %
            %   sol.seq_update_labels();
            %       
            %       Updates the labels of all observations sequentially.
            %       
            %   sol.seq_update_labels(i);
            %
            %       Updates the labels according to the order given by i.
            %
            %   This function will create new atoms.
            %
            
            if isempty(i)
                return;
            end
            
            n = sol.nobs;
            if nargin < 2
                i = int32(1:n);
            else
                i = int32(i);
            end
            ni = numel(i);
            
            rnums = rand(1, ni);
            
            if sol.labels(i(1)) == 0 % draw the first one, if not already
                new_atom(sol, i(1));                
                changed = true;
                cp = 2;
            else                
                changed = false;
                cp = 1;
            end
            
            while cp < ni
                
                % update labels until we have to draw a new one
                [cp, ch] = dpmm_update_labels( ...
                    sol.alpha, ...
                    sol.logliks0, ...
                    sol.logliks, ...
                    [], ...
                    sol.counts, ...
                    sol.labels, ...
                    i, rnums, cp);
                
                changed = changed | ch;
                
                % draw a new one when necessary
                if cp <= ni
                    new_atom(sol, i(cp));
                    changed = true;
                    cp = cp+1;                     
                end
            end
            
            if changed
                sol.groups = [];  % invalidate the old grouping                
            end
        end
        
    end
    
    
    %% Private implementation
    
    methods(Access='private')
    
        function a = new_atom(sol, I)
            
            mdl = obj.model;
            a = mdl.create_atom(I);
            
            K = sol.natoms;
            capa = numel(sol.atoms);
            
            if K == capa  % need to grow
                sol.atoms{1, 2 * K} = [];
                sol.counts(1, 2 * K) = 0;
                sol.logliks(2 * K, end) = 0;
            end
            
            K = K + 1;
            sol.atoms{K} = a;
            sol.counts(K) = numel(I);
            sol.natoms = K;
            sol.logliks(K,:) = mdl.evaluate_loglik(a); 
            
            sol.labels(I) = K;            
        end        
    
    end
    
    
end
