classdef dpmm_sol
    % The class that captures a DPMM solution
    %
    
    % Created by Dahua Lin, on Sep 26, 2011
    %
    
    
    %% properties
    
    properties(GetAccess='public', SetAccess='private')
        nobs;       % the number of observations
        
        iatoms;     % the vector of the atom-indices associated with the observations
        
        groups;     % the cell array of grouped obs indices based on labels
    end
    
    properties(GetAccess='private', SetAccess='private')
        asys;       % an internal struct maintaining a set of atoms
    end
    
    
    properties(Dependent)
        natoms;         % the number of atoms (K)
        max_atom_id;    % the maximum identifier of the atoms
        capacity;       % the maximum #atoms that can be accomodated without growing
        
        atom_ids;       % the vector of atom identifiers (1 x K)
        atoms;          % the cell array of atoms (1 x K)
        priweights;     % the vector of prior weights (1 x K)
        atom_counts;    % the current counts of atoms (1 x K)
        logliks;        % the matrix of log-likelihood values (K x n)
        base_logliks;   % the vector of log-lik values w.r.t. the base
        
        labels;         % the vector of labels (ids of atoms associated with obs)
    end
    
    methods
        function v = get.natoms(sol)
            v = sol.asys.natoms;
        end
        
        function v = get.max_atom_id(sol)
            v = sol.asys.max_aid;
        end
        
        function v = get.capacity(sol)
            v = sol.asys.capacity;
        end
        
        function v = get.atom_ids(sol)
            v = sol.asys.atom_ids(1:sol.natoms);
        end
        
        function v = get.atoms(sol)
            v = sol.asys.atoms(1:sol.natoms);
        end
        
        function v = get.priweights(sol)
            v = sol.asys.priweights(1:sol.natoms);
        end
        
        function v = get.atom_counts(sol)
            v = sol.asys.counts(1:sol.natoms);
        end
        
        function v = get.logliks(sol)
            v = sol.asys.logliks(1:sol.natoms, :);
        end
        
        function v = get.base_logliks(sol)
            v = sol.asys.logliks0;
        end
        
        function v = get.labels(sol)
            I = sol.iatoms;
            n = sol.nobs;
            v = zeros(1, n);
            v(I > 0) = sol.atom_ids(I(I > 0));
        end
    end
    
    %% Construction
    
    methods
        
        function sol = dpmm_sol(amdl, obs, base_id, r0)
            % Constructs an empty DPMM solution
            %
            %   sol = dpmm_sol(amdl, obs);
            %   sol = dpmm_sol(amdl, obs, base_id);
            %   sol = dpmm_sol(amdl, obs, base_id, r0);
            %
            %       Constructs an empty DPMM solution.
            %
            %       amdl is the underlying atom model, and obs is
            %       the observation array.
            %
            %       base_id is the base from which the atom id is
            %       incremented. default = 0.
            %
            %       If r0 is specified, the solution has an initial
            %       capacity that is enough to host at least r0 atoms.
            %
            
            if ~isa(amdl, 'atom_model_base')
                error('dpmm_sol:invalidarg', ...
                    'amdl should be an instance of a class derived from atom_model_base.');
            end
            
            if nargin >= 3
                if ~(isnumeric(base_id) && isscalar(base_id) && ...
                        base_id >= 0 && base_id == fix(base_id))
                    error('dpmm_sol:invalidarg', ...
                        'max_id should be a non-negative integer scalar.');
                end
            else
                base_id = 0;
            end
            
            if nargin >= 4
                if ~(isnumeric(r0) && isscalar(r0) && r0 >= 0 && r0 == fix(r0))
                    error('dpmm_sol:invalidarg', ...
                        'r0 should be a non-negative integer scalar.');
                end
                c0 = 2^(ceil(log2(r0)));
            else
                c0 = 8;
            end
            
            n = amdl.get_num_samples(obs);
            
            sol.nobs = n;
            sol.iatoms = zeros(1, n);
            sol.groups = [];
            
            llik0 = amdl.evaluate_loglik([], obs);
            sol.asys = dpmm_sol.init_atomsys(c0, n, base_id, llik0);
            
        end
        
    end
    
    
    %% Methods for updating
    
    methods
        
        function sol = update_atoms(sol, amdl, obs, H, ainds)
            % Updating atoms based on grouped observations
            %
            %   sol = sol.update_atoms(amdl, obs);
            %   sol = sol.update_atoms(amdl, obs, H);
            %   sol = sol.update_atoms(amdl, obs, H, ainds);
            %
            %       updates all (or selected) atoms based on the current
            %       labeling.
            %
            %       Input arguments:
            %       - amdl:     the underlying atom model
            %       - obs:      the observation array
            %       - H:        the struct capturing the inherited atoms
            %                   (it can be omitted when there is no
            %                   inheritance)
            %       - ainds:    the selected indices to atoms to be
            %                   updated. (update all atoms if omitted).
            %
            
            % parse input arguments
            
            if nargin < 4
                H = [];
                Kp = 0;
            else
                if isempty(H)
                    Kp = 0;
                else
                    Kp = H.num;
                end
            end
            
            AS = sol.asys;
            
            if nargin < 5
                ainds = 1 : AS.natoms;
            end
            na = numel(ainds);
            
            if na == 0
                return;
            end
            
            % do updates
            
            grps = sol.groups;
            if isempty(grps)
                error('dpmm_sol:rterror', ...
                    'The observations have to been grouped before doing atom updates.');
            end
            
            for j = 1 : na
                k = ainds(j);
                gk = grps{k};
                
                if k > Kp
                    a = amdl.posterior_atom(obs, gk);
                else
                    if ~isempty(gk)
                        a = amdl.posterior_atom(obs, gk, H.atoms{k});
                        AS.priweights(k) = H.pricounts(k);
                    else
                        a = H.atoms{k};
                        AS.priweights(k) = H.pricounts(k) * H.q(k);
                    end
                end
                AS.atoms{k} = a;
                AS.logliks(k, :) = amdl.evaluate_loglik(a, obs);
            end
            
            sol.asys = AS;
        end
        
        
        function sol = update_labels(sol, amdl, obs, alpha, inds)
            % Update the labels associated with the observations
            %
            %   sol = sol.update_labels(amdl, obs, alpha);
            %   sol = sol.update_labels(amdl, obs, alpha, inds);
            %
            %       Updates the labels of all or selected observations.
            %
            %       Input arguments:
            %       - amdl:     the underlying atom model
            %       - obs:      the observation array
            %       - alpha:    the concentration parameter
            %       - inds:     the indices of the observations selected
            %                   to be updated (if omitted, all labels
            %                   are to be updated)
            %
            %       Note:
            %       - inds should not contain repeated indices.
            %       - the method will create new atoms when necessary.
            %
            
            if nargin < 5
                ni = sol.nobs;
                s = 1:ni;
            else
                s = inds;
            end
            
            AS = sol.asys;
            
            if sol.natoms == 0
                a = amdl.posterior_atom(obs, s(1));
                llik = amdl.evaluate_loglik(a, obs);
                AS = dpmm_sol.add_atom(AS, a, llik);
                
                sol.iatoms(s(1)) = 1;
                AS.counts(1) = 1;
                s(1) = [];
            end
            
            niters = 0;
            lalpha = log(double(alpha));
            
            while ~isempty(s)
                
                niters = niters + 1;
                
                % re-draw labels
                
                K = AS.natoms;
                
                lpri = log(AS.priweights(1:K) + AS.counts(1:K)).';
                E = bsxfun(@plus, lpri, AS.logliks(1:K, s));
                ev0 = lalpha + AS.logliks0(1, s);
                
                z = dpmm_redraw_labels(E, ev0, rand(1, numel(s)));
                
                % update labels to the solution
                
                sol.iatoms(s) = z;
                if niters == 1
                    AS.counts(1:K) = intcount(K, sol.iatoms);
                else
                    AS.counts(1:K) = AS.counts(1:K) + intcount(K, z);
                end
                
                % create new atom if necessary
                
                s = s(z == 0);
                
                if ~isempty(s)
                    a = amdl.posterior_atom(obs, s(1));
                    llik = amdl.evaluate_loglik(a, obs);
                    AS = dpmm_sol.add_atom(AS, a, llik);
                    
                    K = AS.natoms;
                    sol.iatoms(s(1)) = K;
                    AS.counts(K) = 1;
                    s(1) = [];
                end
                
                
            end
            
            sol.asys = AS;
            sol.groups = intgroup(K, sol.iatoms);
            
            % assert(dpmm_sol.verify_sol(sol.asys, sol.iatoms));
        end
        
        
        function sol = prune(sol, Kp, tolratio)
            % Prune the solution by removing dead atoms
            %
            %   sol = prune(sol, Kp, tolratio);
            %
            %       Prunes the solution by removing dead atoms from
            %       the solution, when the ratio of dead atoms is
            %       above the tolratio.
            %
            %       Kp is the number of atoms inherited from the prior,
            %       which should not be removed, even if they are not
            %       seen in the current samples.
            %
            
            K = sol.natoms;
            
            is_dead = sol.atom_counts == 0;
            if Kp > 0
                is_dead(1:Kp) = false;
            end
            
            if nnz(is_dead) > K * tolratio
                
                deads = find(is_dead);
                [sol.asys, rlmap] = dpmm_sol.prune_atoms(sol.asys, deads);
                
                if ~isempty(sol.groups)
                    sol.groups(deads) = [];
                end
                sol.iatoms = rlmap(sol.iatoms);
                
                % assert(dpmm_sol.verify_sol(sol.asys, sol.iatoms));
            end
        end
        
    end
    
    
    %% private implementation
    
    methods(Static, Access='private')
        
        function S = init_atomsys(K0, n, max_id, logliks0)
            % Initialize an internal atom system
            
            S.natoms = 0;
            S.max_aid = max_id;
            S.capacity = K0;
            S.atom_ids = zeros(1, K0);
            S.atoms = cell(1, K0);
            
            S.priweights = zeros(1, K0);
            S.counts = zeros(1, K0);
            S.logliks = zeros(K0, n);
            S.logliks0 = logliks0;
            
        end
        
        
        function S = add_atom(S, a, llik)
            % Adds a new atom
            %
            
            K = S.natoms;
            
            if K == S.capacity  % grow the capacity
                new_capa = 2^(ceil(log2(K * 2)));
                
                S.capacity = new_capa;
                S.atom_ids(1, new_capa) = 0;
                S.atoms{1, new_capa} = [];
                S.priweights(1, new_capa) = 0;
                S.counts(1, new_capa) = 0;
                S.logliks(new_capa, end) = 0;
            end
            
            K = K + 1;
            
            S.natoms = K;
            S.max_aid = S.max_aid + 1;
            S.atom_ids(K) = S.max_aid;
            S.atoms{K} = a;
            S.logliks(K, :) = llik;
        end
        
        
        function [S, rlmap] = prune_atoms(S, ainds)
            % Prune the atom system by removing specified atoms
            %
            
            nd = numel(ainds);  % the number of atoms to be deleted
            
            K = S.natoms;
            S.natoms = K - nd;
            S.atom_ids(ainds) = [];
            S.atoms(ainds) = [];
            S.priweights(ainds) = [];
            S.counts(ainds) = [];
            S.logliks(ainds, :) = [];
            
            S.capacity = numel(S.atoms);
            
            % relabeling map
            is_retained = true(1, K);
            is_retained(ainds) = false;
            rlmap = zeros(1, K);
            rlmap(is_retained) = 1 : (K - nd);
        end
        
        
        function tf = verify_sol(AS, z)
            % A function for verifying the validity of the solution (for DEBUG)
            %
            
            K = AS.natoms;
            c = AS.counts(1:K);
            cr = intcount(K, z);
            tf = isequal(c, cr);
        end
        
    end
    
end

