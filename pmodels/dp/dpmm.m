classdef dpmm < handle
    % The class to represent a DP mixture model
    %
    
    % Created by Dahua Lin, on Sep 17, 2011
    %
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='private')   
        np_model;   % the associated non-parametric model
        alpha;      % the concentration parameter
        
        obs;        % the observation
        nobs;       % the number of observations (n) 
        
        natoms;     % the number of active atoms (K)                     
        max_atomid; % the maximum identifier of atom        
        inherits;   % the struct capturing the inherited atoms        
    end
    
    
    properties(GetAccess='private', SetAccess='private')
        atom_ids;   % the (growable) vector of atom identifiers
        atoms;      % the (growable) cell array of atoms
        priweights; % the (growable) vector of prior-weights
        counts;     % the (growable) vector of atom appearance counts 
        logliks0;   % the vector of log-likelihood w.r.t. base
        logliks;    % the (growable) matrix of log-likelihood values
        
        labels;     % the labels of observations (1 x n)
                    % the label is the index of the associated atom
                    % whose id is atom_ids(label)
        groups;     % the cell array of grouped indices
    end
    
    
    methods
        
        function r = get_atom_ids(mdl)
            % Gets the identifiers of all active atoms
            r = mdl.atom_ids(1:mdl.natoms);
        end
        
        function r = get_atoms(mdl)
            % Gets the cell array of all active atoms
            r = mdl.atoms(1:mdl.natoms);
        end
    
        function r = get_atom_counts(mdl)
            % Gets the counts of all active atoms
            r = mdl.counts(1:mdl.natoms);
        end
        
        function r = get_logliks(mdl)
            % Gets the log-likelihoods of all observations w.r.t. all
            % active atoms
            r = mdl.logliks(1:mdl.natoms, :);
        end
    
        function r = get_obs_labels(mdl)
            % Gets the labels (indices of atoms) associated with all 
            % observations
            r = mdl.labels(1:mdl.nobs);
        end
        
        function r = get_obs_atom_ids(mdl)
            % Gets the identifiers of the atoms associated with all
            % observations
            n = mdl.nobs;
            I = mdl.labels(1:n);
            r = zeros(1, n);
            r(I > 0) = mdl.atom_ids(I(I > 0));
        end
        
        function r = get_num_inherits(mdl)
            % Gets the number of inherited atoms            
            h = mdl.inherits;
            if isempty(h)
                r = 0;
            else
                r = h.num;
            end
        end
        
        
        function r = get_capacity(mdl)
            % Gets the capacity of current model
            r = numel(mdl.atoms);
        end       
    end
    
    
    %% Constructor
    
    methods
        
        function mdl = dpmm(model, alpha, obs, varargin)
            % Creates an empty DPMM object
            %
            %   sol = dpmm(model, alpha, obs, ...);
            %
            %       creates an empty solution, based on the input
            %       non-parametric model.
            %
            %       alpha is the concentration parameter, and obs is
            %       the observation set.
            %
            %       One can further specify the following options in form
            %       of name/value pairs
            %
            %       - 'inherits':   a struct specifies the inherited
            %                       atoms, which can be made using
            %                       dpmm_inherits.
            %
            %       - 'reserve_n':  The initial capacity of the model
            %                       should be able to host at least
            %                       reserve_n atoms.
            %
            
            % verify inputs
            
            if ~isa(model, 'nonparam_model')
                error('dpmm:invalidarg', ...
                    'The input model should be of class nonparam_model.');
            end
            
            if ~(isfloat(alpha) && isscalar(alpha) && isreal(alpha) && alpha > 0)
                error('dpmm:invalidarg', ...
                    'alpha should be a positive real number.');
            end                        
            
            H = [];
            rn = [];
            
            if ~isempty(varargin)
                
                onames = varargin(1:2:end);
                ovals = varargin(2:2:end);
                
                if ~(numel(onames) == numel(ovals) && iscellstr(onames))
                    error('dpmm:invalidarg', ...
                        'The option list is invalid.');
                end
                
                for i = 1 : numel(onames)
                    
                    cn = onames{i};
                    cv = ovals{i};
                    
                    switch lower(cn)
                        case 'inherits'
                            if ~(isstruct(cv) && isfield(cv, 'tag') && ...
                                    strcmp(cv.tag, 'dpmm_inherits'))
                                error('dpmm:invalidarg', ...
                                    'The inherits are invalid.');
                            end
                            H = cv;
                        case 'reserve_n'
                            if ~(isnumeric(cv) && isscalar(cv) && cv > 0)
                                error('dpmm:invalidarg', ...
                                    'The reserve_n is invalid.');
                            end
                            rn = cv;
                        otherwise
                            error('dpmm:invalidarg', ...
                                'Unknown option name %s', cn);
                    end
                    
                end                
            end
            
            % set fields
            
            n = model.get_num_samples(obs);            
            
            mdl.np_model = model;
            mdl.alpha = double(alpha);
            
            mdl.obs = obs;
            mdl.nobs = n; 
            
            mdl.natoms = 0;                        
            if isempty(H)
                mdl.max_atomid = 0;
            else
                mdl.max_atomid = H.max_id;
            end
            mdl.inherits = H;
            
            if isempty(rn)
                if isempty(H)
                    rn = 8;
                else
                    rn = max(8, H.num);
                end
            else
                if ~isempty(H)
                    rn = max(rn, H.num);
                end
            end                   
            rK = 2^(ceil(log2(rn)));
            
            mdl.atom_ids = zeros(1, rK);
            mdl.atoms = cell(1, rK);     
            mdl.priweights = zeros(1, rK);
            mdl.counts = zeros(1, rK);
            mdl.logliks0 = model.evaluate_loglik([], obs);
            mdl.logliks = zeros(rK, n);
            mdl.labels = zeros(1, n);
                        
            if ~isempty(H) && H.num > 0
                hn = H.num;
                
                mdl.natoms = hn;
                mdl.atom_ids(1:hn) = H.atom_ids;
                mdl.atoms(1:hn) = H.atoms;
                mdl.priweights(1:hn) = H.pricounts .* H.q;
                
                for i = 1 : H.num
                    a = H.atoms{i};
                    mdl.logliks(i, :) = model.evaluate_loglik(a, obs);
                end
            end
        end
        
    end
       
    
    %% Samping steps
    
    methods
        
        function update_atoms(mdl, ainds)
            % Updates the values of atoms
            %
            %   sol.update_atoms();
            %       Re-draws all atoms conditioned on their associated
            %       observations.
            %
            %   sol.update_atoms(ainds);
            %       Re-draws the atoms selected by ainds (an index vector)
            %       conditioned on their associated observations.
            %   
            
            Kp = get_num_inherits(mdl);
            H = mdl.inherits;
            K = mdl.natoms;           
            g = mdl.groups;
            np_mdl = mdl.np_model;
            X = mdl.obs;
            
            if nargin < 2
                ainds = 1 : K;
            end
            na = numel(ainds);
            
            for j = 1 : na
                k = ainds(j);
                gk = g{k};
                
                if k > Kp
                    a = np_mdl.posterior_atom(X, gk);
                else
                    if ~isempty(gk)
                        a = np_mdl.posterior_atom(X, gk, H.atoms{k});
                        mdl.priweights(k) = H.pricounts(k);
                    else
                        a = H.atoms{k};
                        mdl.priweights(k) = H.pricounts(k) * H.q(k);
                    end
                end
                mdl.atoms{k} = a;
                mdl.logliks(k, :) = np_mdl.evaluate_loglik(a, X);                
            end                            
        end
               
                
        function update_labels(mdl, inds)
            % Update the labels associated with the observations 
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
                
                lpri = log(mdl.priweights(1:K) + mdl.counts(1:K)).';                
                E = bsxfun(@plus, lpri, mdl.logliks(1:K, s));
                ev0 = log(mdl.alpha) + mdl.logliks0(1, s);
                
                z = dpmm_redraw_labels(E, ev0, rand(1, numel(s)));
                
                % update labels to the solution
                                
                mdl.labels(s) = z;
                if first_iter
                    mdl.counts(1:K) = intcount(K, mdl.labels);
                else
                    mdl.counts(1:K) = mdl.counts(1:K) + intcount(K, z);
                end                
                % assert(isequal(get_atom_counts(mdl), intcount(mdl.natoms, mdl.labels)));              
                
                % create new atom if necessary
                
                s = s(z == 0);
                
                if ~isempty(s)
                    mdl.new_atom(s(1));
                    s(1) = [];
                    K = K + 1;
                end
                
                first_iter = false;                
                % assert(isequal(get_atom_counts(mdl), intcount(mdl.natoms, mdl.labels)));
            end
            
            mdl.groups = intgroup(K, mdl.labels);
        end
        
        
        
        function prune_model(mdl, tolratio)
            % Prune the model by removing dead atoms
            %
            %   prune_model(mdl, tolratio);
            %
            
            K = mdl.natoms;
            nh = get_num_inherits(mdl);
            if K <= nh
                return;
            end
            
            is_dead = ((1:K) > nh) & (mdl.counts(1:K) == 0);
            if nnz(is_dead) > K * tolratio
                
                deads = find(is_dead);
                ndeads = numel(deads);
                
                mdl.natoms = K - ndeads;
                mdl.atom_ids(deads) = [];
                mdl.atoms(deads) = [];
                mdl.counts(deads) = [];
                mdl.logliks(deads, :) = [];
                
                if ~isempty(mdl.groups)
                    mdl.groups(deads) = [];
                end
                
                lmap = zeros(1, K);
                lmap(~is_dead) = 1 : (K - ndeads);
                mdl.labels = lmap(mdl.labels);
                
                % assert(isequal(get_atom_counts(mdl), intcount(mdl.natoms, mdl.labels)));
            end
            
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
            X = mdl.obs;
            
            a = np_mdl.posterior_atom(X, i);
            llik = np_mdl.evaluate_loglik(a, X);
            
            if K == get_capacity(mdl)
                new_capa = 2^(ceil(log2(K * 2)));
                mdl.atom_ids(1, new_capa) = 0;
                mdl.atoms{1, new_capa} = [];
                mdl.priweights(1, new_capa) = 0;
                mdl.counts(1, new_capa) = 0;
                mdl.logliks(new_capa, end) = 0;
            end
            
            K = K + 1;                        
            mdl.natoms = K;
            mdl.max_atomid = mdl.max_atomid + 1;
            
            mdl.atom_ids(K) = mdl.max_atomid;
            mdl.atoms{K} = a;            
            mdl.counts(K) = 1;
            mdl.logliks(K, :) = llik;
            
            mdl.labels(i) = K;
        end
    
    end
    
    
end


