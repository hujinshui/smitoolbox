classdef topiclda_em < smi_state
    % Topic Latent Dirichlet allocation E-M algorithm state
    %
    
    % Created by Dahua Lin, on Feb 18, 2012
    %
    
    %% properties
    
    properties(GetAccess='public', SetAccess='private')        
        vsize;          % the size of vocabulary (V)        
        pricount = 0;   % the prior count for topic estimation
        
        estep_iters = 20;   % The (maximum) iterations within each E-step
    end
   
    
    %% observations
    
    properties(GetAccess='public', SetAccess='private')
        ndocs;          % the number of documents (n)
        counts;         % the matrix of word counts [V x n]        
        weights;        % the document weights [empty or 1 x n]
    end
    
    %% dynamic states
    
    properties
        sol;            % the current solution (with following fields)   
                        % K:        the number of topics
                        % alpha:    Dirichlet prior parameters [K x 1]
                        % Beta:     topic-wise word distributions [V x K]
                        % Gamma:    per-document topic-distribution [K x n]
                        % W:        accumulated per-topic word counts [V x K]
        
        objVs;          % a struct with the itemized objective values
                        % ell_theta:    expected log-lik of theta
                        % ell_z:        expected log-lik of z
                        % ell_w:        expected log-lik of words
                        % ent_theta:    entropy of theta
                        % ent_z:        entropy of z
                        % lpri_beta:    log-prior of beta
    end
    
    %% methods
    
    methods
        
        function obj = topiclda_em(V, eta)
            % Construct a topiclda E-M state object
            %
            %   obj = topiclda_em(V);
            %   obj = topiclda_em(V, eta);
            %
            %       V is the size of vocabulary.
            %
            %       Here, eta is the prior-word-count for each topic
            %       (by default, eta = 0)
            %       
            
            if ~(isnumeric(V) && isreal(V) && isscalar(V) && V >= 1 && V == fix(V))
                error('topiclda_em:invalidarg', ...
                    'V should be a positive integer.');
            end          
            obj.vsize = V;
            
            if nargin >= 2
                if ~(isfloat(eta) && isreal(eta) && isscalar(eta) && eta >= 0)
                    error('topiclda_em:invalidarg', ...
                        'eta should be a non-negative scalar.');
                end
                obj.pricount = double(eta);
            end            
        end
        
        
        function obj = initialize(obj, C, w, B0, a0)
            % Initialize the state with observations and initial params
            %
            %   obj = obj.initialize(C, [], K);
            %   obj = obj.initialize(C, w, K);
            %   obj = obj.initialize(C, w, B0);
            %   obj = obj.initialize(C, w, B0, a0);
            %
            %       Inputs:
            %       - C:        The per-document word count matrix [V x n]
            %       - w:        The document weights [empty or 1 x n]
            %       - K:        The number of topics
            %       - B0:       The initial word-distributions [V x K]
            %       - a0:       The initial Dirichlet params [K x 1]
            %
            %       If B0 or a0 is not provided, they are initialized
            %       randomly.
            %
            
            V = obj.vsize;
            
            % check inputs
            
            if ~(isfloat(C) && isreal(C) && ismatrix(C) && size(C,1) == V)
                error('topiclda_em:invalidarg', ...
                    'C should be a real matrix with V rows.');
            end
            n = size(C, 2);
            if ~isa(C, 'double'); C = double(C); end
            
            if ~isempty(w)
                if ~(isfloat(w) && isreal(w) && isvector(w) && numel(w) == n)
                    error('topiclda_em:invalidarg', ...
                        'w should be a real vector of length n.');
                end
                if size(w, 1) > 1; w = w.'; end
                if ~isa(w, 'double'); w = double(w); end
            end
                
            if isscalar(B0)
                K = B0;
                if ~(isnumeric(K) && isreal(K) && isscalar(K) && ...
                        K == fix(K) && K >= 1)
                    error('topiclda_em:invalidarg', ...
                        'K should be a positive integer.');
                end
                K = double(K);
                
                % initialize B0
                B0 = rand(V, K);
                B0 = bsxfun(@times, B0, 1 ./ sum(B0, 1));                
            else
                if ~(isreal(B0) && isfloat(B0) && ismatrix(B0) && size(B0, 1) == V)
                    error('topiclda_em:invalidarg', ...
                        'B0 should be a real matrix with V rows.');
                end
                K = size(B0, 2);
                if issparse(B0)
                    B0 = full(B0);
                end
                if ~isa(B0, 'double'); B0 = double(B0); end
            end
                
            if nargin < 5 || isempty(a0)
                a0 = 1 + rand(K, 1);
            else
                if ~(isreal(a0) && isfloat(a0) && ~issparse(a0) && isequal(size(a0), [K 1]))
                    error('topiclda_em:invalidarg', ...
                        'a0 should be a real vector of size K x 1.');
                end
                if ~isa(a0, 'double'); a0 = double(a0); end
            end
            
            % set fields
                                        
            obj.ndocs = n;
            obj.counts = C;
            if ~isempty(w)
                obj.weights = w;
            end
            
            s.K = K;
            s.alpha = a0;
            s.Beta = B0;
            s.Gamma = [];
            s.W = [];
            
            obj.sol = s;                
        end
        
        
        function obj = update(obj)
            % Updates the state (one E-M iteration)
            %
            %   obj = obj.update();
            %
                                              
            [obj.sol, obj.objVs] = topiclda_em_update( ...
                obj.counts, obj.weights, obj.pricount, obj.sol, ...
                obj.estep_iters, obj.objVs);
        end
        
        
        function R = output(obj)
            % Output the current solution
            %
            %   R = obj.output();
            %
            %       R is a struct with the following fields:
            %       - Beta:     word distributions for topics
            %       - alpha:    Dirichlet parameters
            %       - Gamma:    per-document topic distributions
            %
            
            s = obj.sol;
            
            R.Beta = s.Beta;
            R.alpha = s.alpha;
            R.Gamma = s.Gamma;            
        end
        
        
        function b = is_ready(obj)
            % Tests whether the object is ready for running
            %
            %   b = obj.is_ready();
            %
            
            b = ~isempty(obj.sol);            
        end
        
        
        function objv = evaluate_objv(obj)
            % Evaluate the objective function of the current state
            %
            %   objv = obj.evaluate_objv();
            %
            
            o = obj.objVs;
            objv = o.ell_theta + o.ell_z + o.ell_w + ...
                o.ent_theta + o.ent_z + o.lpri_beta;
        end
    
    end
    
end