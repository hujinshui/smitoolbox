classdef gmm_std < smi_prg
    % Standard Gaussian mixture model (GMM) program
    %
    %   The GMM program is an SMI program that implements the inference
    %   over the following formulation.
    %
    %       u_k ~ Gauss(mu, Cu);        for k = 1, ..., K
    %       x_i ~ Gauss(u_{z_i}, Cx);   for i = 1, ..., N
    %
    %   Cx is either given, or from an inverse Wishart distribution.
    %
    
    %   History
    %   -------
    %       - Created by Dahua Lin, on Aug 31, 2011
    %       - Modified by Dahua Lin, on Sep 3, 2011
    %           - based on new smi_prg base.
    %
    
    %% Properties
    
    properties(GetAccess='public', SetAccess='private')
        
        dim;        % the dimension of the underlying space
        
        gmodel;     % the Gaussian generative model
        
        Cu;         % The prior covariance of u (pdmat)
        Cx;         % Cx itself or its prior distribution
        
        est_Cx;     % whether to estimate Cx  
    end
    
    
    %% Construction
    
    methods
        
        function prg = gmm_std(d, Cx, Cu)
            % Creates a standard Gaussian mixture model (GMM) program
            %                        
            %   prg = gmm_std(d, [], Cx);
            %   prg = gmm_std(d, upri, Cx);
            %
            %       Construsts an SMI program that implements the inference of
            %       GMM parameters.
            %
            %       Input arguments:
            %       - d:        the underlying space dimension
            %       - upri:     the Gaussian prior of the variable u.
            %       - Cx:       the covariance for generating x.
            %       
                        
            % verify input arguments
            
            if ~(isnumeric(d) && isscalar(d) && d == fix(d) && d > 0)
                error('gmm_std:invalidarg', 'd should be a positive integer.');
            end
                        
            if nargin >= 2 && ~isempty(Cx)
                if is_pdmat(Cx)
                    if ~(Cx.n == 1 && Cx.d == d)
                        error('gmm_std:invalidarg', ...
                            'Cx should have Cx.n == 1 and Cx.d == d.');
                    end
                else
                    error('gmm_std:invalidarg', ...
                        'The 2nd arg (Cx or Cx_pri) is invalid.');
                end
            else
                Cx = [];
            end
            
            if nargin >= 3 && ~isempty(Cu)
                if ~(is_pdmat(Cu) && Cu.n == 1 && Cu.d == d)
                    error('gmm_std:invalidarg', ...
                        'Cu is not a valid pdmat struct.');
                end
            else
                Cu = [];
            end                         
            
            % create object
            
            prg.dim = double(d);
            
            if ~isempty(Cx)
                if ~isempty(Cu)
                    prg.gmodel = gaussgm(Cx, Cu);
                else
                    prg.gmodel = gaussgm(Cx);
                end
            else
                if ~isempty(Cu)
                    prg.gmodel = gaussgm(d, Cu);
                else
                    prg.gmodel = gaussgm(d);
                end
            end
            
            prg.Cu = Cu;
            prg.Cx = Cx;
            prg.est_Cx = isempty(Cx) || ~is_pdmat(Cx);  
            
            if prg.est_Cx
                error('gmm_std:notimplement', ...
                    'Estimation of Cx is not yet implemented.');
            end
        end
        
    end
    
    
    %% Implementation
    
    methods
    
        function [Sd, Sc] = initialize(prg, obs, L0, optype)
            % Initialize the states
            %
            %   [Sd, Sc] = prg.initialize(obs, L0, 'sample');
            %   [Sd, Sc] = prg.initialize(obs, L0, 'varinfer');
            %
            %       Input arguments:
            %       - obs:  a struct with the following fields:
            %               - X:    the sample matrix (d x n),
            %                       each column is a sample
            %               - K:    the number of mixture components
            %               - mu:   the prior mean (0 or d x 1 vector)
            %                       (if no such field or mu is empty,
            %                        it means no prior for u is used)
            %
            %       - L0:   It can be given in three forms:
            %               - empty:  use default method to initialize
            %               - 1 x n label vector
            %               - K x n label distribution matrix            
            %
            %       Note that these two arguments are packed in a 
            %       1 x 2 cell array as a single input.
            %
            %       Output arguments:
            %       - Sd:   dynamic state struct, with fields: 
            %               - U:    d x K matrix of component means
            %               - Cx:   component covariance(s) in pdmat
            %               - Lik:  K x n pairwise likelihood table
            %
            %       If optype is 'sample', Sd also has
            %               - Z:    1 x n label vector
            %               - grps: 1 x K cell array of grouped indices
            %
            %       If optype is 'varinfer', Sd also has
            %               - Q:    K x n matrix: posterior label
            %                       distributions.
            %               
            %       - Sc:   fixed state struct, with fields:
            %               - X:        the observation matrix
            %               - n:        the number of samples in X
            %               - K:        the number of mixture components
            %               - mu:       the prior mean
            %               - do_samp:  whether it is doing sampling
            %
            
            % verify inputs
            
            d = prg.dim;
            
            if ~(isstruct(obs) && isscalar(obs) && ...
                    isfield(obs, 'X') && isfield(obs, 'K'))
                error('gmm_std:invalidarg', 'obs is not a valid struct.');
            end
            
            X = obs.X;
            K = obs.K;
            mu = [];
            
            if ~(isfloat(X) && ndims(X) == 2 && isreal(X) && size(X,1) == d)
                error('gmm_std:invalidarg', 'X should be a real matrix.');
            end
            n = size(X, 2);
            
            if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 1)
                error('gmm_std:invalidarg', 'K should be an integer.');
            end
            K = double(K);
            
            if isfield(obs, 'mu') && ~isempty(obs.mu)
                mu = obs.mu;
                
                if ~(isfloat(mu) && isreal(mu) && ...
                        (isequal(mu, 0) || isequal(size(mu), [d 1])))
                    error('gmm_std:invalidarg', ...
                        'mu should be either 0 or a d x 1 real vector.');
                end
                
                if isequal(mu, 0)
                    mu = zeros(d, 1);
                end
            end                        
            
            if ~isempty(L0)
                if ~(isnumeric(L0) && ndims(L0) == 2 && ...
                        size(L0,2) == n && (size(L0,1) == 1 || size(L0,1) == K));
                    error('gmm_std:invalidarg', ...
                        'L0 should be a 1 x n row vector.');                
                end
                
                if size(L0, 1) == 1
                    if ~all(L0 == fix(L0) && L0 >= 1 && L0 <= K)
                        error('gmm_std:invalidarg', ...
                            'Some labels in L0 is invalid.');
                    end
                end                                        
            end
                            
            switch optype
                case 'sample'
                    do_samp = true;
                case 'varinfer'
                    do_samp = false;
                otherwise
                    error('gmm_std:invalidarg', ...
                        'Unsupported operation type %s', optype);
            end
                    
            % create states
            
            Sc.X = X;
            Sc.n = size(X, 2);
            Sc.K = K;
            Sc.mu = mu;
            Sc.do_samp = do_samp;            
            
            Sd.U = [];
            if prg.est_Cx
                Sd.Cx = [];
            else
                Sd.Cx = prg.Cx;
            end
            Sd.Lik = [];
                        
            if do_samp
                if isempty(L0)
                    Sd.Z = rand_label(K, n, 'v');
                elseif size(L0, 1) == 1
                    Sd.Z = L0;
                else
                    [~, Sd.Z] = max(L0, [], 1);
                end
                Sd.grps = intgroup(K, Sd.Z);
            else
                if isempty(L0)
                    Sd.Q = rand_label(K, n, 'b');
                elseif size(L0, 1) == 1
                    Sd.Q = l2mat(K, L0);
                else
                    Sd.Q = L0;
                end
            end                        
            
        end    
        
        
        function [Sd, Sc] = update(prg, Sd, Sc)
            % Updates the states (one iteration)
            %
            %   [Sd, Sc] = prg.update(Sd, Sc);
            %
            
            % estimate model parameters (M-step)
            
            gm = prg.gmodel;
            X = Sc.X;
            K = Sc.K;
            mu = Sc.mu;
            
            if Sc.do_samp
                Sd.U = gm.sample_u(X, Sd.Cx, mu, Sd.grps);
            else
                Sd.U = gm.mapest_u(X, Sd.Cx, mu, Sd.Q);
            end
            
            % infer labels (E-step)
            
            Lik = gm.loglik(Sd.U, X);            
            Sd.Lik = Lik;
            
            E = Lik;
            
            if Sc.do_samp
                [~, Sd.Z] = max(E, [], 1);
                Sd.grps = intgroup(K, Sd.Z);
            else
                Sd.Q = nrmexp(E, 1);
            end
            
        end
          
                
        function Sp = make_output(prg, Sd, Sc)
            % Makes output struct from given states
            %
            %   Sp = prg.make_output(Sd, Sc);
            %       
            %       The output struct has the following fields:
            %       - U:    the estimated/sampled component means [d x K]
            %       - Cx:   the estimated/sampled component covariance
            %               (this is output only when est_Cx is true).
            %
            %       If it is doing sampling, Sp also has:
            %       - Z:    the label vector [1 x n]
            %
            %       If it is doing variational EM, Sp also has
            %       - Q:    the posterior distribution of labels [K x n]
            %
            
            Sp.U = Sd.U;
            
            if prg.est_Cx
                Sp.Cx = Sd.Cx;
            end
            
            if Sc.do_samp
                Sp.Z = Sd.Z;
            else
                Sp.Q = Sd.Q;
            end        
        end
        
        
        function objv = evaluate_objective(prg, Sd, Sc)
            % Evaluates the objective value w.r.t the given states
            %
            %   objv = prg.evaluate_objective(Sd, Sc);
            %
            
            gm = prg.gmodel;
            K = Sc.K;
            n = Sc.n;
            
            % log-prior of component params
            
            if ~isempty(prg.Cu) && ~isempty(Sc.mu)
                lpri_u = gm.logpri(Sc.mu, Sd.U);
                lpri_u = sum(lpri_u);
            else
                lpri_u = 0;
            end
            
            % log-prior of labeling 
            
            lpri_z = 0;         
            
            % log-likelihood
            
            if Sc.do_samp                
                llik = Sd.Lik(sub2ind([K, n], Sd.Z, 1:n));
                llik = sum(llik);                
            else
                llik = dot(Sd.Lik, Sd.Q, 1);
                llik = sum(llik);
            end                        
            
            % labeling entropy
            
            if Sc.do_samp
                ent_z = 0;
            else
                ent_z = sum(ddentropy(Sd.Q));
            end
            
            % combine
            
            objv = lpri_u + lpri_z + llik + ent_z;
            
        end
            
    end
    
end




