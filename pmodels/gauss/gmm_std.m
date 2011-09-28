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
                
        gpri;       % The Gaussian prior distribution
        Cx;         % Cx itself or its prior distribution
        
        est_Cx;     % whether to estimate Cx  
    end
    
    
    properties
        dalpha = 1; % the concentration parameter of the Dirichlet prior of pi
                    % dird_alpha >= 1
                    % one can set dird_alpha to inf, in that case
                    % pi is fixed to [1/k, ..., 1/k].        
    end
    
    methods
        function prg = set.dalpha(prg, v)
            if ~(isfloat(v) && isreal(v) && v >= 1)
                error('gmm_std:invalidarg', ...
                    'The value of dalpha should be a real numberin [1, inf]');
            end
            prg.dalpha = v;
        end
    end
    
    
    %% Construction
    
    methods
        
        function prg = gmm_std(d, gpri, Cx)
            % Creates a standard Gaussian mixture model (GMM) program
            %                        
            %   prg = gmm_std(d);
            %   prg = gmm_std(d, gpri);
            %   prg = gmm_std(d, gpri, Cx);
            %
            %       Construsts an SMI program that implements the inference of
            %       GMM parameters.
            %
            %       Input arguments:
            %       - d:        the underlying space dimension
            %
            %       - gpri:     the Gaussian prior of the variable u.
            %                   (If gpri is omitted or empty, then no 
            %                    prior is used)
            %
            %       - Cx:       the covariance for generating x.
            %                   (If Cx is omitted or empty, then Cx is
            %                    to be estimated)
            %       
                        
            % verify input arguments
            
            if ~(isnumeric(d) && isscalar(d) && d == fix(d) && d > 0)
                error('gmm_std:invalidarg', 'd should be a positive integer.');
            end
            
            if nargin >= 2 && ~isempty(gpri)
                if ~(isa(gpri, 'gaussd') && gpri.num == 1 && gpri.dim == d)
                    error('gmm_std:invalidarg', ...
                        'The Gaussian prior gpri is invalid.');
                end
            else
                gpri = [];
            end
                        
            if nargin >= 3 && ~isempty(Cx)
                if ~(is_pdmat(Cx) && Cx.n == 1 && Cx.d == d)
                    error('gmm_std:invalidarg', ...
                        'The 2nd arg (Cx) is invalid.');
                end
            else
                Cx = [];
            end                                   
            
            % create object
            
            prg.dim = double(d);  
            prg.gpri = gpri;
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
            %                       (This field exists only when Cx is to
            %                        be estimated)
            %               - Pi:   the prior probabilities of components
            %                       (It is a 1 x K row vector, initialized
            %                       to be a uniform distribution)
            %               - Lik:  K x n pairwise likelihood table
            %
            %           If optype is 'sample', Sd also has
            %               - Z:    1 x n label vector
            %               - grps: 1 x K cell array of grouped indices
            %
            %           If optype is 'varinfer', Sd also has
            %               - Q:    K x n matrix: posterior label
            %                       distributions.
            %               
            %       - Sc:   fixed state struct, with fields:
            %               - X:        the observation matrix
            %               - n:        the number of samples in X
            %               - K:        the number of mixture components
            %               - Cx:       the covariance of each component
            %               - Jx:       the inverse covariance
            %               - do_samp:  whether it is doing sampling
            %
            %               Note Sc.Cx and Sc.Jx exist only when Cx is 
            %               fixed.
            %
            
            % verify inputs
            
            d = prg.dim;
            
            if ~(isstruct(obs) && isscalar(obs) && ...
                    isfield(obs, 'X') && isfield(obs, 'K'))
                error('gmm_std:invalidarg', 'obs is not a valid struct.');
            end
            
            X = obs.X;
            K = obs.K;
            
            if ~(isfloat(X) && ndims(X) == 2 && isreal(X) && size(X,1) == d)
                error('gmm_std:invalidarg', 'X should be a real matrix.');
            end
            n = size(X, 2);
            
            if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 1)
                error('gmm_std:invalidarg', 'K should be an integer.');
            end
            K = double(K);                                    
            
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
            Sc.do_samp = do_samp;              
            if ~prg.est_Cx
                Sc.Cx = prg.Cx;
                Sc.Jx = pdmat_inv(Sc.Cx);
            end
            
            Sd.U = [];
            if prg.est_Cx
                Sd.Cx = [];
            end
            Sd.Pi = constmat(1, K, 1 / K);
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
            
            g0 = prg.gpri;
            X = Sc.X;
            K = Sc.K;

            % assuming Cx is fixed, estimate U
            
            if Sc.do_samp
                [hp, Jp] = gaussgm_pos(g0, X, Sc.Jx, [], Sd.grps); 
                U = zeros(prg.dim, K);
                for k = 1 : K
                    U(:,k) = gsample(hp(:,k), pdmat_pick(Jp, k), 1, 'ip');
                end                
                Sd.U = U;
            else
                [~, ~, Sd.U] = gaussgm_pos(g0, X, Sc.Jx, [], Sd.Q); 
            end
            
            % infer labels (E-step)            
            
            H = pdmat_mvmul(Sc.Jx, Sd.U);
            cg = gaussd.from_ip(H, Sc.Jx);
            
            Lik = cg.logpdf(X);
            Sd.Lik = Lik;
            
            da = prg.dalpha;            
            if isfinite(da)            
                E = bsxfun(@plus, log(Sd.Pi).', Lik);
            end
            Q = nrmexp(E, 1);
            
            if Sc.do_samp
                Sd.Z = ddsample(Q, 1);
                Sd.grps = intgroup(K, Sd.Z);
            else
                Sd.Q = Q;
            end
            
            % estimate Pi
            
            if isfinite(da)     % if dalpha is inf, Pi is fixed
                
                if Sc.do_samp
                    tw = intcount(K, Sd.Z);                    
                    Sd.Pi = dirichlet_sample(K, tw(:) + da, 1).';                    
                else
                    tw = sum(Q, 2).';
                    if da > 1
                        tw = tw + (da - 1);
                    end
                    Sd.Pi = tw ./ sum(tw);
                end
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
            %       - Pi:   the estimated/sampled prior probability of
            %               components. [1 x K]
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
            
            Sp.Pi = Sd.Pi;
            
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
            
            g0 = prg.gpri;
            K = Sc.K;
            n = Sc.n;
            
            % log-prior of component params
            
            if ~isempty(g0)
                lpri_u = g0.logpdf(Sd.U);
                lpri_u = sum(lpri_u);
            else
                lpri_u = 0;
            end
            
            % log-prior of labeling and Pi
            
            logPi = log(Sd.Pi);
            
            if Sc.do_samp
                tw = intcount(K, Sd.Z).';
            else
                tw = sum(Sd.Q, 2);
            end
            lpri_z = logPi * tw;  
            
            da = prg.dalpha;
            if isfinite(da) && da > 1
                lpri_pi = sum(logPi) * (da - 1);
            else
                lpri_pi = 0;
            end
            
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
            
            objv = lpri_u + lpri_z + lpri_pi + llik + ent_z;
            
        end
            
    end
    
end




