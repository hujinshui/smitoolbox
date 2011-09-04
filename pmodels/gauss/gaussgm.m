classdef gaussgm
    % The class to implement a Linear Gaussian generative model 
    %
    % The model is formulated as follows
    %
    %   generate params:
    %   -----------------
    %   u  ~ Gaussian(B * mu, Cu) or Gaussian(mu, Cu) when B = 1
    %   Cx ~ Inverse Wishart(Phi, invws_deg); 
    %
    %   generate samples:
    %   ------------------
    %   x  ~ Gauss(A * u,  Cx) or Gauss(x,  Cx) when A = 1
    %
    %   The model have two different settings:
    %
    %   (1) Covariance Cx is fixed at design stage.
    %
    %       Output variables:   
    %       - x:    [xdim x 1 double vector]
    %       
    %       Model Parameters:
    %       - u:    [udim x 1 double vector]
    %
    %       Hyper parametes:
    %       - mu:   [dim0 x 1 double vector]
    %
    %       Design parameters:    
    %       - Cx:   [pdmat struct with d == xdim]
    %       - Cu:   [pdmat struct with d == udim]
    %       - A:    [xdim x udim double matrix or a scalar]
    %       - B:    [udim x dim0 double matrix or a scalar]
    %       
    %   (2) Covariance Cx servs as a parameter, which can be a single
    %       covariance matrix shared by all models, or multiple matrices,
    %       each used in a distinct model.
    %
    %   Note that if Cx is fixed, some computation can be done in advance
    %   in constructor, which could speedup the run-time evaluation.
    %    
    
    % Created by Dahua Lin, on Aug 26, 2011
    %
    
    %% properties
    
    properties(GetAccess='public', SetAccess='private')
        
        dim0;       % the dimension of prior mean (mu)
        udim;       % the dimension of the mean parameter (u)
        xdim;       % the dimension of the output variables (x)
                
        Cx;         % the covariance in generating x
                    % (used only in Cx-fixed setting)
        
        Cu;         % the prior covariance in generating u
        A;          % the u --> x transform matrix or a scalar
        B;          % the mu --> u transform matrix or a scalar
        
        Gx;         % a Gaussian distribution: N(0, Cx)
                    % (used only in Cx-fixed setting)
                    
        Gu;         % a Gaussian distribution: N(0, Cu)
        
        fixed_Cx;   % whether Cx is a fixed design parameter
    end
            
    
    %% Construction
    
    methods
        
        function model = gaussgm(Cx, Cu, A, B)
            % Constructs a linear Gaussian generative model
            %
            %   model = gaussgm(d);
            %   model = gaussgm(d, Cu);
            %
            %   model = gaussgm(d, du, A);
            %   model = gaussgm(d, Cu, A);
            %   model = gaussgm(d, du, A, B);
            %   model = gaussgm(d, Cu, A, B);
            %       
            %       Constructs a linear Gaussian generative model, of
            %       which both u and Cx are considered as parameters. 
            %
            %       Here, d is the dimension of output space (i.e. 
            %       the variable x's underlying space).
            %
            %       If no prior is to be used, then one can omit Cu
            %       or just input du instead, which is the dimension
            %       the vector u.
            %
            %       Note A or B can be omitted, in which case, it is
            %       assumed to be identity (i.e. scalar 1).
            %
            %   model = gaussgm(Cx, Cu);
            %   model = gaussgm(Cx, Cu, A);
            %   model = gaussgm(Cx, Cu, A, B);
            %
            %       Constructs a linear Gaussian generative model with
            %       Cx fixed as a design parameter, which should be
            %       input as a pdmat struct.
            %            
            
            % verify arguments
            
            if isnumeric(Cx) && isscalar(Cx)
                xd = Cx;
                Cx = [];
            elseif is_pdmat(Cx)
                xd = Cx.d;
                if ~isa(Cx.v, 'double'); Cx.v = double(Cx.v); end
            else
                error('gaussgm:invalidarg', 'The first arg is invalid.');
            end            
            
            if nargin < 2
                ud = xd;
                Cu = [];
            elseif isnumeric(Cu) && isscalar(Cu)
                ud = Cu;
                Cu = [];
            elseif is_pdmat(Cu)
                ud = Cu.d;
                if ~isa(Cu.v, 'double'); Cu.v = double(Cu.v); end
            else
                error('gaussgm:invalidarg', 'The second arg is invalid.');                
            end            
                                                
            if nargin < 3
                A = 1;                
            else
                if ~( isfloat(A) && ...
                        (isscalar(A) || isequal(size(A), [xd, ud])) )
                    
                    error('gaussgm:invalidarg', ...
                        ['A should be either ', ...
                        'a numeric matrix of size xdim x udim or a scalar.']);
                end                
                if ~isa(A, 'double'); A = double(A); end
            end
            
            if isscalar(A) && (xd ~= ud)
                error('gaussgm:invalidarg', ...
                    'A cannot be a scalar when xdim ~= udim.');
            end            
            
            if nargin < 4
                B = 1;
                d0 = ud;
            else
                if ~( isfloat(B) && ...
                        (isscalar(B) || (ndims(B)==2 && size(B,1)==xd)) )
                    
                    error('gaussgm:invalidarg', ...
                        ['B should be either ', ...
                        'a numeric matrix of size ydim x dim or a scalar.']);
                end                
                d0 = size(B, 2);
                if ~isa(B, 'double'); B = double(B); end
            end
            
            % make object
            
            model.dim0 = d0;
            model.udim = ud;
            model.xdim = xd;
            
            model.Cx = Cx;
            model.Cu = Cu;
            model.A = A;
            model.B = B;
            
            if ~isempty(Cx)
                model.Gx = gaussd.from_mp(0, Cx, 'ip');
            end
            if ~isempty(Cu)
                model.Gu = gaussd.from_mp(0, Cu, 'ip');
            end
            
            model.fixed_Cx = ~isempty(Cx);
        end
        
    end
    
    %% Evaluation methods
    
    methods
            
        function lpri = logpri(model, mu, U, hmap)
            % Evaluates the log-prior of given parameters
            %
            %   lpri = model.logpri(mu, u);
            %
            %       Evaluates the log-prior of the parameters in X 
            %       with respect to the hyper-parameter mu.      
            %
            %   lpri = model.logpri(Mu, U, hmap);
            %
            %       Evaluates the log-prior of the parameters in X
            %       with respect to multiple hyper-parameter in Mu.
            %
            %       In particular, the log-prior value of the k-th
            %       parameter in X should be evaluated based on 
            %       hmap(k)-th prior mean in Mu.
            %
            %   For both cases, lpri is a vector of size 1 x m, where
            %   m is the number of columns in U, such that lpri(k) is
            %   the log-prior value for U(:,k).
            %
            
            % verify inputs
            
            if isempty(model.Gu)
                error('gaussgm:invalidarg', ...
                    'prior is not supported with this model object.');
            end
                        
            mu = chk_mu(model, mu);         
            U = chk_u(model, U);
            
            if nargin < 4 || isempty(hmap)
                if size(mu, 2) ~= 1
                    error('gaussgm:invalidarg', ...
                        'mu should have one column when hmap is not given.');
                end
                hmap = [];
                zmu = all(mu == 0);
            else
                nx = size(X, 2);
                if ~(isnumeric(hmap) && isequal(size(hmap), [1 nx]))
                    error('gaussgm:invalidarg', ...
                        'hmap should be a numeric row vector of length size(X,2).');
                end
                zmu = 0;
            end            
            
            % compute
            
            if zmu
                V = U;      % V = U - mu
            else 
                B_ = model.B;
                if isequal(B_, 1)
                    bmu = mu;
                else
                    bmu = B_ * mu;
                end
                                
                if ~isempty(hmap)
                    bmu = bmu(:, hmap);
                end
                
                if size(U, 2) == size(bmu, 2)
                    V = U - bmu;
                else
                    V = bsxfun(@minus, U, bmu);
                end
            end
            
            lpri = model.Gu.logpdf(V);   
        end
        
                
        function Llik = loglik(model, U, X)
            % Evaluates the generative log-likelihood 
            %
            %   Llik = model.loglik(U, X);
            %
            %       Evaluates the log-likelihood of Y with respect
            %       to the model with parameters given by X, in a
            %       pairwise manner.
            %
            %       Suppose m = size(U,2) and n = size(X,2), then
            %       Llik is a matrix of size m x n, such that 
            %       Llik(k, i) is the likelihood value of the sample
            %       X(:,i) with respect to the parameter U(:,k).
            %
            
            % verify inputs
            
            U = chk_u(model, U);
            X = chk_x(model, X);
            
            % compute
            
            A_ = model.A;
            if isequal(A_, 1)
                AU = U;
            else
                AU = A_ * U;
            end
            
            m = size(U, 2);
            
            if m == 1
                gx0 = model.Gx;
                Llik = gx0.logpdf(bsxfun(@minus, X, AU));
            else
                gx = gaussd.from_mp(AU, model.Cx, 'ip');
                Llik = gx.logpdf(X);
            end            
        end
        
    end
    
    
    %% Sampling & Inference methods
    
    methods
    
        function X = generate_x(model, U, n, g)
            % Generates output variable x given parameters and labels
            %
            %   X = model.generate_x(u, n);
            %
            %       generates n sample from the model with parameter u.
            %
            %   X = model.generate_x(U, n, g);
            %
            %       generates samples from multi-models with parameters 
            %       given as columns of U, based on the grouping given 
            %       in g, a cell array of index vectors.
            %
            %       In particular, the samples with indices in g{k} should
            %       be generated by the k-th model (U(:,k)).
            %
            %       Here, n is the total number of samples.
            %
            
            % verify inputs
            
            U = chk_u(model, U);
            
            if ~(isnumeric(n) && isscalar(n) && n == fix(n) && n >= 0)
                error('gaussgm:invalidarg', ...
                    'n should be a non-negative integer scalar.');
            end                                
            
            if nargin < 4 || isempty(g)
                if size(U, 2) ~= 1
                    error('gaussgm:invalidarg', ...
                        'U must have only one column when no grouping.');
                end
                is_grp = 0;
            else
                if ~(iscell(g) && numel(g) == size(U, 2))
                    error('gaussgm:invalidarg', ...
                        ['g should be a cell array of index vectors', ...
                        'with K cells (K = size(X,2)).']);
                end
                is_grp = 1;
            end
               
            % compute
            
            A_ = model.A;
            if isequal(A_, 1)
                AU = U;
            else
                AU = A_ * U;
            end
                        
            X = model.Gx.sample(n);
            
            if ~is_grp                
                X = bsxfun(@plus, Y, AU);                
            else
                for k = 1 : numel(g)
                    gk = g{k};
                    X(:, gk) = bsxfun(@plus, X(:, gk), AU(:,k));
                end
            end
            
        end
        
        
        function U = mapest_u(model, X, CovX, Mu, Z, hmap)
            % Maximum-a-posteriori estimation of mean parameter u
            %
            %   u = model.mapest_u(X);
            %   u = model.mapest_u(X, Cx);            
            %       
            %       Performs maximum-likelihood estimation (MLE) of the
            %       mean parameter u given a set of samples X, and the
            %       covariance in generating X from u.
            %
            %       When Cx has already been a fixed design parameter,
            %       it need not be given. However, if it is given, the
            %       input Cx will be used instead of the one fixed in
            %       the model.
            %
            %   u = model.mapest_u(X, [], mu);      
            %   u = model.mapest_u(X, Cx, mu);
            %
            %       Performs maximum-a-posteriori (MAP) estimation of
            %       the mean parameter u given a set of samples X, and
            %       the prior mean mu.                   
            %
            %   U = model.mapest_u(X, Cx, mu, W);
            %   U = model.mapest_u(X, Cx, mu, W);
            %
            %       Performs ML or MAP estimation based on the given 
            %       weighted samples. 
            %
            %       Specifically, if W is a matrix of size K x n (where
            %       n = size(X,2)), then K parameters will be estimated,
            %       and returned as columns in U, a udim x K matrix.
            %
            %       In particular, if W is a row vector, then only one
            %       parameter is estimated.
            %
            %       Here, Cx can be set to empty if it is not needed,
            %       and mu can also be set to empty, when MLE is to be
            %       performed.            
            %
            %   U = model.mapest_u(X, Cx, mu, g);
            %   U = model.mapest_u(X, Cx, mu, g);
            %
            %       Performs ML or MAP estimation based on the given 
            %       grouped samples.
            %
            %       Specifically, g is a cell array of index vectors.
            %       Let K = numel(g), then K parameters are to be 
            %       estimated, and returned as columns in U, where
            %       U(:, k) is estimated based on the samples whose
            %       indices are in g{k}.      
            %
            %       Here, Cx can be set to empty if it is not needed,
            %       and mu can also be set to empty, when MLE is to be
            %       performed.
            %
            %   U = model.mapest_u(X, Cx, Mu, W, hmap);
            %   U = model.mapest_u(X, Cx, Mu, g, hmap);
            %
            %       Performs MAP estimation with prior-map (hmap).
            %       
            %       To estimate K parameters, hmap should be a 1 x K
            %       index vector. Then U(:,k) is estimated based on
            %       the prior mean given by Mu(:, hmap(k)).
            %
            
            if nargin < 3; CovX = []; end
            if nargin < 4; Mu = []; end
            if nargin < 5; Z = []; end
            if nargin < 6; hmap = []; end
            
            [H, J] = compute_pos_ip(model, X, CovX, Mu, Z, hmap);    
            U = pdmat_lsolve(J, H);
        end
        
        
        function X = sample_u(model, X, CovX, Mu, Z, hmap)
            % Sample from the posterior distribution of u
            %            
            %   u = model.sample_u(X, [], mu);      
            %   u = model.sample_u(X, Cx, mu);
            %
            %       Samples from the posterior distribution of u, with
            %       observed samples given by X, and prior mean given
            %       by mu.
            %
            %       When Cx has already been a fixed design parameter,
            %       it need not be given. However, if it is given, the
            %       input Cx will be used instead of the one fixed in
            %       the model.    
            %
            %   U = model.sample_u(X, Cx, mu, W);
            %   U = model.sample_u(X, Cx, mu, W);
            %
            %       Samples from the posterior distribution of u, based 
            %       on the given weighted samples. 
            %
            %       Specifically, if W is a matrix of size K x n (where
            %       n = size(X,2)), then K parameters will be sampled,
            %       and returned as columns in U, a udim x K matrix.
            %
            %       In particular, if W is a row vector, then only one
            %       parameter is sampled.
            %
            %       Here, Cx can be set to empty if it is not needed,           
            %
            %   U = model.mapest_u(X, Cx, mu, g);
            %   U = model.mapest_u(X, Cx, mu, g);
            %
            %       Samples from the posterior distribution of u, based 
            %       on the given grouped samples.
            %
            %       Specifically, g is a cell array of index vectors.
            %       Let K = numel(g), then K parameters are to be 
            %       sampled, and returned as columns in U, where
            %       U(:, k) is sampled based on the samples whose
            %       indices are in g{k}.      
            %
            %       Here, Cx can be set to empty if it is not needed,
            %
            %   U = model.mapest_u(X, Cx, Mu, W, hmap);
            %   U = model.mapest_u(X, Cx, Mu, g, hmap);
            %
            %       Samples from the posterior distribution of u, with
            %       a prior-map (hmap).
            %       
            %       To estimate K parameters, hmap should be a 1 x K
            %       index vector. Then U(:,k) is sampled based on
            %       the prior mean given by Mu(:, hmap(k)).
            %
            
            if nargin < 3; CovX = []; end
            if nargin < 4; Mu = []; end
            if nargin < 5; Z = []; end
            if nargin < 6; hmap = []; end
            
            if isscalar(hmap)
                ns = hmap;
                hmap = [];
            else
                ns = 1;
            end
            
            [H, J] = compute_pos_ip(model, X, CovX, Mu, Z, hmap);    
            U = pdmat_lsolve(J, H);  
            C = pdmat_inv(J);

            K = size(H, 2);
            if K == 1               
                X = gsample(U, C, ns);
            else
                % due to the computation of posterior
                % C is not shared across different posterior distributios                
                X = zeros(model.xdim, K);
                for k = 1 : K
                    Ck = pdmat_pick(C, k);
                    X(:, k) = gsample(U(:,k), Ck, 1);
                end
            end
                        
        end
        
    end
    
    
    %% Core inference implementation
    
    methods(Access='private')
        
        function [H, J] = compute_pos_ip(model, X, CovX, Mu, Z, hmap)
            % compute information parameter of posterior of u
            
            %% verify input arguments
            
            xd = model.xdim;
            ud = model.udim;
            d0 = model.dim0;
            
            X = chk_x(model, X);
            n = size(X, 2);
            
            if isempty(CovX)
                if ~model.fixed_Cx
                    error('gaussgm:invalidarg', ...
                    'Cx should be explicitly given when it is not fixed.');
                end
            else
                if ~(is_pdmat(CovX) && CovX.d == xd)                    
                    error('gaussgm:invalidarg', 'Cx is invalid.');
                end
            end
            
            if ~isempty(Mu)
                if ~(isfloat(Mu) && ndims(Mu) && size(Mu,1) == d0)
                    error('gaussgm:invalidarg', ...
                        'Mu should be a numeric matrix with dim0 rows.');
                end                
            end
            
            if ~isempty(Z)
                if isfloat(Z)
                    W = Z;
                    K = size(W, 1);
                    is_grp = false;
                elseif iscell(Z)
                    g = Z;
                    K = numel(g);
                    is_grp = true;
                else
                    error('gaussgm:invalidarg', ...
                        'Invalid argument for weighting/grouping.');
                end
            else
                W = [];
                K = 1;
                is_grp = false;
            end
            
            if ~isempty(CovX)
                if ~(CovX.n == 1 || CovX.n == K)
                    error('gaussgm:invalidarg', ...
                        'Cx.n should equal either 1 or K.');
                end
            end
            
            if ~isempty(hmap)
                if ~(isnumeric(hmap) && isequal(size(hmap), [1, K]))
                    error('gaussgm:invalidarg', ...
                        'hmap should be a numeric vector of size 1 x K.');
                end
            else                
                if ~isempty(Mu) && size(Mu,2) ~= 1
                    error('gaussgm:invalidarg', ...
                        'Without hmap, Mu should be a column vector.');
                end
            end
            
            use_pri = ~isempty(Mu);            
            if use_pri && isempty(model.Cu)
                error('gaussgm:invalidarg', ...
                    'Prior is not supported with this model object.');
            end
            
            %% do computation
            
            % formulas to implement:
            % 
            % h_k <- inv(Cx_k) * (sum_i w_{ki} x_i)
            % J_k <- A' * (sum_i w_i) * inv(Cx) * A
            %
            
            % get Jx = inv(Cx)
            
            if isempty(CovX)
                Jx = model.Gx.J;
            else
                Jx = pdmat_inv(CovX);
            end
            
            % sums up stuff
            
            if ~is_grp
                if isempty(W)
                    sw = n;
                    sx = sum(X, 2);
                else
                    sw = sum(W, 2)';
                    sx = X * W';
                end                
            else
                sw = zeros(1, K);
                sx = zeros(xd, K);
                
                for k = 1 : K
                    gk = g{k};
                    sw(k) = numel(gk);
                    sx(:,k) = sum(X(:,gk), 2);
                end
            end                           
            
            % compute H
            
            H = pdmat_mvmul(Jx, sx);
                                  
            % compute J
            
            A_ = model.A;
            if isscalar(A_)                    
                J = pdmat_scale(Jx, sw * (A_^2));                
            else
                if Jx.n == 1
                    Jmat = pdmat_quad(Jx, A_);
                    J = pdmat_scale(pdmat(Jmat), sw);
                    
                else % Jx.n == K > 1
                    if ud == 1
                        Jvs = zeros(1, K);
                        for k = 1 : K
                            Jx_k = pdmat_pick(Jx, k);
                            Jvs(k) = pdmat_quad(Jx_k, A_);
                        end
                        J = pdmat_scale('s', 1, Jvs);
                    else
                        Jms = zeros(ud, ud, K);
                        for k = 1 : K
                            Jx_k = pdmat_pick(Jx, k);
                            Jms(:,:,k) = pdmat_quad(Jx_k, A_);
                        end
                        J = pdmat_scale('f', ud, Jms);
                    end
                                        
                end
            end
            
            % incorporate prior
            
            if use_pri
                
                Ju = model.Gu.J;
                
                if ~(size(Mu, 2) == 1 && all(Mu == 0)) % need add to H
                    
                    B_ = model.B;
                    if isequal(B_, 1)
                        Bmu = Mu;
                    else
                        Bmu = B_ * Mu;
                    end
                    
                    if ~isempty(hmap)
                        Bmu = Bmu(:, hmap);
                    end
                    
                    dH = pdmat_mvmul(Ju, Bmu);
                    
                    if size(H, 2) == size(dH, 2)
                        H = H + dH;
                    else
                        H = bsxfun(@plus, H, dH);
                    end                    
                end
                
                % add to J          
                J = pdmat_plus(J, Ju);
            end
                        
        end        
    end
    
    
    
    %% Auxiliary methods
    
    methods
        function X = chk_x(model, X)
            xd = model.xdim;
            if ~(isfloat(X) && isreal(X) && ndims(X) == 2 && size(X,1) == xd)
                error('gaussgm:invalidarg', ...
                    'X should be a real matrix with xdim rows.');
            end
            if ~isa(X, 'double')
                X = double(X);
            end
        end        
        
        function U = chk_u(model, U)            
            ud = model.udim;
            if ~(isfloat(U) && isreal(U) && ndims(U) == 2 && size(U,1) == ud)
                error('gaussgm:invalidarg', ...
                    'U should be a real matrix with udim rows.');
            end
            if ~isa(U, 'double')
                U = double(U);
            end
        end
        
        function Mu = chk_mu(model, Mu)
            d0 = model.dim0;
            if ~(isfloat(Mu) && isreal(Mu) && ndims(Mu) == 2 && size(Mu,1) == d0)
                error('gaussgm:invalidarg', ...
                    'Mu should be a real matrix with dim0 rows.');
            end
            if ~isa(Mu, 'double')
                Mu = double(Mu);
            end
        end
    end
    
                
end

