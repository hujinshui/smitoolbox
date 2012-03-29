classdef gauss_lingen < genmodel_base
    % The class that implements a linear Gaussian generative model
    %
    %   The linear Gaussian generative model, with parameter u, is 
    %   formulated as
    %   
    %       x ~ N(A * u, Cx),   or x ~ N(u, Cx) if A is identity.
    %
    %   Here, let d be the dimension of x, and q be that of u, then
    %   A should be a d x q matrix.
    %
    %   Here, the observation can be captured by the conjugate update as
    %
    %       dh = sum_i w_i (A' * Jx * x_i);
    %       dJ = sum_i w_i (A' * Jx * A);
    %
    %   Here, Jx is the inverse of Cx, i.e. the precision matrix.
    %
    
    % Created by Dahua Lin, on Dec 13, 2011
    %
    
    
    %% properties
    
    % basic info
    
    properties(GetAccess='public', SetAccess='private')
        xdim;               % the dimension of observations (d)
        pdim;               % the dimension of parameters (q)
        use_A = false;      % whether transform matrix A is used.        
        Gx;                 % the measurement model
        Gx_cb;              % a constant value (cb) for Gx
    end
    
    % hyper-parameters
    
    properties(GetAccess='public', SetAccess='private')
        Jx;         % the precision matrix of the measurement model         
        A;          % the transform matrix [empty or d x q]
    end    
    
    properties(GetAccess='private', SetAccess='private')
        AtJA;       % the matrix equals A' * Jx * A
    end        
    
    
    methods
        
        function model = set_Jx(model, v)
            % Set the parameter Jx to the model
            %
            %   model = set_Jx(model, Jx);
            %
            
            d = model.xdim;
            if isscalar(v)
                if v > 0
                    model.Jx = pdmat_mat('s', d, v);
                else
                    error('gauss_lingen:invalidarg', ...
                        'The precision (Jx) should be a positive value.');
                end
            elseif is_pdmat(v)
                if ~(v.n == 1 && v.d == d)
                    error('gauss_lingen:invalidarg', ...
                        'Jx should have Jx.n == 1 and Jx.d == d.');
                end
                model.Jx = v;
            else
                error('gauss_lingen:invalidarg', ...
                    'Attempt to set Jx to an invalid value.');
            end
                        
            Jx_ = model.Jx;
            model.Gx.J = Jx_;
            model.Gx_cb = d / 2 - gaussd_entropy(Jx_, 'c');     
            
            if model.use_A
                model.AtJA = gauss_lingen.calc_AtJA(Jx_, model.A);
            end
        end
        
        function model = set_A(model, v)
            % Set the parameter A to the model
            %
            %   model = set_A(model, A);
            %
            
            d = model.xdim;
            q = model.pdim;
            
            if ~(isfloat(v) && isreal(v) && isequal(size(v), [d q]))
                error('gauss_lingen:invalidarg', ...
                    'A should be a real matrix of size d x q.');
            end
            
            model.A = v;
            model.use_A = true;
            if ~isempty(model.Jx)
                model.AtJA = gauss_lingen.calc_AtJA(Jx_, v);
            end
        end
        
    end
    
    
    
    %% methods
    
    methods
        
        %% constructor
        
        function model = gauss_lingen(Jx, A)
            % Construct a Gaussian generative model as formulated above
            %
            %   model = gauss_lingen(d);                       
            %   model = gauss_lingen(d, q);
            %
            %       Creates an empty model with xdim == d and pdim == q.
            %       If q is omitted, then it assumes q == d.
            %
            %       After the construction, Jx and A remains empty.
            %       One has to set Jx (and optionally A) before using
            %       the model.
            %
            %   model = gauss_lingen(Jx);
            %   model = gauss_lingen(Jx, A);
            %       
            %       Creates a model with given precision matrix Jx, 
            %       and (optionally) the transform matrix A.
            %
            
            if isnumeric(Jx) && isscalar(Jx)
                d = Jx;                                
                if ~(isreal(d) && d == fix(d) && d > 0)
                    error('gauss_lingen:invalidarg', 'd should be a positive integer.');
                end
                
                if nargin < 2
                    q = d;
                elseif isnumeric(A) && isscalar(A)
                    q = A;                                    
                    if ~(isreal(q) && q == fix(q) && q > 0)
                        error('gauss_lingen:invalidarg', 'q should be a positive integer.');
                    end
                else
                    error('gauss_lingen:invalidarg', 'The 2nd argument is invalid.');
                end

                model.xdim = double(d);
                model.pdim = double(q);
                                
            else
                
                if ~is_pdmat(Jx)
                    error('gauss_lingen:invalidarg', 'Jx should be a pdmat struct.');
                end
                
                d = Jx.d;
                if nargin < 2
                    q = d;
                    uA = false;
                else
                    if ~(isfloat(A) && isreal(A) && ndims(A) == 2 && size(A,1) == d)
                        error('gauss_lingen:invalidarg', ...
                            'A should be a real matrix with d rows.');
                    end
                    q = size(A, 2);
                    uA = true;
                end
                
                model.xdim = d;
                model.pdim = q;
                model.use_A = uA;
                model.Jx = Jx;
                
                if uA
                    model.A = A;
                    model.AtJA = gauss_lingen.calc_AtJA(Jx, A);
                end
                                
                model.Gx = gaussd('c', 0, Jx);
                model.Gx_cb = d / 2 - gaussd_entropy(Jx, 'c');                
            end            
        end                        
        
        
        %% observation query
        
        function n = query_obs(model, X)
            % Get the number of observation samples
            %
            %   n = model.query_obs(X);
            %       verifies the validity of X as an observation set,
            %       and returns the number of samples in X.
            %
            
            d = model.xdim;
            if ~(isfloat(X) && isreal(X) && ndims(X) == 2 && size(X,1) == d)
                error('gauss_lingen:invalidarg', ...
                    'The observations should be a real matrix with d rows.');
            end
            n = size(X, 2);
        end
        
        
        function n = query_params(model, U)
            % Get the number of parameters
            %
            %   n = model.query_params(U);
            %       verifies the validity of X as an observation set,
            %       and returns the number of samples in X.
            %
            
            q = model.pdim;
            if ~(isfloat(U) && isreal(U) && ndims(U) == 2 && size(U,1) == q)
                error('gauss_lingen:invalidarg', ...
                    'The observations should be a real matrix with q rows.');
            end
            n = size(U, 2);
        end
        
        function Uc = combine_params(model, varargin) %#ok<MANU>
            % Combine multiple parameters
            %
            %   Uc = model.combine_params(U1, U2, ...);
            %
            
            Uc = [varargin{:}];
        end
        
        
        %% log-likelihood evaluation 
        
        function LL = loglik(model, U, X)
            % Evaluate the log-likelihood values at given samples
            %
            %   LL = loglik(model, U, X);
            %       evaluates the log-likelihood at the samples given 
            %       in X, with respect to the parameters given in U.
            %
            
            % verify inputs
            
            q = model.pdim;
            if ~(isfloat(U) && isreal(U) && ndims(U) == 2 && size(U, 1) == q)
                error('gauss_lingen:invalidarg', ...
                    'The params U should be a real matrix with q rows.');
            end
            
            % evaluate
            
            if model.use_A
                U = model.A * U;
            end
                                  
            g = model.Gx;
            g.n = size(U, 2);
            g.h = pdmat_mvmul(g.J, U);
            ca = sum(g.h .* U, 1);
            cb = model.Gx_cb;
            
            LL = gaussd_logpdf(g, X, {ca, cb});            
        end        
        
        
        %% maximum likelihood estimation
        
        function U = mle(model, X, W, I)
            % Performs maximum likelihood estimation of the parameters
            %
            %   U = model.mle(X, W);
            %   U = model.mle(X, W, I);
            %
            %       performs maximum likelihood estimation based on
            %       given (weighted) set of data
            %
            
            % verify inputs
            
            if nargin < 3
                W = [];
            end
                                    
            if nargin >= 4
                X = X(:, I);
                if ~isempty(W)
                    W = W(I, :);
                end
            end    
            
            n = model.query_obs(X);
                        
            % compute
            
            if isempty(W)
                U = sum(X, 2) * (1 / n);
            else
                sw = sum(W, 1);
                U = X * bsxfun(@times, W, 1 ./ sw);
            end
            
            if model.use_A
                H_ = model.AtJA;
                A_ = model.A;
                
                U = H_ \ (A_' * pdmat_mvmul(model.Jx, U));
            end            
        end                
                
        
        %% conjugate update
        
        function S = capture(model, X, W, I)
            % Capture observations into conjugate updates 
            %
            %   S = model.capture(X, W);
            %   S = model.capture(X, W, I);
            %       computes the conjuate updates to the canonical params
            %       of the prior based on given (weighted) set of samples.
            %
            %       Inputs:
            %       - X:        the sample matrix, size: d x n
            %       - W:        the weight matrix, empty or n x K
            %       - I:        the selected indices
            %
            %       Outputs:
            %       - I:
            %
            
            % verify inputs
            
            d = model.xdim;
            q = model.pdim;
            uA = model.use_A;    
            Jx_ = model.Jx;
            
            if nargin < 3
                W = [];
            end
            
            if nargin >= 4
                X = X(:, I);
                if ~isempty(W)
                    W = W(I, :);
                end
            end
            
            n = model.query_obs(X);            
            
            if Jx_.ty == 's'
                Jsca = 1;
                jv = Jx_.v;
            else
                Jsca = 0;
            end
            
            % compute dh
            
            if Jsca
                if isempty(W)
                    dh = sum(X, 2) * jv;
                else
                    dh = X * (jv * W);
                end
            else
                JX = pdmat_mvmul(Jx_, X);
                if isempty(W)
                    dh = sum(JX, 2);
                else
                    dh = JX * W;
                end
            end
            
            if uA
                A_ = model.A;
                dh = A_' * dh;
            end
            
            
            % compute dJ
            
            if isempty(W)
                K = 1;
                tw = n;
            else
                K = size(W, 2);
                tw = sum(W, 1);
            end
            
            
            if ~uA
                if Jsca
                    dJ = pdmat('s', d, tw * jv);
                else
                    dJ = pdmat_scale(Jx_, tw);
                end
            else
                
                dJ0 = model.AtJA;
                
                if isscalar(tw)
                    dJ = pdmat('f', q, dJ0 * tw);
                else
                    dJ = zeros([size(dJ0), K]);
                    for k = 1 : K
                        dJ(:,:,k) = dJ0 * tw(k);
                    end
                    dJ = pdmat('f', q, dJ);
                end
            end
            
            % make S
            
            S.tag = 'gaussd';
            S.ty = 'c';
            S.n = K;
            S.d = q;
            S.h = dh;
            S.J = dJ;            
        end
        
        
    end
    
    
    %% Auxiliary implementation
    
    methods(Static, Access='private')
        
        function H = calc_AtJA(J, A)            
            if J.ty == 's'
                jv = J.v;
                H = jv * (A' * A);
            else
                H = pdmat_pwquad(J, A);
                H = 0.5 * (H + H');
            end            
        end
        
    end
    
    
             
end


