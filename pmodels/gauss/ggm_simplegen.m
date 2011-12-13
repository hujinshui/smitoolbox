classdef ggm_simplegen
    % The class that implements a simple Gaussian generative model
    %
    %   The simple Gaussian generative model, with parameter u, is 
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
                    error('ggm_simplegen:invalidarg', ...
                        'The precision (Jx) should be a positive value.');
                end
            elseif is_pdmat(v)
                if ~(v.n == 1 && v.d == d)
                    error('ggm_simplegen:invalidarg', ...
                        'Jx should have Jx.n == 1 and Jx.d == d.');
                end
                model.Jx = v;
            else
                error('ggm_simplegen:invalidarg', ...
                    'Attempt to set Jx to an invalid value.');
            end
                        
            model.Gx.J = model.Jx;
            model.Gx_cb = d / 2 - gaussd_entropy(model.Jx, 'c');            
        end
        
        function model = set_A(model, v)
            % Set the parameter A to the model
            %
            %   model = set_A(model, A);
            %
            
            d = model.xdim;
            q = model.pdim;
            
            if ~(isfloat(v) && isreal(v) && isequal(size(v), [d q]))
                error('ggm_simplegen:invalidarg', ...
                    'A should be a real matrix of size d x q.');
            end
            
            model.A = v;
            model.use_A = true;
        end
        
    end
    
    
    
    %% methods
    
    methods
        
        %% constructor
        
        function model = ggm_simplegen(Jx, A)
            % Construct a Gaussian generative model as formulated above
            %
            %   model = ggm_simplegen(d);                       
            %   model = ggm_simplegen(d, q);
            %
            %       Creates an empty model with xdim == d and pdim == q.
            %       If q is omitted, then it assumes q == d.
            %
            %       After the construction, Jx and A remains empty.
            %       One has to set Jx (and optionally A) before using
            %       the model.
            %
            %   model = ggm_simplegen(Jx);
            %   model = ggm_simplegen(Jx, A);
            %       
            %       Creates a model with given precision matrix Jx, 
            %       and (optionally) the transform matrix A.
            %
            
            if isnumeric(Jx) && isscalar(Jx)
                d = Jx;                                
                if ~(isreal(d) && d == fix(d) && d > 0)
                    error('ggm_simplegen:invalidarg', 'd should be a positive integer.');
                end
                
                if nargin < 2
                    q = d;
                elseif isnumeric(A) && isscalar(A)
                    q = A;                                    
                    if ~(isreal(q) && q == fix(q) && q > 0)
                        error('ggm_simplegen:invalidarg', 'q should be a positive integer.');
                    end
                else
                    error('ggm_simplegen:invalidarg', 'The 2nd argument is invalid.');
                end

                model.xdim = double(d);
                model.pdim = double(q);
                                
            else
                
                if ~is_pdmat(Jx)
                    error('ggm_simplegen:invalidarg', 'Jx should be a pdmat struct.');
                end
                
                d = Jx.d;
                if nargin < 2
                    q = d;
                    uA = false;
                else
                    if ~(isfloat(A) && isreal(A) && ndims(A) == 2 && size(A,1) == d)
                        error('ggm_simplegen:invalidarg', ...
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
                error('ggm_simplegen:invalidarg', ...
                    'The observations should be a real matrix with d rows.');
            end
            n = size(X, 2);            
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
                error('ggm_simplegen:invalidarg', ...
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
        
        
        %% conjugate update
        
        function [dh, dJ] = capture(model, X, w)
            % Capture observations into conjugate updates 
            %
            %   [dh, dJ] = model.capture(X, w);
            %       computes the conjuate updates to the canonical params
            %       of the prior based on given (weighted) set of samples.
            %
            %       Inputs:
            %       - X:        the sample matrix, size: d x n
            %       - w:        the weights of the samples, size: m x n, 
            %                   or a scalar. If m > 1, then multiple 
            %                   updates are to be evaluated, each 
            %                   corresponding to a group in w.
            %
            %       Outputs:
            %       - dh:       the updates to the potential vector [q x m]
            %       - dJ:       a pdmat struct with dJ.n == m and 
            %                   dJ.d == q.
            %
            
            % verify inputs
            
            d = model.xdim;
            q = model.pdim;
            uA = model.use_A;    
            Jx_ = model.Jx;
            
            n = model.query_obs(X);
            
            if isempty(w)
                w = 1;
            else
                if ~(isfloat(w) && isreal(w) && ...
                        (isscalar(w) || (ndims(w) == 2 && size(w, 2) == n)))
                    error('ggm_simplegen:invalidarg', ...
                        'w should be a real scalar or a matrix with n columns.');
                end
                m = size(w, 1);
            end
            
            if Jx_.ty == 's'
                Jsca = 1;
                jv = Jx_.v;
            else
                Jsca = 0;
            end
            
            % compute dh
            
            if Jsca
                if isscalar(w)
                    dh = sum(X, 2) * (jv * w);
                else
                    dh = X * (jv * w)';
                end
            else
                JX = pdmat_mvmul(Jx_, X);
                if isscalar(w)
                    dh = sum(JX, 2) * w;
                else
                    dh = JX * w';
                end
            end
            
            if uA
                A_ = model.A;
                dh = A_' * dh;
            end
            
            
            % compute dJ
            
            if isscalar(w)
                tw = w * n;
            else
                tw = sum(w, 2).';
            end
            
            
            if ~uA
                if Jsca
                    dJ = pdmat('s', d, tw * jv);
                else
                    dJ = pdmat_scale(Jx_, tw);
                end
            else
                
                if Jsca
                    dJ0 = jv * (A_' * A_);
                else
                    dJ0 = pdmat_pwquad(Jx_, A_);
                    dJ0 = 0.5 * (dJ0 + dJ0');
                end
                
                if isscalar(tw)
                    dJ = pdmat('f', q, dJ0 * tw);
                else
                    dJ = zeros([size(dJ0), m]);
                    for k = 1 : m
                        dJ(:,:,k) = dJ0 * tw(k);
                    end
                    dJ = pdmat('f', q, dJ);
                end
            end
        end
        
        
    end
    
             
end


