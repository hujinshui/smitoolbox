classdef gammad
    % The class to represent (one or multi-dimensional) Gamma distribution
    %
    %   Each object of this class can contain one or multiple Gamma
    %   distributions. A Gamma distribution is characterized by two
    %   parameters:
    %
    %   - alpha:    the shape parameter 
    %   - beta:     the scale parameter (inverse of scale)
    %
    %   The pdf is given by
    %
    %   f(x) = (alpha - 1) log(x) - x / beta 
    %        - (alpha log(beta) + gammaln(alpha))
    %  
    %   For an object with n distributions over d-dimensional space:    
    %   - shape is a matrix of size 1 x n or d x n.
    %   - scale can be a 1 x n row vector or a scale (shared)    
    %
    
    %   History
    %   -------
    %       - Created by Dahua Lin, on Aug 31, 2011
    %
    
    %% properties
    
    properties(GetAccess='public', SetAccess='private')
        
        dim;    % the dimension of the underlying space
        num;    % the number of distributions contained in the object
        
        alpha;  % the shape parameter(s)
        beta;   % the scale parameter(s)
        
        lpconst;   % a constant term in logpdf
                   % = -(alpha * log(beta) + gammaln(alpha))
                   % 1 x n row vector or 0, or empty (if not computed)
    end
        
    methods
        
        function v = mean(obj)
            % Evaluates the mean(s) of the distribution
            %
            %   v = mean(obj);
            %
            %       v will be an d x n matrix, with v(:,i) corresponding 
            %       to the i-th distribution in obj.
            %
        
            a = obj.alpha;
            b = obj.beta;
            d = obj.dim;
            
            if isscalar(b)
                if b == 1
                    v = a;
                else
                    v = a * b;
                end
            else
                v = bsxfun(@times, a, b);
            end
            
            if size(a, 1) ~= d
                v = v(ones(d, 1), :);
            end                
        end
                
        function v = var(obj)
            % Evaluates the variance(s) of the distribution
            %
            %   v = var(obj);
            %
            %       v will be an d x n matrix, with v(:,i) corresponding 
            %       to the i-th distribution in obj.
            %
            
            a = obj.alpha;
            b = obj.beta;
            d = obj.dim;
            
            if isscalar(b)
                if b == 1
                    v = a;
                else
                    v = a * (b^2);
                end
            else
                v = bsxfun(@times, a, (b.^2));
            end
            
            if size(a, 1) ~= d
                v = v(ones(d, 1), :);
            end
        end
            
        
        function v = mode(obj)
            % Gets the mode of the distribution
            %
            %   v = mode(obj);
            %
            %       v will be an d x n matrix, with v(:,i) corresponding 
            %       to the i-th distribution in obj.
            %
            
            a = obj.alpha;
            b = obj.beta;
            d = obj.dim;
            
            if isscalar(b)
                if b == 1
                    v = (a - 1);
                else
                    v = (a - 1) * b;
                end
            else
                v = bsxfun(@times, a - 1, b);
            end
            
            if size(a, 1) ~= d
                v = v(ones(d, 1), :);
            end 
            v = max(v, 0);
        end
        
        function v = entropy(obj)
            % Evaluates the entropy value(s) of the distribution
            %
            %   v = entropy(obj);
            %
            %       v will be a 1 x n matrix, with v(:,i) corresponding 
            %       to the i-th distribution in obj.
            %            
            
            a = obj.alpha;
            b = obj.beta;
            d = obj.dim;
            
            v = gamma_entropy(a);
            if ~isequal(b, 1)
                if isscalar(b) 
                    v = v + log(b);
                else
                    v = bsxfun(@plus, v, log(b));
                end
            end
            
            if size(v, 1) == 1
                if d > 1
                    v = v * d;
                end
            else
                v = sum(v, 1);
            end
        end
        
    end
    
    %% Construction
    
    methods 
        
        function obj = gammad(alpha, beta, d, op)
            % Constructs a Gamma distribution object
            %
            %   obj = gammad(alpha);
            %   obj = gammad(alpha, beta);
            %
            %       constructs a Gamma distribution object given 
            %       the parameters.
            %
            %       Inputs:
            %       - alpha:    the shape parameter of size d x n.
            %
            %       - beta:     the scale parameter, which can be
            %                   either a 1 x n row vector or a 
            %                   scale (if shared by all distributions).
            %
            %   obj = gammad(alpha, beta, d);
            %   
            %       To construct multi-dimensional gamma distributions,
            %       where each dimension has the same shape param, then
            %       one can use this syntax.
            %
            %       Here, shape can be input a 1 x n row vector, and
            %       use the 3rd argument to specify the dimension.
            %
            %   obj = gammad(alpha, beta, [], 'pre');
            %   obj = gammad(alpha, beta, d, 'pre');
            %
            %       Do pre-computation of the term 
            %       alpha * log(beta) + gammaln(alpha), which might
            %       speed-up the evaluation of logpdf or pdf later.
            %
                        
            % verify inputs
            
            if ~(isfloat(alpha) && ndims(alpha) == 2)
                error('gammad:invalidarg', ...
                    'alpha should be a numeric matrix.');
            end            
            [d_, n] = size(alpha);
            
            if ~( isfloat(beta) && (isscalar(beta) || ...
                    isequal(size(beta), [1 n])) )
                error('gammad:invalidarg', ...
                    'beta should be a scalar or a 1 x n row vector.');
            end
            
            if nargin < 3 || isempty(d)
                d = d_;
            else
                if ~(isnumeric(d) && isscalar(d) && d == fix(d) && d >= 1)
                    error('gammad:invalidarg', ...
                        'd should be a positive integer scalar.');
                end
                if d_ > 1 && d ~= d_
                    error('gammad:invalidarg', ...
                        'The size of shape is inconsistent with d.');
                end
            end
            
            if nargin < 4
                lpc = [];
            else
                if ~(ischar(op) && strcmpi(op, 'pre'))
                    error('gammad:invalidarg', ...
                        'The 4th argument can only be ''pre''.');
                end
                lpc = gammad.calc_lpconst(alpha, beta, d);
            end
                                    
            % create object
            
            obj.dim = d;
            obj.num = n;
            obj.alpha = alpha;
            obj.beta = beta;
            obj.lpconst = lpc;
        end            
    end
    
    
    %% Evaluation
    
    methods
       
        function L = logpdf(obj, X, si)
            % Evaluates log-pdf of given samples
            %
            %   L = obj.logpdf(X);
            %
            %       computes the logarithm of pdf at the samples given
            %       as columns of X.
            %
            %       Suppose there are m distributions in obj, and n 
            %       columns in X, then L is a matrix of size m x n, 
            %       where L(k, i) is the log-pdf of the i-th sample
            %       with respect to the k-th distribution.
            %
            %   L = obj.logpdf(X, si);
            %
            %       computes the logarithm of pdf with respect to the
            %       distributions selected by si.
            %
            
            d = obj.dim;
            if ~(isfloat(X) && ndims(X) == 2 && size(X,1) == d)
                error('gammad:invalidarg', ...
                    'X should be a numeric matrix with size(X,1) == dim.');
            end
            
            a = obj.alpha;
            b = obj.beta;                        
            lpc = obj.lpconst;
            
            if nargin >= 3 && ~isempty(si)
                a = a(:, si);
                if ~isscalar(b)
                    b = b(:, si);
                end
                if ~isempty(lpc)
                    lpc = lpc(1, si);
                end
            end
            
            if isempty(lpc)
                lpc = gammad.calc_lpconst(a, b, d);                
            end
                        
            if isequal(a, 1)
                T1 = 0;
            else
                if size(a, 1) == d
                    T1 = (a - 1)' * log(X);
                else
                    T1 = (a - 1)' * sum(log(X), 1);
                end
            end
            
            if d == 1
                sX = X;
            else
                sX = sum(X, 1);
            end
            if isequal(b, 1)
                T2 = sX;
            else
                T2 = (1./b)' * sX;
            end
            
            if size(T1, 1) == size(T2, 1)
                L = T1 - T2;
            else
                L = bsxfun(@minus, T1, T2);
            end
                
            if ~isequal(lpc, 0)
                L = bsxfun(@plus, L, lpc.');                                
            end                        
        end
        
        
        function L = pdf(obj, X, si)
            % Evaluates log-pdf of given samples
            %
            %   L = obj.pdf(X);
            %
            %       computes the pdf values at the samples given
            %       as columns of X.
            %
            %       Suppose there are m distributions in obj, and n
            %       columns in X, then L is a matrix of size m x n,
            %       where L(k, i) is the log-pdf of the i-th sample
            %       with respect to the k-th distribution.
            %
            %   L = obj.pdf(X, si);
            %
            %       computes the pdf values with respect to the
            %       distributions selected by si.
            %
            
            if nargin < 3
                L = exp(logpdf(obj, X));
            else
                L = exp(logpdf(obj, X, si));
            end            
        end
                
    end    
    
    
    methods(Static, Access='private')
        
        function v = calc_lpconst(a, b, d)
            % Calculates the log-pdf constant term
            %
            
            if isequal(b, 1)
                c1 = 0;
            else
                if d == 1
                    sa = a;
                else
                    if size(a, 1) == 1
                        sa = a * d;
                    else
                        sa = sum(a, 1);
                    end
                end
                c1 = sa .* log(b);
            end
            
            if isequal(a, 1)
                c2 = 0;
            else
                c2 = gammaln(a);
                if d > 1
                    if size(a, 1) == 1
                        c2 = c2 * d;
                    else
                        c2 = sum(c2, 1);
                    end
                end
            end
                
            v = -(c1 + c2);
        end    
    end
    
    
    %% Sampling
    
    methods
        
        function X = sample(obj, n, i) 
            % Samples from the Gamma distribution
            %            
            %   X = obj.sample();
            %   X = obj.sample(n);
            %       draws n samples from the Gamma distribution. 
            %       (It must be obj.num == 1 for this syntax).
            %
            %       When n is omitted, it is assumed to be 1.
            %
            %   X = obj.sample(n, i);
            %       draws n samples from the i-th Gamma distribution.            
            %
            %       When i is a vector, then n can be a vector of the 
            %       same size. In this case, it draws n(j) samples from
            %       the i(j)-th distribution.                     
            %
            
            if nargin < 2; n = 1; end
            
            a = obj.alpha;
            b = obj.beta;
            d = obj.dim;
            
            if nargin < 3 || isempty(i)
                if obj.num > 1
                    error('gammad:sample:invalidarg', ...
                        'i is needed when obj contains multiple models.');
                end
                
                X = gamma_sample(a, b, n, d);
                
            else
                if ~(isvector(n) && isnumeric(n))
                    error('gammad:pos_sample:invalidarg', ...
                        'n must be a numeric vector');
                end
                if ~(isvector(i) && isnumeric(i))
                    error('gammad:pos_sample:invalidarg', ...
                        'i must be a numeric vector');
                end
                if numel(n) ~= numel(i)
                    error('gammad:pos_sample:invalidarg', ...
                        'The sizes of n and i are inconsistent.');
                end
                
                a = a(:, i);
                if ~isscalar(b)
                    b = b(i);
                end                    
                K = size(a, 2);
                
                if K == 1
                    X = gamma_sample(a, b, n, d);
                else
                    
                    N = sum(n);
                    X = zeros(d, N, class(a));
                    ek = 0;
                    if isscalar(b)
                        b = b(ones(1, K));
                    end
                    
                    for k = 1 : K
                        sk = ek + 1;
                        ek = ek + n(k);
                        X(:, sk:ek) = gamma_sample(a(:,k), b(k), n(k), d);
                    end
                end
                
            end % end if there are n & i
            
        end % end sample function
    
    end
    
end




