classdef dirichletd
    % The class to represent a Dirichlet distribution
    %
    %   Each object of this class can contain one or multiple
    %   Dirichlet distributions. A Dirichlet distribution over 
    %   a (K-1)-dimensional simplex is characterized by a 
    %   K-dimensional vector alpha.
    %
    %   The pdf is given by
    %
    %   f(x) = G(sum alpha) / sum G(alpha) * prod_i x(i)^(alpha(i)-1)
    %
    %   Here, each x is a K-dimensional vector with sum(x) == 1, and
    %   G denotes the gamma function.
    %
    %   For an object with n distribution over a (K-1)-dimensional
    %   probability simplex, alpha can be either a K x n matrix, or
    %   a 1 x n row vector (for symmetric Dirichlet distribution).
    %   
    
    % Created by Dahua Lin, on Sep 27, 2010
    %
    
    %% Properties
    
    properties        
        K;          % the dimension of the simplex is (dim - 1)
        num;        % the number of distributions contained in the object
        
        alpha;      % the parameter(s): K x n or 1 x n
                
        logB;       % a constant term in logpdf: 1 x n
                    % = sum(gammaln(alpha)) - gammaln(sum(alpha)).
                    % or is empty (if not pre-computed)
    end
    
    %% Statistics
    
    methods
        
        function v = mean(obj)
            % Evaluates the mean(s) of the distributions
            %
            %   v = mean(obj);
            %
            %       v will be an K x n matrix, with v(:,i) corresponding
            %       to the i-th distribution in obj.
            %
            
            k = obj.K;
            n = obj.num;                        
            
            if k > 1
                a = obj.alpha;                                
                if size(a, 1) == 1  % symmetric
                    v = constmat(k, n, 1 / k);
                else                % non-symmetric
                    s = sum(a, 1);
                    if n == 1
                        v = a ./ s;
                    else
                        v = bsxfun(@times, a, 1 ./ s);
                    end
                end
            else
                v = ones(1, n);
            end
        end
        
        
        function v = var(obj)
            % Evaluates the variance(s) of the distributions
            %
            %   v = var(obj);
            %
            %       v will be an K x n matrix, with v(:,i) corresponding
            %       to the i-th distribution in obj.
            %
            
            k = obj.K;
            n = obj.num;
            
            if k > 1
                a = obj.alpha;
                if size(a, 1) == 1  % symmetric
                    p = 1 / k;
                    v = p * (1 - p) ./ (k * a + 1);
                    v = v(ones(1, k), :);
                else                % non-symmetric
                    s = sum(a, 1);
                    if n == 1
                        p = a ./ s;
                        v = p .* (1 - p) * (1 / (s + 1));
                    else
                        p = bsxfun(@times, a, 1 ./ s);
                        v = bsxfun(@times, p .* (1 - p), 1 ./ (s + 1));
                    end                    
                end                
            else
                v = zeros(1, n);
            end
        end
        
        
        function C = cov(obj)
            % Evaluates the covariance matrix of the distribution
            %
            %   C = cov(obj);
            %
            %       This function applies only when num == 1.
            %       The output C is a covariance matrix of size K x K.
            %
            
            k = obj.K;
            
            if k > 1
                a = obj.alpha;
                if isscalar(a)      % symmetric
                    s = k * a;
                    p = 1 / k;
                    vo = - p^2 / (s + 1);
                    C = adddiag(constmat(k, k, vo), p / (s + 1));
                else                % non-symmetric
                    s = sum(a);
                    p = a / s;
                    C = - (p * p') * (1/(s+1));
                    C = adddiag(C, p ./ (s + 1));
                end
            else
                C = 0;
            end
        end
                
        
        function v = mode(obj)
            % Evaluates the mode(s) of the distributions
            %
            %   v = mode(obj);
            %
            %       v will be an K x n matrix, with v(:,i) corresponding
            %       to the i-th distribution in obj.
            %
            
            k = obj.K;
            n = obj.num;
            
            if k > 1
                alpha_ = obj.alpha;
                
                if size(alpha_, 1) == 1  % symmetric                    
                    if ~all(alpha_ > 0)
                        error('dirichletd:invalidarg', ...
                            'The mode is valid only when all alpha >= 1 and some alpha > 1.');
                    end
                    v = constmat(k, n, 1 / k);
                    
                else                % non-symmetric
                    a = alpha_ - 1;
                    s = sum(a, 1);
                    if ~(all(a(:) >= 0) && all(s > 0))
                        error('dirichlet:invalidarg', ...
                            'The mode is valid only when all alpha >= 1 and some alpha > 1.');
                    end                    
                    
                    if n == 1
                        v = a ./ s;
                    else
                        v = bsxfun(@times, a, 1 ./ s);
                    end
                end
            else
                v = ones(1, n);
            end
        end
        
        
        function v = entropy(obj)
            % Evaluates the entropy of the distribution(s)
            %
            %   v = entropy(obj);
            %
            %       v will be an 1 x n vector, with v(i) corresponding
            %       to the i-th distribution in obj.
            %

            k = obj.K;
            a = obj.alpha;
            
            if ~isempty(obj.logB)
                lnB = obj.logB;
            else
                lnB = dirichletd.calc_logB(k, a);
            end
            
            k = obj.K;
            if size(a, 1) == 1
                v = lnB + k * (a - 1) .* (psi(k * a) - psi(a));
            else
                n = obj.num;
                s = sum(a, 1);
                if n == 1
                    dpsi = psi(s) - psi(a);
                else
                    dpsi = bsxfun(@minus, psi(s), psi(a));
                end
                v = lnB + dot(a - 1, dpsi);
            end
        end
        
    end
    
    
    %% Construction
    
    methods
        
        function obj = dirichletd(K, alpha, op)
            % Constructs a Dirichlet distribution object
            %
            %   obj = dirichletd(K, alpha);
            %
            %       constructs a Dirichlet distribution object comprised
            %       of Dirichlet distribution(s) over the 
            %       (K-1)-dimensional probability simplex, i.e. 
            %       each sample of the distribution is a K-dimensional
            %       vector that sums to 1.
            %
            %       alpha is the parameter of the distribution, which 
            %       can be a K x n matrix or a 1 x n vector (for symmetric
            %       Dirichlet).
            %
            %   obj = dirichletd(K, alpha, 'pre');
            %
            %       Do pre-computation of logB, which might speed-up the
            %       evaluation of various functions, such as entropy and
            %       logpdf.
            %
            
            % verify inputs
            
            if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 1)
                error('dirichletd:invalidarg', 'K should be a positive integer.');
            end
            
            if ~(isfloat(alpha) && isreal(alpha) && ndims(alpha) == 2 && ...
                    (size(alpha,1) == 1 || size(alpha, 1) == K))
                error('dirichletd:invalidarg', ...
                    'alpha should be a real row vector or a real matrix of n rows.');
            end
            
            n = size(alpha, 2);
            
            obj.K = double(K);
            obj.num = n;
            obj.alpha = alpha;
                        
            if nargin >= 3
                if ~(ischar(op) && strcmpi(op, 'pre'))
                    error('gammad:invalidarg', ...
                        'The 4th argument can only be ''pre''.');
                end
                obj.logB = dirichletd.calc_logB(K, alpha);
            end
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
            
            % verify input
            
            k = obj.K;
            m = obj.num;
            
            if ~(isfloat(X) && ndims(X) == 2 && size(X, 1) == k)
                error('dirichletd:invalidarg', ...
                    'X should be a numeric matrix with size(X,1) == K.');
            end
            
            % prepare
            
            a = obj.alpha;
            lnB = obj.logB;
                        
            if nargin >= 3 && ~isempty(si)
                a = a(:, si);
                if ~isempty(lnB)
                    lnB = lnB(1, si);
                end
            end
            
            if isempty(lnB)
                lnB = dirichletd.calc_logB(k, a);
            end
            
            % compute
            
            if size(a, 1) == 1
                V = (a - 1)' * sum(log(X), 1);
            else            
                V = (a - 1)' * log(X);
            end
            
            if m == 1
                L = V - lnB;
            else
                L = bsxfun(@minus, V, lnB.');
            end            
        end
        
        
        function P = pdf(obj, X, si)
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
                P = exp(logpdf(obj, X));
            else
                P = exp(logpdf(obj, X, si));
            end
        end
        
    end
    
    
    %% Sampling
    
    methods
        
        function X = sample(obj, n, i)
            % Samples from the Dirichlet distribution
            %            
            %   X = obj.sample();
            %   X = obj.sample(n);
            %       draws n samples from the Dirichlet distribution. 
            %       (It must be obj.num == 1 for this syntax).
            %
            %       When n is omitted, it is assumed to be 1.
            %
            %   X = obj.sample(n, i);
            %       draws n samples from the i-th Dirichlet distribution.            
            %
            %       When i is a vector, then n can be a vector of the 
            %       same size. In this case, it draws n(j) samples from
            %       the i(j)-th distribution.                     
            %            
            
            if nargin < 2; n = 1; end
            
            k = obj.K;
            a = obj.alpha;
            
            if nargin < 3 || isempty(i)
                if obj.num > 1
                    error('dirichletd:invalidarg', ...
                        'i needs to be given when obj contains multiple distributions.');
                end
                
                X = dirichlet_sample(k, a, n);
                
            else
                if ~(isvector(n) && isnumeric(n))
                    error('dirichletd:invalidarg', ...
                        'n must be a numeric vector');
                end
                if ~(isvector(i) && isnumeric(i))
                    error('dirichletd:invalidarg', ...
                        'i must be a numeric vector');
                end
                if numel(n) ~= numel(i)
                    error('gammad:invalidarg', ...
                        'The sizes of n and i are inconsistent.');
                end
                a = a(:, i);
                m = size(a, 2);
                
                if m == 1
                    X = dirichlet_sample(k, a, n);
                else
                    N = sum(n);
                    X = zeros(k, N);
                    ej = 0;
                    for j = 1 : m
                        sj = ej + 1;
                        ej = ej + n(j);
                        X(:, sj:ej) = dirichlet_sample(k, a(:,j), n(j));
                    end
                end
            end
        end
        
    end
    
    
    %% Private implementation
    
    methods(Static, Access='private')
        
        function v = calc_logB(K, a)
            
            if K == 1
                v = zeros(1, size(a, 2));
            else
                if size(a, 1) == 1
                    v = gammaln(a) * K - gammaln(K * a);
                else                                
                    v = sum(gammaln(a), 1) - gammaln(sum(a, 1));
                end
            end
        end
        
    end
    
    
end



