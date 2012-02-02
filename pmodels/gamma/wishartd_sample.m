function W = wishartd_sample(S, df, n)
% Samples from a Wishart distribution
%
%   W = wishartd_sample(d, df);
%       draws a sample from a d-dimensional df-degree standard Wishart
%       distribution (whose scale matrix is an identity matrix).
%
%       The result is returned in form of a d x d matrix.
%
%   W = wishartd_sample(d, df, n);
%       draws n samples from a d-dimensional df-degree standard Wishart
%       distribution.
%
%       The results are returned as an d x d x n array, with each page
%       being a sample.
%
%   W = wishartd_sample(S, df);
%   W = wishartd_sample(S, df, n);
%       draws n samples from an df-degree Wishart distribution with
%       scale matrix S, which should be a pdmat struct (S.n == 1).
%       If n is omitted, it draws one sample
%

%   History
%   -------
%       - Created by Dahua Lin, on Spe 2, 2011
%

%% verify input arguments

if isnumeric(S) && isscalar(S)
    d = S;
    if ~(d == fix(d) && d >= 1)
        error('wishartd_sample:invalidarg', ...
            'd should be a positive integer scalar.');
    end
    use_S = 0;
    
elseif is_pdmat(S)
    if S.n ~= 1
        error('wishartd_sample:invalidarg', ...
            'S should contain only one matrix.');
    end
    d = S.d;
    use_S = 1;
    
else
    error('wishartd_sample:invalidarg', ...
        'The first argument to wishartd_sample is invalid.');
end

if ~(isfloat(df) && isreal(df) && isscalar(df) && df >= d - 1)
    error('wishartd_sample:invalidarg', ...
        'df should be a real scalar with m >= d - 1.');
end

if nargin < 3
    n = 1;
else
    if ~(isnumeric(n) && isscalar(n) && n == fix(n) && n >= 0)
        error('wishartd_sample:invalidarg', ...
            'n should be a non-negative integer scalar.');
    end
end

%% main

if d == 1
    W = randg(df / 2, 1, n) * 2;
    
    if use_S
        W = W * S.v;
    end
    
    if n > 1
        W = reshape(W, 1, 1, n);
    end    
    
else
    gv = 0.5 * (df - (0:d-1).');
    css = randg(gv(:, ones(1,n))) * 2;
    
    if d == 2
        nrms = randn(1, n);
    else
        nrms = randn(d*(d-1)/2, n);
    end
    
    B = wishartd_sample_cimp(css, nrms);
    
    if (~use_S) || (S.ty == 's' || S.ty == 'd')
        
        if n == 1
            W = B * B';
        else
            W = zeros(d, d, n);
            for i = 1 : n
                cB = B(:,:,i);
                W(:,:,i) = cB * cB';
            end
        end
        
        if use_S
            if S.ty == 's'
                W = W * S.v;
            else
                sqv = sqrt(S.v);
                W = W .* (sqv * sqv');
            end            
        end
        
    else
        if d == 2
            L = chol2x2(S.v);
        else
            L = chol(S.v, 'lower'); 
        end
        
        if n == 1
            LB = L * B;
            W = LB * LB';
        else
            W = zeros(d, d, n);
            for i = 1 : n
                cLB = L * B(:,:,i);
                W(:,:,i) = cLB * cLB';
            end
        end
    end
end
    
