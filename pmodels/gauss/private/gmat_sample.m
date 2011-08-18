function X = gmat_sample(cf, mu, C, n, i)
% sample from Gaussian distributon(s)
%
%   X = gsample(cf, mu, C, n);
%   X = gsample(cf, mu, C, n, i);
%
%   Input arguments:
%   - cf:   the covariance form
%   - mu:   the mean vector(s)
%   - C:    the covariance
%   - n:    the number of samples to draw
%   - i:    the indices of models
%
%   Several cases:
%   - i omitted (only one model): draws n samples from the model
%   - otherwise: draws n(k) samples from the i(k)-th model 
%

%% verify input arguments

if ~(ischar(cf) && isscalar(cf))
    error('gsample:invalidarg', 'cf should be a char scalar.');
end

[d, m] = size(mu);

mc = gmat_num(cf, C);
if ~(mc == 1 || mc == m)
    error('gsample:invalidarg', 'The size of C is inconsistent with mu.');
end

if nargin <= 4
    if m >= 2
        error('gsample:invalidarg', ...
            'There should be only one model when indices are omitted.');
    end    
    
else
    if ~isequal(size(i), size(n))
        error('gsample:invalidarg', ...
            'The sizes of i and n are inconsistent.');
    end
    
    if m == 1
        if ~isequal(i, 1)
            error('gsample:invalidarg', ...
                'i can only be 1 or omitted when there is only one model.');
        end
    end
end

%% main

if m == 1    
    X = gsample(cf, mu, C, n);
        
elseif isscalar(i)
    X = gsample(cf, mu(:,i), gmat_sub(cf, C, i), n);
    
else
    N = sum(n);
    X = zeros(d, N, class(mu));
        
    ej = 0;
    for k = 1 : numel(i)
        sj = ej + 1;
        ej = ej + n(k);
        
        if mc == 1
            Ck = C;
        else
            Ck = gmat_sub(cf, C, i(k));
        end            
                        
        X(:, sj:ej) = gsample(cf, mu(:,i(k)), Ck, n(k));        
    end
end
    

