function Q = fmm_inferQ(Pi, logliks)
% Infer the Q matrix of finite mixture model
%
%   Q = fmm_inferQ(Pi, logliks);
%       infers the Q matrix of a finite mixture model, given 
%       prior probabilities (Pi), and log-likelihoods. 
%
%       Input arguments:
%       - Pi:       the prior probabilities 
%                   (K x 1 vector, K x n matrix, or empty)
%
%       - logliks:  the table of log-likelihood values (K x n matrix).
%
%       Here, K is the number of components and n is the number of samples.
%

% Created by Dahua Lin, on Sep 28, 2011
%

%% verify input

if ~(isfloat(logliks) && isreal(logliks) && ndims(logliks) == 2)
    error('fmm_inferQ:invalidarg', 'logliks should be a real matrix.');
end
[K, n] = size(logliks);

if ~( isfloat(Pi) && isreal(Pi) && ...
        (isempty(Pi) || isequal(size(Pi), [K 1]) || isequal(size(Pi), [K n])) )    
    error('fmm_inferQ:invalidarg', ...
        'Pi should be either a scalar, a K x 1 vector, or a K x n matrix.');
end

%% main

if isempty(Pi)    
    E = logliks;    
else
    if size(Pi, 2) == n
        E = logliks + log(Pi);
    else
        E = bsxfun(@plus, logliks, log(Pi));
    end
end

Q = nrmexp(E);
    
    
