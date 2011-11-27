function r = kmpick_rand(X, K, cost, cfunc) 
% Randomly pick K samples whose cost to existing centers are positive
%
%   r = kmpick_rand(X, K, cost, cfunc);
%   r = kmpick_rand(n, K, cost, cmat);
%
%       Randomly pick K samples whose cost to existing centers are 
%       positive. If cost is empty, then simply pick K distinct
%       samples from X.
%
%       Input arguments:
%       - X:        the sample matrix 
%       - n:        the number of samples
%       - K:        the number of new samples to pick
%       - cost:     the minimum matching cost based on existing centers,
%                   which can be empty.
%       - cfunc:    the cost computation function
%       - cmat:     the pairwise cost matrix of size n x n
%
%       Output argument:
%       - r:        the indices of the samples selected to be new centers
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 27, 2010
%

%% verify input

if isscalar(X) && isnumeric(cfunc)
    n = X;
elseif isnumeric(X) && isa(cfunc, 'function_handle')
    n = size(X, 2);    
else
    error('kmpick_rand:invalidarg', 'The inputs are invalid.');
end


%% main

if isempty(cost)
    r = randpick(n, K);
else
    I = find(cost > 0);
    r = randpick(numel(I), K);
    r = I(r);
end

