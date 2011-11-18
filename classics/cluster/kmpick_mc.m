function r = kmpick_mc(X, K, cost, cfunc)
% Progressively pick K samples with maximum cost to existing centers
%
%   r = kmpick_mc(X, K, cost, cfunc);
%   r = kmpick_mc(n, K, cost, cmat);
%
%       progressively pick K new centers with maximum cost.
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
    use_cmat = true;
    
elseif isnumeric(X) && isa(cfunc, 'function_handle')
    n = size(X, 2);
    use_cmat = false;
    
else
    error('kmpick_mc:invalidarg', 'The inputs are invalid.');
end
    
    
%% main

r = zeros(1, K);

for k = 1 : K
    
    % pick a new center
    
    if isempty(cost)
        i = randi(n);
    else
        [mc, i] = max(cost); %#ok<ASGLU>
    end
    r(k) = i;
    
    % update cost
    
    if use_cmat
        cost_i = cfunc(i, :);
    else
        cost_i = cfunc(X(:,i), X);
    end
    
    if isempty(cost)
        cost = cost_i;
    else
        cost = min(cost, cost_i);
    end                        
end
    

