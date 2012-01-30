function params = mm_estimate(gm, pri, X, w, Z, samp)
% Estimate mixture model given sample assignment
%
%   params = mm_estimate(gm, pri, X, w, Z);
%   params = mm_estimate(gm, pri, X, w, Z, samp);
%
%       Estimates component model parameters given sample assignment.
%
%       Input arguments:
%       - gm:       the generative model of observations
%       - pri:      the prior model of parameters
%       - X:        the observations
%       - w:        the weights of observations
%       - Z:        the sample assignment, which can be in either of 
%                   the following forms:
%                   - soft assignment matrix of size K x n
%                   - {K, label_vector}.
%       - samp:     whether to do sampling (default = false)
%
%   Remarks
%   -------
%       - This function is to help the implementation of various mixture
%         model estimation, and is not supposed to be directly called by
%         end users.
%

% Created by Dahua Lin, on Dec 27, 2011
%

%% main

if nargin < 6
    samp = false;
end

% form weight matrix

if isnumeric(Z)
    if isempty(w)
        W = Z;
    else
        W = bsxfun(@times, Z, w);
    end
    K = size(W, 1);
else
    K = Z{1};
    L = Z{2};
    
    if K < 10
        W = l2mat(K, L);
    else
        W = l2mat(K, L, 'sparse');
    end
    
    if ~isempty(w)
        W = bsxfun(@times, W, w);
    end
end
    
% estimate

if isempty(pri)
    params = gm.mle(X, W);
else      
    if samp
        params = cell(1, K);
        for k = 1 : K
            cap = gm.capture(X, W(k,:));            
            params{k} = pri.pos_sample(cap, 1);
        end
        params = gm.combine_params(params{:});
    else
        cap = gm.capture(X, W);  
        params = pri.mapest(cap);
    end
end
    
    
