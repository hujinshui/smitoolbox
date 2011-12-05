function LP = gaussd_logpdf(G, X, ca_cb)
% Evaluate the log-pdf of Gaussian models 
%
%   LP = gaussd_logpdf(G, X);
%   LP = gaussd_logpdf(G, X, {ca, cb})
%
%       computes the log-pdf of the samples given as columns of X, 
%       with respect to the Gaussian distributions in G.
%
%       Inputs:
%       - G:        A gaussd struct
%       - X:        the sample matrix: each column is a column.
%       - ca_cb:    the two constants evalulated by gaussd_const.
%                   If not provided, the function will call 
%
%       Outputs:
%       - LP:       the evaluated result matrix. Suppose K = G.n and
%                   there are n columns in X. Then the size of LP is
%                   K x n. In particular, LP(k, i) is the log-pdf
%                   at X(:,i) w.r.t. the k-th model in G.
%

% Created by Dahua Lin, on Dec 5, 2011
%


%% main

if nargin < 3
    [ca, cb] = gaussd_const(G);
else
    if ~(iscell(ca_cb) && numel(ca_cb) == 2)
        error('gaussd_logpdf:invalidarg', ...
            'The 3rd argument to gaussd_logpdf should be a cell array with two cells.');
    end
    ca = ca_cb{1};
    cb = ca_cb{2};
end
    
D = gaussd_sqmahdist(G, X, ca);

if isscalar(cb)
    LP = D + cb;
else
    LP = bsxfun(@plus, D, cb);
end


