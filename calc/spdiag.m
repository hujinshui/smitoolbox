function M = spdiag(v1, v2)
% Creates a sparse diagonal matrix
%
%   M = spdiag(dv);
%       creates a diagonal matrix with the diagonal elements given by dv.
%       Suppose dv is a vector of length n, then M is an n x n matrix,
%       with M(i,i) equaling dv(i).
%
%   M = spdiag(n, v);
%       creates an n x n diagonal matrix with each diagonal element being
%       v. Here v is a scalar.
%

%   History
%   -------
%       - Created by Dahua Lin, on Apr 17, 2010
%

%% main

if nargin == 1
    dv = v1;
    
    if ~isvector(dv)
        error('spdiag:invalidarg', 'dv should be a vector.');
    end
    
    n = numel(dv);
    
    if size(dv, 1) > 1
        dv = dv.';
    end
    
    M = sparse(1:n, 1:n, dv, n, n);
    
elseif nargin == 2
    
    n = v1;
    v = v2;
    
    if ~isscalar(v)
        error('spdiag:invalidarg', 'v should be a scalar.');
    end
    
    M = sparse(1:n, 1:n, v, n, n);
end
    


