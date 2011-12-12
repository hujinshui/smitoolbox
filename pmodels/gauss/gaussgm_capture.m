function [dh, dJ] = gaussgm_capture(X, w, Jx, A)
% Capture the observations for Gaussian generative model
%
%   The Gaussian generative model, with parameter u, is formulated as
%   
%       x ~ N(A * u, Cx),   or x ~ N(u, Cx) if A is identity.
%
%   Here, let d be the dimension of x, and q be that of u.
%
%   Here, the observation can be captured by the conjugate update as
%
%       dh = sum_i w_i (A' * Jx * x_i);
%       dJ = sum_i w_i (A' * Jx * A);
%
%   Here, Jx is the inverse of Cx, i.e. the precision matrix.
%
%   [dh, dJ] = gaussgm_capture(X, w, Jx);
%   [dh, dJ] = gaussgm_capture(X, w, Jx, A);
%
%       Evaluates the conjugate updates from the given data.
%
%       Inputs:
%       - X:        the sample matrix, size: d x n
%       - w:        the weights of the samples, size: m x n, or a scalar.
%                   If m > 1, then multiple updates are to be evaluated,
%                   each corresponding to a group in w.
%       - Jx:       the conditional precision matrix, which can be 
%                   either a pdmat struct, or a scalar.
%       - A:        the transform matrix. If omitted, it is identity.
%
%       Outputs:
%       - dh:       the updates to the potential vector [q x m]
%       - dJ:       a pdmat struct with dJ.n == m and dJ.d == q.
%

% Created by Dahua Lin, on Dec 10, 2011
%

%% verify inputs

if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
    error('gaussgm_capture:invalidarg', 'X should be a real matrix.');
end
[d, n] = size(X);

if isempty(w)
    w = 1;
else
    if ~(isfloat(w) && isreal(w) && ...
            (isscalar(w) || (ndims(w) == 2 && size(w, 2) == n)))
        error('gaussgm_capture:invalidarg', ...
            'w should be a real scalar or a matrix with n columns.');
    end
    m = size(w, 1);
end

if isfloat(Jx) && isreal(Jx) && isscalar(Jx)
    jv = Jx;
    Jsca = 1;
elseif is_pdmat(Jx)
    if ~(Jx.n == 1 && Jx.d == d)
        error('gaussgm_capture:invalidarg', ...
            'Jx should have Jx.n == 1 and Jx.d == d.');
    end
    if Jx.d == 1 || Jx.ty == 's'
        jv = Jx.v;
        Jsca = 1;
    else
        Jsca = 0;
    end
end    

if nargin < 4
    use_A = 0;
else
    if ~(isfloat(A) && isreal(A) && ndims(A) == 2 && size(A,1) == d)
        error('gaussgm_capture:invalidarg', ...
            'A should be a real matrix with size(A,1) == d.');
    end
end

%% main

% compute Jx * X

if Jsca
   
    if isscalar(w)
    else
    end
        
else
    
    
    
end
    
    
    
    


    


