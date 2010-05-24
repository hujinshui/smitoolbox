function [gt, U] = gstransform(g, varargin)
% Compute the sum of transform with respect to a Gaussian matrix
%
%   gt = gstransform(g, As, W, op);
%       computes the sum of the transforms of the Gaussian matrix in g,
%       according to the following formula:
%
%           R_k = \sum_{i=1}^n w_{ki} A_{ki}^T M A_{ki}
%
%       Here, M is the positive definite matrix represented by g, and
%       R_k is the k-th matrix in gt. The weights are provided in the
%       rows of matrix W, and the matrices A_{ki} are given in As.
%
%       If k takes the values in 1,...,m, and i takes values in 1,...,n,
%       then W can be a matrix of size m x n. If m == 1, W can be 
%       specified as an empty array, in which case, w_{ki} always equals 1.
%
%       Suppose, g.dim == d, then each A_{ik} can be a d x q matrix, 
%       and thus, gt.dim == q.
%
%       The matrices A_{ik} are given by the argument As. 
%       Depending on the value of op, As can be in different form:
%
%       If op is 'N',  each A_{ik} is given in a normal matrix form.
%       Then, the size of As can be either of the following:
%           - d x q x n x m:   A_{ik} given by As(:,:,i,k)
%           - d x q x n:       every A_{ik} given by the same As(:,:,i)
%       If op is 'T',  each A_{ik} is given in a transpose matrix form.
%       Then, the size of As can be either of the following:
%           - q x d x n x m:   A_{ik} given by As(:,:,i,k)'
%           - q x d x n:       every A_{ik} given by the same As(:,:,i)'
%       If op is 'R', each A_{ik} is a row vector.
%       Then, the size of As can be either of the following:
%           - q x n x m:       A_{ik} given by As(:,i,k)'
%           - q x n:           every A_{ik} given by the same As(:,:,i)'
%       If op is 'C', each A_{ik} is a column vector.
%       Then, the size of As can be either of the following:
%           - d x n x m:       A_{ik} given by As(:,i,k)
%           - d x n:           every A_{ik} given by the same As(:,:,i)
%
%       Note that op = 'R' only applies to the case with g.dim == 1.
%
%   [gt, U] = gstransform(g, As, Y, W, op);
%       additionally computes the the following quantity:
%
%           u_k = \sum_{i=1}^n w_{ki} A_{ki}^T M y_i
%
%       The vectors x_i is given by X(:,i), and in the output, u_k
%       is given by U(:,k). 
%
%       X should be a numeric matrix of size d x n, and U is a numeric
%       matrix of size q x m.
%
%   Remarks
%   -------
%       - This function is designed specially for posterior updating
%         of Gaussian prior for linear models. The updating can be 
%         done by adding R_k and u_k respectively to the canonical 
%         parameters theta2 and theta1, when g gives the inverse 
%         covariance of the likelihood model.
%

% Created by Dahua Lin, in Apr 5th, 2010
%

%% verify input arguments

error(nargchk(4, 5, nargin));

if nargin == 4
    As = varargin{1};
    W = varargin{2};
    op = varargin{3};   
    solveU = false;
else
    As = varargin{1};
    Y = varargin{2};
    W = varargin{3};
    op = varargin{4};
    solveU = nargout >= 2;
end

assert(isa(g, 'gmat') && g.num == 1, 'gstransform:invalidarg', ...
    'g should be an object of class gmat that contains a single matrix.');

d = g.dim;

assert(ischar(op) && isscalar(op) && any(op == 'NTRC'), ...
    'gstransform:invalidarg', ...
    'op should be either ''N'' or ''T''.');

assert(isfloat(As), 'gstransform:invalidarg', ...
    'As should be a numeric array.');

if op == 'N' || op == 'T'    
    assert(ndims(As) <= 4, 'gstransform:invalidarg', ...
        'The size of As is incorrect.');      
    m1 = size(As, 4);
    n = size(As, 3);
else    
    assert(ndims(As) <= 3, 'gstransform:invalidarg', ...
        'The size of As is incorrect.');        
    n = size(As, 2);
    m1 = size(As, 3);
end

switch op
    case 'N'
        assert(size(As, 1) == d, 'gstransform:invalidarg', ...
            'The size of As is incorrect.');
        q = size(As, 2);
    case 'T'
        assert(size(As, 2) == d, 'gstransform:invalidarg', ...
            'The size of As is incorrect.');
        q = size(As, 1);
    case 'R'       
        assert(d == 1, 'gstransform:invalidarg', ...
            'op == ''R'' only applies when the g.dim == 1.');           
        q = size(As, 1);
    case 'C'
        assert(size(As, 1) == d, 'gstransform:invalidarg', ...
            'The size of As is incorrect.');   
        q = 1;
end

if m1 > 1
    assert(~isempty(W), 'gstransform:invalidarg', ...
        'W cannot be empty when m > 1.');
end

if ~isempty(W)
    assert(isfloat(W) && ndims(W) == 2, 'gstransform:invalidarg', ...
        'W should be either empty or a numeric matrix.');
    
    m2 = size(W, 1);    
    assert((m1 == 1 || m2 == m1) && size(W, 2) == n, ...
        'gtransform:invalidarg', 'The size of W is illegal.');
    
    m = max(m1, m2);
else
    m = 1;
end

%% main

if d == 1    
    s2 = fullform(g);
    s2 = s2(1);
    
    if q == 1
        A = reshape(As, n, m1);
        A2 = A .^ 2;
        if isempty(W)
            R = s2 * sum(A2, 1);
        else
            if m1 == 1 || m2 == 1
                R = s2 * (W * A2);
            else 
                R = s2 * sum(A2 .* W.', 1);            
            end
        end
        
        if solveU
            if isempty(W)
                U = s2 * (Y * A);
            elseif m2 == 1
                U = s2 * ((W .* Y) * A);
            elseif m1 == 1
                U = s2 * (W * (Y.' .* A)).';
            else
                U = s2 * sum(A .* bsxfun(@times, W, Y).', 1);
            end               
        end
        
    else
        if op ~= 'R'
            A = reshape(As, q, n, m1);
        else
            A = As;
        end
        
        if m == 1
            if isempty(W)
                R = s2 * (A * A');
                
                if solveU
                    U = s2 * (A * Y');
                end
            else
                R = s2 * (bsxfun(@times, A, W) * A');
                
                if solveU
                    U = s2 * (A * (W .* Y)');
                end
            end
                        
        else
            R = zero_R(q, m, As, W);
            
            if solveU
                U = zero_U(q, m, As, W, Y);
                WY = bsxfun(@times, W, Y)';
            end
            
            for k = 1 : m
                if m1 == 1
                    Ak = A;
                else
                    Ak = A(:,:,k);
                end
                
                R(:,:,k) = s2 * (bsxfun(@times, Ak, W(k,:)) * Ak');
                
                if solveU
                    U(:,k) = s2 * (Ak * WY(:,k));
                end                
            end
        end                
    end  
    
else  % d > 1
            
    if q == 1
        A = reshape(As, d, n * m1);
                
        V = reshape(sum(A .* (g * A), 1), n, m1);
        
        if isempty(W)
            R = sum(V, 1);
        elseif m1 == 1 || m2 == 1
            R = W * V;
        else
            R = sum(W.' .* V, 1);
        end
        
        if solveU
            H = g * Y;
            if m1 == 1
                F = sum(A .* H, 1).';
            else
                F = sum(bsxfun(@times, reshape(A, [d n m1]), H), 1);
                F = reshape(F, n, m1);
            end
            
            if isempty(W)
                U = sum(F, 1);
            elseif m2 == 1
                U = W * F;
            elseif m1 == 1
                U = (W * F).';
            else
                U = sum(W.' .* F, 1);
            end
        end        
        
    else
        
        if op == 'T'
            I = reshape(1:q*d, q, d).';
            As = reshape(As, q * d, n, m1);
            As = reshape(As(I(:), :, :), [d q n m1]);
        end
                
        if m1 == 1                      
            R = mmsum(As, reshape(g * As(:,:), [d q n]), W, 'TN');
            
            if solveU
                U = mvsum(As, g * Y, W, 'T');
            end
        else
            R = zero_R(q, m, As, W);            
            if solveU                                                        
                U = zero_U(q, m, As, W, Y);
                gY = g * Y;
            end
            
            for k = 1 : m
                Ak = As(:,:,:,k);
                R(:,:,k) = mmsum(Ak, reshape(g * Ak(:,:), [d q n]), W(k,:), 'TN');
                if solveU
                    U(:,k) = mvsum(Ak, gY, W(k,:), 'T');
                end
            end                            
        end
                            
    end    
end


% output

if q == 1
    gt = gmat_iso(q, R(:).');
else
    if m == 1
        R = 0.5 * (R + R');
    else
        for k = 1 : m
            cR = R(:,:,k);
            R(:,:,k) = 0.5 * (cR + cR');
        end
    end
    gt = gmat_comp(R);
end


%% Auxiliary functions

function U = zero_U(q, m, As, W, Y)

if isempty(W)
    U = zeros(q, m, class(As(1) * Y(1)));
else
    U = zeros(q, m, class(As(1) * W(1) * Y(1)));
end


function R = zero_R(q, m, As, W)

if isempty(W)
    R = zeros(q, q, m, class(As(1)));
else
    R = zeros(q, q, m, class(As(1) * W(1)));
end


