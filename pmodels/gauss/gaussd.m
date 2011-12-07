function G = gaussd(op, param1, param2)
% Constructs a Gauss distribution struct
%
%   G = gaussd('m', mu, C);
%
%       constructs a Gaussian distribution struct with mean and 
%       covariance.
%
%       Suppose you are to construct a struct G comprised of n 
%       Gaussian distributions on a d-dimensional space.
%
%       Inputs:
%       - mu:       the mean vector(s), which can be either a d x n matrix
%                   or just zero (indicating it has zero-mean)
%       - C:        C can be given in either of the following form:
%                   - a scalar s: indicating a covariance s * I
%                   - a d x 1 vector: indicating a diagonal covariance
%                   - a d x d full covariance matrix
%                   - a pdmat struct, with C.d == d. 
%                     Here, C.n can be 1 or n. When C.n == 1, it means
%                     all distributions shared the same covariance.
%
%       Outputs:
%       - G:        a gaussd struct using mean parameters.
%       
%   G = gaussd('c', h, J);
%
%       constructs a Gaussian distributin struct with potential vector
%       and precision matrix.
%
%       Suppose you are to construct a struct G comprised of n 
%       Gaussian distributions on a d-dimensional space.
%
%       Inputs:
%       - h:        the potential vector(s), which can be either a d x n 
%                   matrix or just zero (indicating it has zero-mean)
%       - J:        J can be given in either of the following form:
%                   - a scalar s: indicating a matrix s * I
%                   - a d x 1 vector: indicating a diagonal precision 
%                     matrix
%                   - a d x d full precision matrix
%                   - a pdmat struct, with J.d == d. 
%                     Here, J.n can be 1 or n. When J.n == 1, it means
%                     all distributions shared the same precision matrix.
%
%       Outputs:
%       - G:        a gaussd struct using canonical parameters.
%
%   G = gaussd('m', G0);
%   G = gaussd('c', G0);
%
%       converts the gaussd struct G0 to a specific type.
%   

% Re-created by Dahua Lin, on Dec 5, 2011
%

%% main skeleton

if ischar(op) && isscalar(op)
            
    if op == 'm'
        
        if isnumeric(param1)
            [d, n, mu, C] = verify_args(param1, param2, 'mu', 'C');        
            G.tag = 'gaussd';
            G.ty = 'm';
            G.d = d;
            G.n = n;
            G.mu = mu;
            G.C = C;
        
        elseif is_gaussd(param1)
            ty = param1.ty;
            if ty == 'm'
                G = param1;
            elseif ty == 'c'
                G = cvt_c2m(param1);
            end
            
        else
            error('gaussd:invalidarg', 'The inputs are invalid.');
        end
        
    elseif op == 'c'
        
        if isnumeric(param1)
            [d, n, h, J] = verify_args(param1, param2, 'h', 'J');
            
            G.tag = 'gaussd';
            G.ty = 'c';
            G.d = d;
            G.n = n;
            G.h = h;
            G.J = J;
            
        elseif is_gaussd(param1)
            ty = param1.ty;
            if ty == 'c'
                G = param1;
            elseif ty == 'm'
                G = cvt_m2c(param1);
            end
                        
        else
            error('gaussd:invalidarg', 'The inputs are invalid.');
        end
            
        
    else
        error('gaussd:invalidarg', ...
            'The 1st argument to gaussd can only be either ''m'' or ''c''.');
    end
else
    error('gaussd:invalidarg', ...
        'The 1st argument to gaussd can only be either ''m'' or ''c''.');
end


%% verify and parse inputs

function [d, n, a1, a2] = verify_args(a1, a2, a1_name, a2_name)

if ~(isfloat(a1) && isreal(a1) && ndims(a1) == 2)
    error('%s should be a real matrix.', a1_name);
end

if isfloat(a2) && isreal(a2)
    if isequal(a1, 0)
        a2 = pdmat(a2);
        d = a2.d;
        n = 1;        
    else
        [d, n] = size(a1);
        if isscalar(a2)
            a2 = pdmat('s', d, a2);
        else
            a2 = pdmat(a2);
            if a2.d ~= d
                error('The size of %s and %s is inconsistent.', a1_name, a2_name);
            end            
        end
    end
    
elseif is_pdmat(a2)
    if isequal(a1, 0)
        d = a2.d;
        if a2.n ~= 1
            error('%s.n must be one when %s is a zero scalar.', a2_name, a1_name);
        end
        n = 1;        
    else
        [d, n] = size(a1);
        if ~(a2.d == d && (a2.n == n || a2.n == 1))
            error('The size of %s and %s is inconsistent.', a1_name, a2_name);
        end        
    end    
else
    error('The form of %s is invalid.', a2_name);
end
    

%% conversion functions

function G = cvt_c2m(G0)
% convert from c-type to m-type
%
%   mu = inv(J) * h;
%   C = inv(J);
%

C = pdmat_inv(G0.J);
if isequal(G0.h, 0)
    mu = 0;
else
    mu = pdmat_mvmul(C, G0.h);
end

G.tag = 'gaussd';
G.ty = 'm';
G.d = G0.d;
G.n = G0.n;
G.mu = mu;
G.C = C;


function G = cvt_m2c(G0)
% convert from m-type to c-type
%
%   h = inv(C) * mu;
%   J = inv(C);
%

J = pdmat_inv(G0.C);
if isequal(G0.mu, 0)
    h = 0;
else
    h = pdmat_mvmul(J, G0.mu);
end

G.tag = 'gaussd';
G.ty = 'c';
G.d = G0.d;
G.n = G0.n;
G.h = h;
G.J = J;

