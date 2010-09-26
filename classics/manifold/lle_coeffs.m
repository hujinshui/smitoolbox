function A = lle_coeffs(X, nbs, cons, reg)
% Solves the coefficients of locally linear embedding
%
%   A = lle_coeffs(X, nbs);
%       solves the coefficients of locally linear embedding with respect to 
%       the give neighborhood system.
%
%       X is a d x n matrix giving the sample points by columns. 
%       nbs is a cell array representing the neighborhood system. 
%       nbs{i} is the vector of indices of the neighboring samples of 
%       the i-th sample.
%
%       In the output, A is a n x n sparse matrix, where A(i, :) is
%       the vector of coefficients for constructing X(:,i) using other
%       samples.
%   
%   A = lle_coeffs(X, nbs, cons, reg);
%       solves the locally linear embedding coefficients with constraints and
%       regularization.
%
%       cons can be either empty (indicating no constraint is enforced) or
%       'affine', which requires that the coefficients sum to 1 for each
%       sample.
%
%       reg is the regularization coefficient, which must be non-negative.
%       If it is omitted, no regularization is enforced. 
%

% History
% -------
%   - Created by Dahua Lin, on Nov 26.
%

%% parse and verify input arguments

assert(isfloat(X) && ndims(X) == 2, 'lle_coeffs:invalidarg', ...
    'X should be a numeric matrix.');

[d, n] = size(X);

assert(iscell(nbs) && numel(nbs) == n, 'lle_coeffs:invalidarg', ...
    'nbs should be a cell array with n cells.');

if nargin < 3
    cons = [];
else
    assert(isempty(cons) || (ischar(cons) && strcmp(cons, 'affine')), ...
        'lle_coeffs:invalidarg', 'cons is invalid.');
end

if nargin < 4
    reg = 0;
else
    assert(isfloat(reg) && isreal(reg) && isscalar(reg) && reg >= 0, ...
        'lle_coeffs:invalidarg', 'reg should be a non-negative scalar.');
end


%% main

I = cell(n, 1);
J = cell(n, 1);
V = cell(n, 1);

for i = 1 : n    
    j = nbs{i};
    j = j(:);
    
    I{i} = i(ones(size(j)));
    J{i} = j;                
end


if isempty(cons)    
    for i = 1 : n
        v = linreg(X(:, J{i}), X(:,i), 1, reg);
        V{i} = v;
    end
    
elseif strcmp(cons, 'affine')    
    for i = 1 : n 
        v = affreg(X(:, J{i}), X(:,i), 1, reg);
        V{i} = v;
    end        
end

I = vertcat(I{:});
J = vertcat(J{:});
V = vertcat(V{:});

A = sparse(I, J, V, n, n);


