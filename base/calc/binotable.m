function M = binotable(n)
% Generate a table of binomial coefficients
%
%   M = binotable(n);
%       generates a table M of size (n+1) x (n+1), where
%       when i >= j, M(i, j) = nchoose(i-1, j-1), and
%       when i < j, M(i, j) = 0.
%
%       M is a lower triangle matrix.
%

%   History
%   -------
%       - Created by Dahua Lin, on June 6, 2010
%

%% verify input

if ~(isnumeric(n) && isscalar(n) && n >= 0)
    error('binotable:invalidarg', ...
        'n should be a non-negative integer scalar.');
end


%% main

I = repmat((1:n+1).', 1, n+1);
J = repmat(1:n+1, n+1, 1);

si = I >= J;
I = I(si);
J = J(si);

F = [1, cumprod(1:n)];

V = F(I) ./ (F(J) .* F(I - J + 1));
M = zeros(n+1, n+1);
M(I + (n+1) * (J-1)) = V;

    
