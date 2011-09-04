function S = rand_pdmat(ty, d, n, evrgn)
% Generates a random pdmat struct 
%
%   S = rand_pdmat(ty, d, n, [a, b]);
%       generates a random pdmat struct.
%
%       Input arguments:
%       - ty:       The type of matrix form ('s', 'd', or 'ty')
%       - d:        the dimension (a matrix has size d x d)
%       - n:        the number of matrices in the struct
%       - [a, b]:   the range of eigenvalues. 
%

% Created by Dahua Lin, on Sep 3, 2011
%

%% verify input

if ~(ischar(ty) && isscalar(ty))
    error('rand_pdmat:invalidarg', 'ty should be a char scalar.');
end

if ~(isnumeric(d) && isscalar(d) && d == fix(d) && d >= 1)
    error('rand_pdmat:invalidarg', 'd should be a positive integer.');
end

if ~(isnumeric(n) && isscalar(n) && n == fix(n) && n >= 1)
    error('rand_pdmat:invalidarg', 'n should be a positive integer.');
end

if ~(isfloat(evrgn) && numel(evrgn) == 2 && isreal(evrgn))
    error('rand_pdmat:invalidarg', ...
        'The 4th argument should be a pair of real numbers.');
end

%% main

a = evrgn(1);
b = evrgn(2);

switch ty
    case 's'
        v = a + rand(1, n) * (b-a);
    case 'd'
        v = a + rand(d, n) * (b-a);
    case 'f'        
        if n == 1
            v = make_C(d, a, b);
        else
            v = zeros(d, d, n);
            for i = 1 : n
                v(:,:,i) = make_C(d, a, b);
            end
        end
        
    otherwise
        error('test_pdmat:rterror', 'Unknown form type %s', ty);
end

S = pdmat(ty, d, v);


%% Sub routines

function C = make_C(d, a, b)

if d == 1
    C = a + rand() * (b - a);
else
    r = orth(randn(d, d));
    ev = a + rand(d, 1) * (b-a);
    C = r' * bsxfun(@times, ev, r);
    C = (C + C') * 0.5;
end
