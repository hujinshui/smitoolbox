function devcheck(title, x1, x2, tol)
% Check whether two arrays are close enough, if not raise a warning
%
%   devcheck(title, x1, x2, tol);
%

assert(isfloat(x1));
assert(isfloat(x2));
assert(isequal(size(x1), size(x2)));

maxdev = max(abs(x1(:) - x2(:)));

if maxdev > tol
    warning('devcheck:largedev', ...
        'Large deviation encountered in verifying %s (dev = %g)', ...
        title, maxdev);
end
