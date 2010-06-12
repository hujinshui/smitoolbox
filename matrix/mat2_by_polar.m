function C = mat2_by_polar(a, b, theta)
% Create 2x2 positive definite matrices from the polar representation
%
%   C = mat2_by_polar(a, b, theta);
%       creates a 2 x 2 positive definite matrix as 
%
%           C = R * diag([a b]) * R'.
%
%       Here, R is a rotation matrix given by R = rotmat2(theta).
%
%       a, b, and theta can be vectors of n elements, then C is an 
%       array of size 2 x 2 x n, with C(:,:,i) derived from 
%       a(i), b(i), and theta(i).
%

%   History
%   -------
%       - Created by Dahua Lin, on Jun 11, 2010
%

%% main

if ~(isfloat(a) && isfloat(b) && isfloat(theta))
    error('mat2_by_polar:invalidarg', ...
        'a, b, and theta should be numeric scalars or vectors.');
end

if isscalar(a) && isscalar(b) && isscalar(theta)
    c = cos(theta);
    s = sin(theta);
    cc = c * c;
    ss = s * s;
    cs = c * s;
    C = [a * cc + b * ss, (a-b) * cs; (a-b) * cs, a * ss + b * cc];
    
elseif isvector(a) && isequal(size(a), size(b), size(theta))
    c = cos(theta);
    s = sin(theta);
    cc = c .^ 2;
    ss = s .^ 2;
    cs = c .* s;
    
    v1 = a .* cc + b .* ss;
    v2 = (a - b) .* cs;
    v3 = a .* ss + b .* cc;   
    
    C = reshape([v1; v2; v2; v3], [2 2 numel(a)]);
else
    error('mat2_by_polar:invalidarg', ...
        'a, b, and theta should be numeric scalars or vectors of the same size.');
end
   
    
    



