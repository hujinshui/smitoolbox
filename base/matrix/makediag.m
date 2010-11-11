function R = makediag(v)
% Construct (multiple) diagonal matrices from diagonal elements
%
%   R = makediag(v);
%       Suppose v is a matrix of size d x n, then it creates a d x d x n
%       array, such that R(:,:,i) is diag(v(:,i)) for each i from 1 to n.
%

% Created by Dahua Lin, on April 14, 2010
%

%% main

if ~(isnumeric(v) && ndims(v) == 2)
    error('makediag:invalidarg', ...
        'v should be a numeric matrix.');
end

[d, n] = size(v);

if d == 1
    R = reshape(v, [1 1 n]);
elseif n == 1
    R = diag(v);
else
    R = zeros(d * d, n, class(v));
    i = 1 + (0:d-1) * (d+1);
    R(i, :) = v;
    R = reshape(R, [d d n]);
end

