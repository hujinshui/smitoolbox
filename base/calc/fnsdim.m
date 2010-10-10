function d = fnsdim(A)
% Get the first non-singleton dimension of an array
%
%   d = fnsdim(A);
%       gets the first non-singleton dimension of the array A.
%

% Created by Dahua Lin, on Sep 13, 2010
%

%% main

nd = ndims(A);
if nd == 2    
    if size(A, 1) > 1
        d = 1;
    else
        d = 2;
    end
else
    s = size(A);
    d = find(s > 1, 1);
    if isempty(d)
        d = nd;
    end
end

    
