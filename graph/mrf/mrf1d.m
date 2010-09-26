function W = mrf1d(n, kernel, offsets)
% Construct a second-order MRF between nodes with 1D index set
%
%   W = mrf1d(n, kernel);
%       returns the affinity matrix of the constructed MRF among
%       n nodes with 1D index set. 
%
%       kernel is a vector specifying the weights between pairs of
%       nodes with different index distance. Suppose the length
%       of kernel is h. In the output, W is a sparse matrix of 
%       size n x n, where W(i, i+k) and W(i, i-k) are set to kernel(k),
%       when k <= h, otherwise they are set to zeros
%       
%   W = mrf1d(n, kernel, offsets)
%       returns the affinity matrix of the MRF constructed as follows:
%
%       Here, kernel specifies the weights between pairs of nodes with 
%       specified offsets. Concretely, it will set W(i, i+offsets(k))
%       and W(i, i-offsets(k)) to kernel(k).
%
%       Note that if the length of kernel is h, then mrf1d(n, kernel)
%       is equivalent to mrf1d(n, kernel, 1:h);
%

%   History
%   -------
%       - Created by Dahua Lin, on Apr 16, 2010
%       - Modified by Dahua Lin, on Sep 21, 2010
%           - 
%

%% verify input arguments

if ~(isscalar(n) && isnumeric(n) && n == fix(n) && n >= 1)
    error('mrf1d:invalidarg', 'n should be a positive integer scalar.');
end

if ~(isfloat(kernel) && isreal(kernel) && isvector(kernel))
    error('mrf1d:invalidarg', 'kernel should be a real vector.');
end
h = numel(kernel);

if nargin < 3
    offsets = 1 : h;
else
    if ~(isnumeric(offsets) && isequal(size(offsets), size(kernel)) && ...
            all(offsets == fix(offsets)))
        error('mrf1d:invalidarg', ...
            'offsets should be an array of integers with same size of kernel.');
    end
end


%% main

kz = kernel == 0;
if any(kz)
    offsets = offsets(~kz);
    kernel = kernel(~kz);
    h = numel(kernel);
end

if h == 1
    I = 1 : n-offsets;
    J = 1+offsets : n;
    V = kernel(ones(1, n-offsets));
        
else
    ns = n * h - sum(offsets);
    I = zeros(1, ns);
    J = zeros(1, ns);
    V = zeros(1, ns);
    
    i = 0;
    for k = 1 : h
        o = offsets(k);
        
        si = i+1 : i+(n-o);        
        I(si) = 1 : n-o;
        J(si) = 1+o : n;
        V(si) = kernel(k);        
        i = i + (n-o);
    end
end

W = sparse([I J], [J I], [V V], n, n);    

