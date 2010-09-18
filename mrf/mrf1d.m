function W = mrf1d(n, kernel)
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

% Created by Dahua Lin, on Apr 16, 2010
%

%% verify input arguments

if ~(isscalar(n) && isnumeric(n) && n == fix(n) && n >= 1)
    error('mrf1d:invalidarg', 'n should be a positive integer scalar.');
end

if ~(isfloat(kernel) && isreal(kernel) && isvector(kernel))
    error('mrf1d:invalidarg', 'kernel should be a real vector.');
end

h = numel(kernel);
if n <= h
    error('mrf1d:invalidarg', 'n should be greater than the kernel size.');
end

%% main

if h == 1
    I = 1:n-1;
    J = 2:n;
    V = kernel(ones(1,n-1));
        
else
    ns = n * h - h * (1+h) / 2;
    I = zeros(1, ns);
    J = zeros(1, ns);
    V = zeros(1, ns);
    
    i = 0;
    for k = 1 : h
        si = i+1 : i+(n-k);        
        I(si) = 1 : n-k;
        J(si) = 1+k : n;
        V(si) = kernel(k);        
        i = i + (n-k);
    end
end

W = sparse([I J], [J I], [V V], n, n);    

