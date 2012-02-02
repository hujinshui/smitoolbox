function [M, inds] = l2mat(K, L, varargin)
% L2MAT Creates a binary table from label map
%
%   M = l2mat(K, L, ...);
%       creates a binary table M from a label map L. Here, L is a 
%       1 x n vector, with L(i) giving the label of the i-th sample.
%       K is the maximum label value. 
%       
%       In, for the i-th column, M(k, i) = 1 if k == L(i) otherwise
%       M(k, i) = 0. Each column has at most one entry equaling 1.
%       If L(i) < 1 or L(i) > K, then M(:,i) is a zero column.
%
%       By default, M is a full matrix of double class of size K x n.
%       However, one can specify additional options to change the default
%       behavior. The function supports the following options:
%       - 'logical':    construct M in logical class
%       - 'sparse':     construct a non-sparse matrix
%
%       One can input no option or one or multiple options in a statement.
%       For example, if you want to construct a sparse matrix of
%       logical type, then you can write:
%
%           M = l2mat(K, L, 'sparse', 'logical');
%
%   [M, inds] = l2mat(K, L, ...);
%       additionally returns the indices of the labels that are within
%       the range of [1, K].
%

%   History
%   -------
%       - Created by Dahua Lin, on Oct 31, 2009
%       - Modified by Dahua Lin, on Apr 7, 2010
%           - change name to l2mat
%           - add support of options: logical and sparse.
%

%% verify input arguments

if ~(isnumeric(L) && isvector(L))
    error('l2mat:invalidarg', 'L should be a numeric vector.');
end

n = length(L);

m_sparse = false;
m_logical = false;

if ~isempty(varargin)   
    for i = 1 : length(varargin)
        op = varargin{i};
        if ~ischar(op)
            error('l2mat:invalidopt', 'Each option should be string.');
        end
        if strcmp(op, 'logical')
            m_logical = true;
        elseif strcmp(op, 'sparse')
            m_sparse = true;
        else
            error('l2mat:invalidopt', 'Unknown option %s', op);
        end
    end    
end


%% main

if size(L, 1) > 1  
    L = L.';
end

if any(L < 1 | L > K)
    inds = find(L >= 1 & L <= K);
    L = L(inds);
else
    inds = 1 : n;
end

if m_sparse    
    if m_logical
        M = sparse(L, inds, true, K, n);
    else
        M = sparse(L, inds, 1, K, n);
    end
else
    idx = L + K * (inds - 1);
    
    if m_logical
        M = false(K, n);
        M(idx) = true;
    else
        M = zeros(K, n);
        M(idx) = 1;
    end
end


