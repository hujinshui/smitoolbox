function [L, cutv] = ug_mincut(G, tws)
% Solve Minimum Cut based on undirected graph
%
%   [L, cutv] = ug_mincut(G, tws);
%       solves the Minimum cut (maxflow) problem based on an undirected
%       graph. 
%
%       In input, G represents the inter-node weighted graph.
%       Let n be the number of nodes of the graph. Then G can be either
%       of an n x n affinity matrix, of a cell array in form of
%       {n, I, J, V} to give the edge entries. 
%
%       The argument tws gives the edge weights between nodes and 
%       terminals (source and sink). tws should be a 2 x n matrix,
%       where tws(1,:) gives the weights for source, and tws(2,:)
%       gives the weights for sink.
%
%       Note that G (or V) and tws should have the same value type, which 
%       can  be either of double, single, or int32. These matrices can be
%       either full or sparse.
%
%       In the output, L is a logical vector of size 1 x n. 
%       L(i) == 0 indicates that the i-th node is assigned to the source
%       part, while L(i) == 1 indicates that the j-th node is assigned
%       to the sink part.
%
%       cutv is the value of the minimum cut.
%
%   Remarks
%   -------
%       - This function is based on Boycov-Kolmogorov algorithm (3.01), 
%         which was originally implemented in C++. The source code is
%         incorporated in the private directory for the convenience
%         of re-compilation.
%
%         Web address of the C++ source: http://vision.csd.uwo.ca/code/
%
%       - The caller should ensure that W is a symmetric matrix. In 
%         the program, only upper triangle part is used (i < j).
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 26, 2010
%

%% verify input

if isnumeric(G)    
    if ~(isreal(G) && ndims(G) == 2 && size(G,1) == size(G,2))
        error('ug_mincut:invalidarg', ...
            'W should be a symmetric real matrix.');
    end
    [I, J, V] = find(G);
    n = size(G, 1);
    
elseif iscell(G) && numel(G) == 4
    n = G{1};
    I = G{2};
    J = G{3};
    V = G{4};
    
    if ~(isnumeric(n) && isscalar(n) && n > 0 && n == fix(n))
        error('ug_mincut:invalidarg', 'n should be a positive integer.');
    end
    
    if ~(isvector(I) && isequal(size(I), size(J), size(V)))
        error('ug_mincut:invalidarg', ...
            'I, J, and V should be vectors of the same size.');
    end
end

if ~(isnumeric(tws) && isreal(tws) && isequal(size(tws), [2 n]))
    error('ug_mincut:invalidarg', ...
        'tws should be a real matrix of size 2 x n.');
end


valc = class(V);
if ~(strcmp(valc, 'double') || strcmp(valc, 'single') || strcmp(valc, 'int32'))
    error('ug_mincut:invalidarg', ...
        'The class of weight values should be either double, single, or int32.');
end

if ~isa(tws, valc)
    error('ug_mincut:invalidarg', ...
        'tws should have the same class as G.');
end


%% main

% pre-process I, J, V

se = I < J;
I = I(se);
J = J(se);
V = V(se);

if ~isa(I, 'int32')
    I = int32(I);
end

if ~isa(J, 'int32')
    J = int32(J);
end

if issparse(tws)
    tws = full(tws);
end


% solve

if nargout < 2
    L = mincut_bk_mex(I, J, V, tws);
else
    [L, cutv] = mincut_bk_mex(I, J, V, tws);
end




