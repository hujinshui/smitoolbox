function G = mgraph(varargin)
% Construct an MGraph struct
%
%   MGraph struct is a common data structure to represent a graph
%   in smi/graph algorithms. The main advantage of using MGraph is
%   its efficiency in interacting with C++ mex implementations.
%
%   A MGraph struct contains following fields:
%   - n:    the number of nodes
%   - I:    the source node indices of all edges (m x 1)
%   - J:    the target node indices of all egdes (m x 1)
%   - W:    the edge weights (m x 1)
%
%   G = mgraph(W);
%       constructs an MGraph struct using edge weight matrix. 
%
%   G = mgraph(n, I, J, W);
%       constructs an MGraph struct using edge information
%
%   G = mgraph(G);
%       verifies whether G is a validity-checked mgraph struct
%

% Created by Dahua Lin, on Oct 1, 2010
%

%% main

if nargin == 1 
    var1 = varargin{1};
    if isnumeric(var1) || islogical(var1)
        W = var1;
        n = size(W,1);
        if ~(ndims(W) == 2 && n == size(W,2))
            error('mgraph:invalidarg', 'W should be a square matrix.');
        end
        [I, J, W] = find(W);
        
        if islogical(W) || isinteger(W)
            W = int32(W);
        end        
        G = struct('n', n, 'I', int32(I)-1, 'J', int32(J)-1, 'W', W, ...
            'mgr_checked', true);
        
    elseif isstruct(var1) && isscalar(var1)
        G = var1;
        
        if ~(isfield(G, 'mgr_checked') && G.mgr_checked)
            error('mgraph:invalidarg', 'G is not a validity-checked mgraph struct.');
        end
    end

elseif nargin == 4
    
    G.n = check_n(varargin{1});
    [G.I, G.J, G.W] = check_IJW(varargin{2}, varargin{3}, varargin{4});    
    G.mgr_checked = true;

else
    error('mgraph:invalidarg', 'The input arguments to mgraph are invalid.');
end


%% sub functions for checking

function n = check_n(n)

if ~(isscalar(n) && n >= 0 && n == fix(n))
    error('mgraph:invalidarg', 'n should be a positive integer.');
end

n = double(n);


function [I, J, W] = check_IJW(I, J, W)

if ~(isvector(I) && isnumeric(I))
    error('mgraph:invalidarg', 'I should be a numeric vector.');
end

if ~(isvector(J) && isnumeric(J))
    error('mgraph:invalidarg', 'J should be a numeric vector.');
end

if ~(isvector(W) && (isnumeric(W) || islogical(W)))
    error('mgraph:invalidarg', 'G.W should be a numeric or logical vector.');
end

if ~isequal(size(I), size(J), size(W))
    error('mgraph:invalidarg', 'The sizes of I, J, and W are not consistent.');
end

if ~isa(I, 'int32')
    I = int32(I);
end

if ~isa(J, 'int32')
    J = int32(J);
end

if islogical(W) || isinteger(W)
    W = int32(W);
end
    
I = I - 1;
J = J - 1;


