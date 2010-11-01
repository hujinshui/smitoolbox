function T = gr_tree(varargin)
% Constructs a directed tree/forest structure
%
%   Each tree/forest is represented by a struct with the following 
%   fields:
%   - tag:      the tag indicating the type of the struct ('gr_tree')
%   - n:        the number of vertices
%   - ntr:      the number of trees (if the entire graph is connected, K = 1)
%   - m:        the number of edges (n - ntr)
%   - rs:       the roots [1 x ntr int32 zero-based]
%   - ps:       the parents of each node [n x 1 int32 zero-based]
%   - cs:       the concatenated array of children [m x 1 int32 zero-based]
%   - ncs:      the number of childrens of each node [n x 1 int32]
%   - os:       the offsets of each node for cs [n x 1 int32 zero-based]
%
%   T = gr_tree(parents);
%       constructs a tree/forest by specifying the parents of all nodes.
%       
%       Suppose there are n nodes, then parents should be a vector 
%       of length n, where parents(i) is the parent of the i-th node.
%       For root nodes, their parents are zeros.
%

% Created by Dahua Lin, on Oct 31, 2010
%

%% main

if nargin == 1 && isnumeric(varargin{1})
    parents = varargin{1};

    T = from_parents(parents);    
    
else
    error('gr_tree:invalidarg', ...
        'The input arguments are invalid.');
end


%% core functions

function T = from_parents(parents)

if ~(isnumeric(parents) && isvector(parents))
    error('gr_tree:invalidarg', 'parents should be a vector.');
end

n = numel(parents);

if issparse(parents); parents = full(parents); end
if size(parents, 2) > 1; parents = parents.'; end

roots = find(parents == 0);

if isempty(roots)
    error('gr_tree:rterror', 'No root is found in the tree.');
end
ntr = numel(roots);

[sp, cs] = sort(parents); %#ok<ASGLU>
cs = cs(ntr+1:end);

ncs = intcount([1, n], parents).';
os = [0; cumsum(ncs(1:end-1))];

T = struct( ...
    'tag', 'gr_tree', ...
    'dty', 'd', ...
    'n', n, ...
    'ntr', ntr, ...
    'm', n - ntr, ...
    'rs', int32(roots) - 1, ...
    'ps', int32(parents) - 1, ...
    'cs', int32(cs) - 1, ...
    'ncs', int32(ncs), ...
    'os', int32(os));






