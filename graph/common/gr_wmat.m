function W = gr_wmat(g, w)
% Compute the edge-weight matrix from a graph
%
%   W = gr_wmat(g);
%   W = gr_wmat(g, w);
%
%       returns an n x n weight (sparse) matrix W, given by
%
%           W(u, v) = w_e   if e = (u, v) is an edge of g
%                   = 0     otherwise.
%
%       If w is omitted, it creates a logical matrix. 
%       
%       Here, g can be a graph struct or degree-bounded graph struct.
%       For the former, w should be a vector of length m, while for 
%       the latter, w should be a matrix of size K x n.
%       
%       In both cases, w can also be a scalar if all edges have the 
%       same weight.
%

% Created by Dahua Lin, on Oct 28, 2011
%

%% verify input

if is_gr(g)
    n = double(g.n);
    if nargin < 2
        w = [];
    else
        m = g.m;
        if ~(isnumeric(w) && (isscalar(w) || (isvector(w) && numel(w) == m)))
            error('gr_wmat:invalidarg', ...
                'w should be a either a numeric scalar or a vector of length m.');
        end
    end
    
    is_gbnd = false;
    
elseif is_gr_bnd(g)
    n = double(g.n);
    if nargin < 2
        w = [];
    else
        K = g.K;        
        if ~(isnumeric(w) && (isscalar(w) || isequal(size(w), [K n])) )
            error('gr_wmat:invalidarg', ...
                'w should be a numeric matrix of size K x n.');
        end
    end
    
    is_gbnd = true;
   
else
    error('gr_wmat:invalidarg', ...
        'g should be a graph struct or degree-bounded graph struct.');    
end


%% main

% preprocess w

if is_gbnd
    
    i = repmat(1:n, g.K, 1);
    j = double(g.nbs);
    
    u = find(g.nbs > 0);
    i = i(u);
    j = j(u);
    
    if ~isempty(w)
        if ~isscalar(w)
            w = w(u);
        end
    end

else
   
    i = double(g.edges(1,:));
    j = double(g.edges(2,:));
    
    if ~isempty(w)        
        if ~isscalar(w)
            if size(w, 1) > 1
                w = w.';
            end
            
            if g.dty == 'u'
                w = [w, w];
            end
        end   
    end
end

% make W

if isempty(w)    
    W = sparse(i, j, true, n, n);        
else
    if ~isa(w, 'double') 
        w = double(w);
    end
    W = sparse(i, j, w, n, n);        
end


