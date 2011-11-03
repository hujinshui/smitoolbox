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
%       Here, w can be either a vector of length m, or a scalar of all
%       edges have the same weight.
%       


% Created by Dahua Lin, on Oct 28, 2011
%

%% verify input

if ~is_gr(g)
    error('gr_wmat:invalidarg', 'g should be a graph struct.');
end

if nargin < 2
    w = [];
else
    m = g.m;
    if ~(isnumeric(w) && (isscalar(w) || (isvector(w) && numel(w) == m)))
        error('gr_wmat:invalidarg', ...
            'w should be a either a numeric scalar or a vector of length m.');
    end
    

end


%% main

% preprocess w

if ~isempty(w)            
    if ~isscalar(w)
        if size(w, 1) > 1 
            w = w.'; 
        end        
        
        if g.dty == 'u'
            w = [w, w];
        end
    end
    
    if ~isa(w, 'double') 
        w = double(w); 
    end
end

% make W

n = double(g.n);
i = double(g.edges(1,:));
j = double(g.edges(2,:));

if isempty(w)    
    W = sparse(i, j, true, n, n);        
else
    W = sparse(i, j, w, n, n);        
end


