function y = pieces(x, edges, dir, alg)
% PIECES Locates the piece that the input values are in
%
%   y = pieces(x, edges);
%   y = pieces(x, edges, 'L');
%   y = pieces(x, edges, 'R');
%
%       determines which bins each value of x is in. Particularly, 
%       let edges be a vector of length m+1 as [e_0, ..., e_m].
%       These edges divide the real line into m+2 bins. 
%
%       If the third argument is 'L' (or omitted), the bins are
%       left-closed, as given by
%       (-inf, e_0), [e_0, e_1), ..., [e_{m-1}, e_m), [e_m, inf).
%
%       If the third argument is 'R', the bins are right-closed, as
%       given by
%       (-inf, e_0], (e_0, e_1], ..., (e_{m-1}, e_m], (e_m, inf).
%
%       The indices of these bins are 0, 1, ..., m+1. If x(i) falls in
%       the bin whose index is k, then y(i) is k.
%
%       x should be a matrix of class double, or single, then y will be
%       a matrix of the same size of class double.
%       
%
%   y = pieces(x, edges, dir, alg);
%       further specifies which algorithm the function should use.
%
%       The value of alg can be either if the following:
%       - 'std':    the standard implementation. The complexity is
%                   O(m n), where m is the number of bins and n the
%                   number of values in x.
%       - 'sorted': the implementation specially for the case when
%                   x is sorted. This should be used only when x is
%                   really sorted. The complexity is O(m + n).
%       - 'auto':   it automatically chooses the algorithm.
%                   In particular, 
%                   if m > 2 * log(n), then it first sort x and then 
%                   use the method designed for sorted values;
%                   otherwise, it uses the standard implementation.
%       
%       If alg is not specified, it uses the default 'auto'.
%   

%   History 
%   -------
%       - Created by Dahua Lin, on June 8, 2010
%

%% verify input

if ~(ndims(x) == 2 && isfloat(x) && ~issparse(x) && isreal(x))
    error('pieces:invalidarg', 'x should be a non-sparse real matrix.');
end

if ~(isvector(edges) && isa(edges, class(edges)) && ~issparse(edges))
    error('pieces:invalidarg', 'edges should be a vector of the same class as x.');
end

if nargin < 3
    is_left = true;
else
    if strcmp(dir, 'L')
        is_left = true;
    elseif strcmp(dir, 'R')
        is_left = false;
    else
        error('pieces:invalidarg', ...
            'dir should be a char with value ''L'' or ''R''.');
    end
end

if nargin < 4
    alg = 'auto';
else
    if ~ischar(alg)
        error('pieces:invalidarg', ...
            'alg should be a string giving the algorithm name.');
    end              
end
        
%% main

switch alg
    
    case 'std'
        y = pieces_cimp(x, edges, is_left, false);
        
    case 'sorted'
        y = pieces_cimp(x, edges, is_left, true);
        
    case 'auto'    
        if numel(edges) > 2 * log(numel(x))
            [sx, si] = sort(x);
            sy = pieces_cimp(sx, edges, is_left, true);
            y = zeros(size(sy));
            y(si) = sy;
        else
            y = pieces_cimp(x, edges, is_left, false);
        end            
end
        

