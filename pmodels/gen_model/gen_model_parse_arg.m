function [K, W, g, has_hmap, has_pri] = gen_model_parse_arg(method, n, nh, Z, hmap)
% A function that helps the parsing of arguments for gen_model methods
%
%   [K, W, g, has_pri] = gen_model_parse_arg(method, n, nh, Z, hmap)
%
%   Input arguments:
%   - method:   'mapest' or 'sample'
%   - n:        the number of samples in input sample array
%   - Z:        the Z argument (as g or W, or empty)
%   - nh:       the number of input hyper parameters
%   - hmap:     the hyper-param map
%
%   For omitted args, the caller should input as an empty array.
%
%   Output arguments:
%   - K:        the number of distinct parameters to estimate.
%   - W:        the weights (empty if not weighted)
%   - g:        the cell array of grouped indices (empty if not grouped) 
%   - has_hmap: whether hyper-param-map is used
%   - has_pri:  whether prior is provided (if false, MLE should be
%               performed)

%% main

switch method
    case 'mapest'
        for_sample = 0;
    case 'sample'
        for_sample = 1;
    otherwise
        error('gen_model_parse_arg:invalidarg', ...
            'Invalid method name %s', method);
end

if isempty(Z)
    K = 1;
    W = [];
    g = [];
    
elseif isnumeric(Z)
    W = Z;
    if ~(isfloat(W) && ndims(W) == 2 && size(W, 2) == n)
        error('gen_model:invalidarg', ...
            'W should be a numeric matrix with n columns.');
    end
    K = size(W, 1);
    g = [];
    
elseif iscell(Z)
    g = Z;
    K = numel(g);
    W = [];
    
else
    error('gen_model:invalidarg', ...
        'The 3rd argument to %s_params is invalid.', method);
end

if nh == 0
    has_pri = 0;
    has_hmap = 0;
    if for_sample
        error('gen_model:invalidarg', ...
            'prior parameter is needed for samping.');
    end
else    
    has_pri = 1;
    if isempty(hmap)
        if nh ~= 1
            error('gen_model:invalidarg', ...
                'Mu should be a column vector when umap is omitted.');
        end
        has_hmap = 0;
        
    elseif for_sample && isscalar(umap)
        if ~(nh == 1 && K == 1)
            error(['gen_model:invalidarg', ...
                'multi-sample drawing is only supported ', ...
                'with single posterior model.']);
        end
        has_hmap = 0;
        
    else
        if ~(isnumeric(umap) && isequal(size(umap), [1 K]))
            error('gen_model:invalidarg', ...
                'umap should be a numeric row vector of length K.');
        end
        has_hmap = 1;
    end
end


