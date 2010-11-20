classdef pca_std
% The class for Standard Principal Component Analysis (PCA)
%   
%   PCA is a model to derive a linear representation in lower 
%   dimensionality that best preserves variational energy
%   of the sample space
%
%   History
%   -------
%       - Created by Dahua Lin, on May 30, 2008
%       - Modified by Dahua Lin, on Jun 6, 2008
%       - Modified by Dahua Lin, on Nov 20, 2010
%

    properties(GetAccess = 'public', SetAccess = 'private')        
        in_dim;             % the dimensionality of input space (d)
        out_dim;            % the dimensionality of output space (q) 
        basis;              % the basis matrix (a d x q orthonormal matrix)
        center;             % the center in the input space
                
        vars;       % the variance along each principal dimension (q x 1)
        residue;    % the residual varianace of discarded subspace (scalar)        
    end
    
    properties(Dependent)        
        principal_var;      % the variance preserved in principal space
        total_var;          % the total variance of entire space (sum(vars) + residue)
        principal_ratio;    % the ratio of variance preserved in principal space
        residue_ratio;      % the ratio of residual variance        
    end
    
    methods
        function v = get.principal_var(obj)
            v = sum(obj.vars);
        end
        
        function v = get.total_var(obj)
            v = sum(obj.vars) + obj.residue;
        end
        
        function r = get.principal_ratio(obj)
            pv = sum(obj.vars);
            r = pv / (pv + obj.residue);
        end
        
        function r = get.residue_ratio(obj)
            pv = sum(obj.vars);
            r = obj.residue / (pv + obj.residue);
        end
    end
    
    
    methods        
        % construction
        
        function obj = pca_std(B, x0, vars, residue)
            % constructs a PCA model with required information
            %
            %   obj = pca_std(B, x0, vars, residue)
            %       constructs and returns a PCA model, with the following
            %       information:
            %           - B:    the basis matrix (d x q)
            %           - x0:   the center vector in input space (d x 1)
            %           - vars: the variances of principal dimensions (q x 1)
            %           - residue: the residue variance
            %
            
            if ~(isfloat(B) && ndims(B) == 2)
                error('pca_std:invalidarg', 'B should be a matrix.');
            end
            
            [d, q] = size(B);
            
            if q > d
                error('pca_std:invalidarg', ...
                    'The dimension of subspace should not exceed that of the input space.');
            end
            
            if ~(isfloat(x0) && isequal(size(x0), [d 1]))
                error('pca_std:invalidarg', ...
                    'x0 should be a numeric vector of size d x 1.');
            end                        
            
            if ~(isfloat(vars) && isequal(size(vars), [q 1]))
                error('pca_std:invalidarg', ...
                    'vars should be a numeric vector of size q x 1.');
            end
            
            if ~(isfloat(residue) && isscalar(residue))
                error('pca_std:invalidarg', ...
                    'residue should be a numeric scalar.');
            end
                                  
        
            obj.in_dim = d;
            obj.out_dim = q; 
            obj.basis = B;
            obj.center = x0;
            
            obj.vars = vars;
            obj.residue = residue;        
        end
        
        
        function sobj = select(obj, inds)
            % creates a PCA model with selected subset of dimensions
            %
            %   sobj = obj.select(inds);
            %       creates a subspace with subset of dimensions specifed
            %       by the indices in inds.
            %
            %   This method overrides select@subspace_base
            %
            %   Remarks
            %       The caller should ensure that there's no repeated 
            %       selection in the inds, otherwise, the variance 
            %       computation may not be correct.
            %                        
            
            svars = obj.vars(inds);
            sres = sum(obj.vars) - sum(svars) + obj.residue;
            
            sobj = pca_std(obj.basis(:, inds), obj.x_anchor, svars, sres);
            
        end
        
    end
            
   
    methods
        % manipulate PCA model
        
        function n = get_truncated_dim(obj, criterion, v)
            % compute the preseved dimension after truncation according to
            % the specified criteria
            %
            %   n = obj.get_truncated_dim('dim', n);            
            %       if n <= out_dim, simply returns n, otherwise, raises an
            %       error.
            %
            %   n = obj.get_truncated_dim('ratio', r);
            %       return n, such that the n leading dimension preserves
            %       at least ratio r of the total variance in the input
            %       space.
            %
            %       If r > principal_ratio, out_dim will be returned, with
            %       a warning issued.
            %
            %   n = obj.get_truncated_dim('tol', r);
            %       return n, such that all dimensions with variance
            %       smaller than r * vars(1) are discarded.
            %
            
            if ~ischar(criterion)
                error('pca_std:get_truncated_dim:invalidarg', ...
                    'criterion should be a striing.');
            end
            
            if ~(isnumeric(v) && isscalar(v))
                error('pca_std:get_truncated_dim:invalidarg', ...
                'The criterion value v should be a numeric scalar.');
            end
            
            switch criterion
                case 'dim'
                    if v <= obj.out_dim
                        n = v;
                    else
                        error('pca_std:get_truncated_dim:invalidarg', ...
                            'The number of dimensions to preserve exceeds out_dim.');
                    end
                    
                case 'ratio'
                    assert(v > 0 && v < 1, 'pca_std:get_truncated_dim:invalidarg', ...
                        'The ratio value should have 0 < v < 1.');
                    
                    pv = sum(obj.vars);
                    tv = pv + obj.residue;
                    
                    if pv >= v * tv
                        n = find(cumsum(obj.vars) >= v * tv, 1);
                    else
                        n = obj.out_dim;
                        
                        warning('pca_std:ratio_too_large', ...
                            'The ratio value exceeds the principal ratio.');
                    end
                    
                case 'tol'
                    assert(v > 0 && v < 1, 'pca_std:get_truncated_dim:invalidarg', ...
                        'The tol-ratio value should have 0 < v < 1.');
                    
                    n = find(obj.vars >= v * obj.vars(1), 1, 'last');     
                    
                otherwise
                    error('pca_std:get_truncated_dim:invalidarg', ...
                        'Unknown truncation criterion %s', criterion);
            end            
            
        end % get_truncated_dim
        
        
        function M = truncated(obj, criterion, v)
            % returns a truncated model based on this model
            %
            %   M = obj.truncated('dim', n);
            %       truncates the model by selecting the first n
            %       dimensions.
            %
            %   M = obj.truncated('ratio', r);
            %       truncates the model by preserving ratio r of 
            %       variance in the principal subspace.
            %
            %   M = obj.truncated('tol', r);
            %       truncates the model by only preserving the 
            %       dimensions of sufficiently large variance,
            %       i.e. variance > r * max(variance)
            %
            
            ns = obj.get_truncated_dim(criterion, v);                
            M = obj.select(1:ns);
            
        end                        
        
        
        function Y = transform(obj, X, s)
            % projects the input vectors into the subspace
            %
            %   Y = obj.transform(X);
            %       transforms the input vectors by projecting them into
            %       the subspace as
            %           y = B' * (x - center)
            %
            %       X should be a d x n matrix with each column being
            %       an input vector, where d == in_dim, and n is the 
            %       number of vectors. And, Y will be a q x n matrix 
            %       with its columns being the output vectors.
            %
            %   Y = obj.transform(X, s);
            %       conducts partial projection, in which the vectors
            %       are projected to a subspace spanned by the subset
            %       of basis indexed by s, as
            %           y = B(:, s)' * (x - center)
            %       
            %       X should be a d x n matrix as in full projection,
            %       while the output Y will be a q' x n matrix, with
            %       q' equaling the number of basis selected by s.
            %
            %       s can be any form of indexing supported by MATLAB,
            %       including integer array or logical array.
            %
                    
            x0 = obj.center;
            if ~all(x0 == 0)
                X = bsxfun(@minus, X, x0);
            end
            
            if nargin == 2
                Y = obj.basis' * X;
            else
                Y = obj.basis(:, s)' * X;
            end            
            
        end
        
        
        function X = reconstruct(obj, Y, s)
            % reconstructs the vector in input space
            %
            %   X = obj.reconstruct(Y);
            %       reconstructs the vectors in the input space from
            %       the transformed vectors using backward transform,
            %       which uses the values in the transformed vectors
            %       as coefficient of linear combination of basis, as
            %           x = B * y + x_anchor
            %
            %       Y should be a q x n matrix with each column being
            %       a transformed vector in the subspace, while in the
            %       output, X will be a d x n matrix. Here, d and q
            %       respectively denote in_dim and out_dim, and n is the
            %       number of vectors.
            %
            %   X = obj.reconstruct(Y, s);
            %       reconstructs the vectors in the input space from
            %       part of the transformed vectors. Only a subset of
            %       basis is utilized, as
            %           x = B(:, s) * y + center
            %
            %       In partial reconstruction, Y should be a q' x n 
            %       matrix, where q' is the number of the dimensions
            %       selected by the indexing array s. While, in the
            %       output, X remains a d x n matrix.
            %
            
            if nargin == 2
                X = obj.basis * Y;
            else
                X = obj.basis(:, s) * Y;
            end
            
            x0 = obj.center;
            if ~all(x0 == 0);
                X = bsxfun(@plus, X, x0);
            end                
            
        end         
        
    end
    
    
    
    
    methods(Static)
        % Training
        
        function obj = train(X, varargin)
            % trains a PCA model based on data
            %
            %   obj = pca_std.train(X);
            %       trains a PCA model based on data in X.
            %
            %       Suppose there are n samples in d-dimensional space,
            %       then X should be a d x n matrix with each column
            %       representing a sample.
            %
            %       The subspace dimension is determined by the rank of X,
            %       which is no larger than min(d, n-1), the maximum
            %       possible rank.
            %
            %   obj = pca_std.train(X, name1, value1, name2, value2, ...);
            %       trains a PCA model with options given in name-value
            %       list.            
            %
            %       The following options are available.
            %
            %       - 'method':     the method to train PCA model, which
            %                       can be either of the following strings:
            %                       - 'auto': automatically decide a proper
            %                                 method. (default)
            %                       - 'std':  the standard method, which
            %                                 solves the eigen-problem of X*X'
            %                       - 'svd':  the method based on SVD
            %                                 decompostion.
            %                       - 'trans': the method based on solves
            %                                  the eigen-problem of X'*X,
            %                                  which is efficient when n < d.
            %
            %       - 'maxdim':     the maximum number of dimensions to 
            %                       preserve. The option value should be 
            %                       a scalar in the range [1, min(d, n-1)].
            %
            %       - 'ratio':      the ratio of variance to preserve in
            %                       the principal subspace, which should be
            %                       a value in the range (0, 1).
            %
            %       - 'tol':        set a criterion that all preserved
            %                       dimensions should have not-too-small
            %                       large variance. 
            %
            %                       The value t should have 0 < t < 1,
            %                       which means that all preserved
            %                       dimensions should have a variance not
            %                       smaller than t * max(variance)
            %
            %       - 'weights':    The weights of samples. If this option
            %                       is specified, then the model will be
            %                       trained on weighted samples. 
            %
            %                       The weights should be a 1 x n row
            %                       vector, with weights(i) being the
            %                       weight of the i-th sample.
            %
            %       - 'center':     The center of the data in the input
            %                       space. 
            %
            %                       It can be either a d x 1 vector, or 0.
            %                       
            %                       If not specified, the sample mean will
            %                       serve as the center.
            %
            %   The subspace dimension is determined as follows:
            %       - If none of the options 'maxdim', 'ratio', and 'tol' 
            %         are specified, then it sets the subspace dimention to
            %         the rank of X.
            %       - If 'ratio' is not specified, but 'maxdim' or 'ratio'
            %         is specified, then the dimension of subspace is
            %         decided to satisfy the conditions.
            %       - If 'ratio' is specified, the dimension is decided as
            %         the smallest dimension that can preserve at least
            %         specified ratio of variance in the principal
            %         subspace, while satisfying the conditions set by
            %         'maxdim' or 'tol' if they are specified.
            %       - If even the maximum dimension satisyfing the
            %         conditions by 'maxdim' and 'tol' cannot preserve the
            %         specified ratio of variance, that maximum dimension
            %         is used, while a warning is issued.
            %
            
            %% parse and verify input arguments
            
            assert(isfloat(X) && ndims(X) == 2, 'pca_std:train:invalidarg', ...
                'X should be a numeric matrix of floatiing point value type.');
            
            [d, n] = size(X);
            dim_ub = min(d, n-1);
            
            
            % default options
            
            method = 'auto';
            maxdim = [];
            ratio = [];
            tol = [];
            weights = [];
            x0 = [];            
            
            if ~isempty(varargin)
                
                names = varargin(1:2:end);
                values = varargin(2:2:end);
                
                nopts = numel(names);
                assert(numel(values) == nopts && iscellstr(names), ...
                    'pca_std:train:invalidsyntax', ...
                    'The name-value list is incorrect.');
                
                for i = 1 : nopts
                    
                    name = names{i};
                    v = values{i};
                    
                    switch name
                        case 'method'
                            if ~ischar(v)
                                error('pca_std:invalidoption', ...
                                'The method option value should be a string.');
                            end
                            
                            switch v
                                case {'auto', 'std', 'svd', 'trans'}
                                    method = v;
                                otherwise
                                    error('pca_std:invalidoption', ...
                                        'Unknown method %s for trainging PCA', v);
                            end
                            
                        case 'maxdim'
                            if ~(isnumeric(v) && isscalar(v) && ...
                                v == fix(v) && v >= 1 && v <= dim_ub)
                                error('pca_std:invalidoption', ...
                                    'The maxdim should be an integer in the range [1, min(d, n-1)].');
                            end                            
                            maxdim = v;
                            
                        case 'ratio'
                            if ~(isnumeric(v) && isscalar(v) && v > 0 && v < 1)                                
                                error('pca_std:invalidoption', ...
                                    'The ratio should be a scalar with 0 < r < 1.');
                            end                            
                            ratio = v;
                            
                        case 'tol'
                            if ~(isnumeric(v) && isscalar(v) && v > 0 && v < 1)
                                error('pca_std:invalidoption', ...
                                    'The tol value should be a scalar with 0 < t < 1.');
                            end                            
                            tol = v;
                            
                        case 'weights'                            
                            if ~(isfloat(v) && ndims(v) == 2 && ...
                                size(v,1) == 1 && size(v,2) == n)
                                error('pca_std:invalidoption', ...
                                    'The weights should be a 1 x n numeric vector.');
                            end                            
                            weights = v;
                                
                        case 'center'
                            if ~(isequal(v, 0) || ( isfloat(v) && ...
                                ndims(v) == 2 && size(v,1) == d && size(v,2) == 1))
                                error('pca_std:invalidoption', ...
                                    'The center should be either 0 or d x 1 numeric vector.');
                            end                            
                            x0 = v;
                            
                        otherwise
                            error('pca_std:invalidoption', ...
                                'Unknown option %s for training PCA', name);
                            
                    end % option-name-switch                    
                
                end % each option
                
            end % has options
                                   
            %% pre-process samples
            
            % shift-center
            if isempty(x0)
                if isempty(weights)
                    x0 = sum(X, 2) / n;
                else
                    x0 = X * (weights' / sum(weights));
                end
            end
            
            if ~all(x0 == 0)
                X = bsxfun(@minus, X, x0);
            end                                                            
            
            % weight samples
            if isempty(weights)
                tw = n;
            else
                tw = sum(weights);
                X = bsxfun(@times, X, weights);
            end
            
            %% perform spectral analysis
            
            % decide method
            if strcmp(method, 'auto')
                if n^2 * (d + n) < d^3
                    method = 'trans';
                else
                    method = 'std';
                end
            end
            
            % compute
            switch method
                case 'std'    
                    M = X * X';
                    M = 0.5 * (M + M');
                    [U, D] = eig(M);
                    clear M;
                    
                    devs = diag(D);
                    clear D;
                    
                case 'trans'
                    M = X' * X;
                    M = 0.5 * (M + M');                    
                    [V, D] = eig(M);
                    clear M;
                    
                    devs = diag(D);
                    clear D;
                    
                    U = X * V;
                    clear V;
                    
                    U = bsxfun(@times, U, 1./sqrt(sum(U.*U, 1)));
                    
                case 'svd'
                    [U, D] = svd(X, 'econ');
                    devs = diag(D) .^ 2;
                    clear D;                    
            end
            
            % re-arrange
            [svars, si] = sort(devs, 1, 'descend');
            svars = svars(1:dim_ub);
            svars(svars < 0) = 0;
            U = U(:, si(1:dim_ub));
            
            tvar = sum(svars);
            
            %% post-process: subspace-dim decision
            
            if isempty(maxdim) && isempty(ratio) && isempty(tol)               
                
                % no condition set, use rank
                subd = find(svars > eps(class(svars)), 1, 'last');
                
            else
                
                subd = dim_ub;
                
                % max dimension
                if ~isempty(maxdim) && subd > maxdim
                    subd = maxdim;
                end
                
                % value tolerance
                if ~isempty(tol) && svars(subd) < tol * svars(1)
                    subd = find(svars >= tol * svars(1), 1, 'last');
                end
                
                % preserve ratio
                if ~isempty(ratio)                    
                    cum_vars = cumsum(svars);                    
                    
                    if cum_vars(subd) >= ratio * tvar
                        subd = find(cum_vars >= ratio * tvar, 1);
                    else
                        warning('pca_std:fail_preserve_variance', ...
                            'The subspace model satisfying the required conditions cannot preserve sufficient variance.');
                    end
                end
                
            end
            
            if subd < dim_ub                
                svars = svars(1:subd);
                U = U(:, 1:subd);
                sres = tvar - sum(svars);
            else
                sres = 0;
            end
            
            svars = svars / tw;
            sres = sres / tw;
                        
            % output
            
            obj = pca_std(U, x0, svars, sres);            
            
        end
        
    end

end 

