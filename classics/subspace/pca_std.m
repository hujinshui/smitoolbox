classdef pca_std < subspace_base
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
%           - the new class is now derived from the subspace_base class
%

    properties(GetAccess = 'public', SetAccess = 'private')       
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
        
        function obj = pca_std(B, xa, vars, residue)
            % constructs a PCA model with required information
            %
            %   obj = pca_std(B, xa, vars, residue)
            %       constructs and returns a PCA model, with the following
            %       information:
            %           - B:    the basis matrix (d x q)
            %           - xa:   the anchor(center) vector in input space
            %           - vars: the variances of principal dimensions (q x 1)
            %           - residue: the residue variance
            %
            
            q = size(B, 2);
            
            assert(isfloat(vars) && size(vars,1) == q && size(vars,2) == 1, ...
                'pca_std:invalidarg', ...
                'vars should be a numeric vector of size q x 1.');
            
            assert(isfloat(residue) && isscalar(residue), ...
                'pca_std:invalidarg', ...
                'residue should be a numeric scalar.');
            
            assert(all(diff(vars) <= 0), 'pca_std:varsnotsorted', ...
                'The vars should be sorted in descending order.');
                                  
            obj = obj@subspace_base(B, xa);                        
        
            obj.vars = vars;
            obj.residue = residue;        
        end
        
        
        function sobj = select(obj, inds)
            % creates a subspace model with selected subset of dimensions
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
            
            assert(ischar(criterion), 'pca_std:get_truncated_dim:invalidarg', ...
                'criterion should be a striing.');
            
            assert(isnumeric(v) && isscalar(v), ...
                'pca_std:get_truncated_dim:invalidarg', ...
                'The criterion value v should be a numeric scalar.');
            
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
            center = [];            
            
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
                            assert(ischar(v), 'pca_std:invalidoption', ...
                                'The method option value should be a string.');
                            
                            switch v
                                case {'auto', 'std', 'svd', 'trans'}
                                    method = v;
                                otherwise
                                    error('pca_std:invalidoption', ...
                                        'Unknown method %s for trainging PCA', v);
                            end
                            
                        case 'maxdim'
                            assert(isnumeric(v) && isscalar(v) && ...
                                v == fix(v) && v >= 1 && v <= dim_ub, ...
                                'pca_std:invalidoption', ...
                                'The maxdim should be an integer in the range [1, min(d, n-1)].');
                            
                            maxdim = v;
                            
                        case 'ratio'
                            assert(isnumeric(v) && isscalar(v) && v > 0 && v < 1, ...
                                'pca_std:invalidoption', ...
                                'The ratio should be a scalar with 0 < r < 1.');
                            
                            ratio = v;
                            
                        case 'tol'
                            assert(isnumeric(v) && isscalar(v) && v > 0 && v < 1, ...
                                'pca_std:invalidoption', ...
                                'The tol value should be a scalar with 0 < t < 1.');
                            
                            tol = v;
                            
                        case 'weights'                            
                            assert(isfloat(v) && ndims(v) == 2 && ...
                                size(v,1) == 1 && size(v,2) == n, ...
                                'pca_std:invalidoption', ...
                                'The weights should be a 1 x n numeric vector.');
                            
                            assert(all(v >= 0), 'pca_std:invalidoption', ...
                                'The weights should be all non-negative.');
                            
                            weights = v;
                                
                        case 'center'
                            assert(isequal(v, 0) || ( isfloat(center) && ...
                                ndims(v) == 2 && size(v,1) == d && size(v,2) == 1), ...
                                'pca_std:invalidoption', ...
                                'The center should be either 0 or d x 1 numeric vector.');
                            
                            center = v;
                            
                        otherwise
                            error('pca_std:invalidoption', ...
                                'Unknown option %s for training PCA', name);
                            
                    end % option-name-switch                    
                
                end % each option
                
            end % has options
                                   
            %% pre-process samples
            
            % shift-center
            if isempty(center)
                if isempty(weights)
                    center = sum(X, 2) / n;
                else
                    center = X * (weights' / sum(weights));
                end
            end
            
            if ~all(center == 0)
                X = bsxfun(@minus, X, center);
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
            
            obj = pca_std(U, center, svars, sres);            
            
        end
        
    end

end 







