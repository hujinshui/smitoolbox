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
        dim;        % the dimension of input space (d)
        pdim;       % the dimension of the principal subspace (q)
        basis;      % the principal basis matrix (a d x q orthonormal matrix)
        center;     % the center in the input space (d x 1 vector, or just 0)
                
        eigvals;    % the eigenvalues for principal dimensions (q x 1)
        
        principal_var;  % the variance preserved in principal subspace
        residue_var;    % the residual variance (not preserved)
        total_var;      % the total variance of the input space
        
    end
    
    properties(Dependent)        
        principal_ratio;    % the ratio of variance preserved in principal space
        residue_ratio;      % the ratio of residual variance        
    end
    
    methods        
        function r = get.principal_ratio(obj)
            r = obj.principal_var / obj.total_var;
        end
        
        function r = get.residue_ratio(obj)            
            r = obj.residue_var / obj.total_var;
        end
    end
    
    
    methods        
        % construction
        
        function obj = pca_std(B, x0, eigvals, residue)
            % constructs a PCA model with required information
            %
            %   obj = pca_std(B, x0, eigvals, residue)
            %       constructs and returns a PCA model, with the following
            %       information:
            %           - B:        the basis matrix (d x q)
            %           - x0:       the center vector in input space (d x 1)
            %           - eigvals:  the eigenvalues of principal dimensions (q x 1)
            %           - residue:  the residue variance
            %
            
            if ~(isfloat(B) && ndims(B) == 2)
                error('pca_std:invalidarg', 'B should be a matrix.');
            end
            
            [d, q] = size(B);
            
            if q > d
                error('pca_std:invalidarg', ...
                    'The dimension of subspace should not exceed that of the input space.');
            end
            
            if ~(isfloat(x0) && (isequal(x0, 0) || isequal(size(x0), [d 1])))
                error('pca_std:invalidarg', ...
                    'x0 should be a numeric vector of size d x 1 or zero.');
            end                        
            
            if ~(isfloat(eigvals) && isequal(size(eigvals), [q 1]))
                error('pca_std:invalidarg', ...
                    'eigvals should be a numeric vector of size q x 1.');
            end
            
            if ~(isfloat(residue) && isscalar(residue))
                error('pca_std:invalidarg', ...
                    'residue should be a numeric scalar.');
            end
                                  
        
            obj.dim = d;
            obj.pdim = q; 
            obj.basis = B;
            obj.center = x0;
            
            obj.eigvals = eigvals;
            
            obj.principal_var = sum(eigvals);
            obj.residue_var = residue;
            obj.total_var = obj.principal_var + obj.residue_var;
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
            
            sevals = obj.eigvals(inds);
            sres = obj.principal_var - sum(sevals) + obj.residue_var;
            
            sobj = pca_std(obj.basis(:, inds), obj.center, sevals, sres);            
        end
        
        
        function sobj = truncate(obj, p)
            % truncates a PCA model to specified dimension
            %
            %   sobj = obj.truncate(p);
            %       truncates the PCA model to p dimension.
            %
            %       Note that p should not exceed obj.pdim;
            %
            
            if p == obj.pdim
                sobj = obj;
            else
                sobj = select(obj, 1:p);
            end                
        end
        
        
    end
            
    
   
    methods
                      
        function Y = transform(obj, X, s)
            % projects the input vectors onto the subspace
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
        
        function obj = estimate(X, pd, varargin)
            % Estimates a PCA model from Data
            %
            %   obj = pca_std.estimate(X, pd, ...);
            %
            %       trains a PCA model based on data in X.
            %
            %       Suppose there are n samples in d-dimensional space,
            %       then X should be a d x n matrix with each column
            %       representing a sample.
            %
            %       The dimension of the principal subspace is determined
            %       in different ways, depending on the 2nd argument pd.
            %       In particular,
            %
            %       - pd can be a positive scalar:
            %
            %           if 0 < pd < 1, then it determines a minimum 
            %           subspace that preserves at least ratio pd of the 
            %           total variance. (e.g if pd = 0.9, then it preserves 
            %           90% of the variance).
            %
            %           if pd >= 1, then the dimension of the principal
            %           subspace is pd. Note that in this case, pd should
            %           have pd <= min(d, n-1), where n is the number of
            %           columns in X.
            %
            %       - pd can be [] (empty) or omitted, then it determines
            %         a subspace that preserves 99.9% of the variance.
            %         (eqivalent to setting pd to 0.999).
            %
            %       - pd can be a function handle that supports the 
            %         usage as follows:
            %
            %           d = pd(eigvals);
            %
            %         It takes as input a sorted (in descending order)
            %         vector of eigenvalues, and returns the dimension
            %         of the principal subspace.
            %       
            %       In addition, one can specify the following options
            %       in name/value pairs (optionally).
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
            %       - 'weights':    The weights of samples. If this option
            %                       is specified, then the model will be
            %                       trained on weighted samples. 
            %
            %                       The weights should be a 1 x n row
            %                       vector, with weights(i) being the
            %                       weight of the i-th sample.
            %
            %       - 'center':     The pre-computed center of the data 
            %                       in the input space. 
            %
            %                       It can be either a d x 1 vector, or 0.
            %                       
            %                       If not specified, the sample mean will
            %                       serve as the center.
            %
            
            %% parse and verify input arguments
            
            % verify X
            if ~(isfloat(X) && ndims(X) == 2)
                error('pca_std:estimate:invalidarg', ...
                    'X should be a numeric matrix.');
            end
            
            [d, n] = size(X);
            
            % verify pd
            
            dim_ub = min(d, n-1);            
            if nargin < 2 || isempty(pd)
                pd = 0.999;
            else
                if (isnumeric(pd) && isscalar(pd) && pd > 0) || ...
                        isa(pd, 'function_handle')
                    
                    if isnumeric(pd) && pd >= 1
                        if ~(pd == fix(pd) && pd <= dim_ub)
                            error('pca_std:estimate:invalidarg', ...
                                'When pd >= 1, pd should be an integer and pd <= min(d, n-1)');
                        end
                    end                    
                else
                    error('pca_std:estimate:invalidarg', ...
                        'The 2nd argument pd is invalid.');
                end
            end

            
            % default options
            
            method = 'auto';
            weights = [];
            x0 = [];       
            
            % parse options
            
            if ~isempty(varargin)
                
                names = varargin(1:2:end);
                values = varargin(2:2:end);
                
                nopts = numel(names);
                if ~(numel(values) == nopts && iscellstr(names))
                    error('pca_std:train:invalidsyntax', ...
                        'The name-value list is invalid.');
                end
                
                for i = 1 : nopts
                    
                    name = lower(names{i});
                    v = values{i};
                    
                    switch name
                        case 'method'
                            if ~(ischar(v) && ismember(v, {'auto','std','svd','trans'}))
                                error('pca_std:invalidoption', ...
                                    'The method option value should be a string.');
                            end
                            method = v;                           
                            
                        case 'weights'                            
                            if ~(isfloat(v) && isvector(v) && numel(v) == n)
                                error('pca_std:invalidoption', ...
                                    'The weights should be a numeric vector of length n.');
                            end         
                            if size(v, 1) > 1; v = v.'; end
                            weights = v;
                                
                        case 'center'
                            if ~(isfloat(v) && ...
                                    (isequal(v, 0) || isequal(size(v), [d 1])))
                                error('pca_std:invalidoption', ...
                                    'The center should be either 0 or d x 1 numeric vector.');
                            end                            
                            x0 = v;
                            
                        otherwise
                            error('pca_std:invalidoption', ...
                                'Unknown option %s for estimating PCA', name);                            
                    end                                    
                end                 
            end 
            
                                   
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
            
            % re-arrange in descending order           
            [sevs, si] = sort(devs, 1, 'descend');                                    
            sevs = max(sevs(1:dim_ub), 0);
            U = U(:, si(1:dim_ub));
            
            sevs = sevs * (1 / tw);
            tvar = sum(sevs);
            
            %% select principal subspace
            
            if isnumeric(pd)
                if pd >= 1
                    p = pd;
                else
                    vr = cumsum(sevs) * (1 / tvar);
                    p = find(vr >= pd, 1);
                    if isempty(p); p = dim_ub; end
                end
            else
                p = pd(sevs);
            end
                                                
            if p < dim_ub                
                pevs = sevs(1:p);
                Up = U(:, 1:p);
                res = tvar - sum(pevs);
            else
                pevs = sevs;
                Up = U;
                res = 0;
            end                        
                        
            % output
            
            obj = pca_std(Up, x0, pevs, res);            
            
        end
        
    end

end 

