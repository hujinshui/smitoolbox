classdef linear_svm
    % The class to implement standard Linar SVM
    %
    
    % Created by Dahua Lin, on Apr 14, 2011
    %
    
    properties(GetAccess='public', SetAccess='private')
        dim;        % the input space dimension
        ns;         % the number of support vectors
        Xs;         % support vectors [d x ns]
        ys;         % the label of support vectors [1 x ns]
        
        w;          % the coefficient vector [d x 1]
        b;          % the offset value [scalar]
    end
    
    
    methods
    
        function obj = linear_svm(Xs, ys, w, b)
            % Constructs a linear SVM model
            %
            %   obj = linear_svm(Xs, ys, w, b);
            %
            %   Input arguments:
            %   - Xs:   the support vectors [d x ns]
            %   - ys:   the label of support vectors [1 x ns or ns x 1]
            %   - w:    the coefficient vector [d x 1]
            %   - b:    the offset [scalar]
            %
            
            % verify inputs
            
            if ~(isfloat(Xs) && ndims(Xs) == 2 && isreal(Xs))
                error('linear_svm:invalidarg', ...
                    'Xs should be a 2D real matrix.');
            end
            if ~isa(Xs, 'double'); Xs = double(Xs); end
            [d, n] = size(Xs);
            
            if ~(isfloat(ys) && isvector(ys) && isreal(ys) && length(ys) == n)
                error('linear_svm:invalidarg', ...
                    'ys should be a real vector of length n.');
            end
            if ~isa(ys, 'double'); ys = double(ys); end
            if size(ys, 1) > 1; ys = ys.'; end
            
            if ~(isfloat(w) && isequal(size(w), [d 1]) && isreal(w))
                error('linear_svm:invalidarg', ...
                    'w should be a real vector of size d x 1.');
            end
            if ~isa(w, 'double'); w = double(w); end
            
            if ~(isfloat(b) && isscalar(b) && isreal(b))
                error('linear_svm:invalidarg', ...
                    'b should be a real scalar.');
            end
            b = double(b);
            
            obj.dim = d;
            obj.ns = n;
            obj.Xs = Xs;
            obj.ys = ys;
            
            obj.w = w;
            obj.b = b;            
        end
        
        
        function r = predict(obj, X)
            % Compute the predictor value on new samples
            %
            %   r = obj.predict(X);
            %
            
            r = obj.w' * X + obj.b;
        end
    end
    
    
    methods(Static)
        
        function obj = from_sol(X, y, sol)
            % Construct a linear SVM from QP solution
            %
            %   obj = kernel_svm.from_sol(X, y, sol);
            %
            %   Input arguments:
            %   - X:    all training samples [d x n]
            %   - y:    the labels of all training samples [length n]
            %   - sol:  the QP solution ((d+1+n) x 1)
            %
            
            % verify inputs
            
            if ~(isfloat(X) && ndims(X) == 2 && isreal(X))
                error('linear_svm:from_sol:invalidarg', ...
                    'X should be a real matrix.');
            end
            [d, n] = size(X);
            
            if ~(isfloat(y) && isvector(y) && length(y) == n && isreal(y))
                error('linear_svm:from_sol:invalidarg', ...
                    'y should be a real vector of length n.');
            end
            if size(y, 1) > 1; y = y.'; end
            
            ds = d + 1 + n;
            if ~(isfloat(sol) && isequal(size(sol), [ds, 1]) && isreal(sol))
                error('linear_svm:from_sol:invalidarg', ...
                    'sol should be a real vector of size (d+1+n) x 1.');
            end
            
            % main
            
            w_ = sol(1:d);
            b_ = sol(d+1);
            
            v = y .* (w_' * X + b_);
            si = v <= (1 + 1e-5);    
            Xs_ = X(:, si);
            ys_ = y(si);
            
            obj = linear_svm(Xs_, ys_, w_, b_);            
        end
        
        
        function obj = train(X, y, c, varargin)
            % Train a linear SVM from data
            %
            %   obj = linear_svm.train(X, y, c);
            %   
            %       trains a linear SVM model on input data.
            %
            %       Input arguments:
            %       - X:    the matrix of input samples [d x n]
            %       - y:    the labels of the input samples [lengt n]
            %               the value of y(i) can be 1 or -1.
            %       - c:    the weight of slack variables.
            %       
            %       One can input additional options in form of name/value
            %       pairs:
            %
            %       - 'verbose': whether to show the progress of training.
            %                    (default = false)
            %
            %       - 'solver':  the function handle to solve the
            %                    qp_problem. (default = '@mstd_qp' with
            %                    algorithm set to 'interior-point-convex').
            %
            
            % verify input arguments
            
            if ~(isfloat(X) && ndims(X) == 2 && isreal(X))
                error('linear_svm:train:invalidarg', ...
                    'X should be a 2D real matrix.');
            end
            n = size(X, 2);
            
            if ~(isfloat(y) && isvector(y) && isreal(y) && length(y) == n)
                error('linear_svm:train:invalidarg', ...
                    'y should be a real vector of length n.');
            end
            if size(y, 1) > 1; y = y.'; end
                        
            if ~(isfloat(c) && isscalar(c) && isreal(c) && c > 0)
                error('linear_svm:train:invalidarg', ...
                    'c should be a real positive scalar.');
            end
            
            % parse options
            
            verbose = false;
            solver = @(P) mstd_qp(P, ...
                optimset('Algorithm', 'interior-point-convex', 'Display', 'off'));
                        
            if ~isempty(varargin)
                
                onames = varargin(1:2:end);
                ovals = varargin(2:2:end);
                
                if ~(iscellstr(onames) && numel(onames) == numel(ovals))
                    error('linear_svm:train:invalidarg', ...
                        'The option list is invalid.');
                end
                
                for i = 1 : numel(onames)
                    cn = onames{i};
                    cv = ovals{i};
                    
                    switch lower(cn)                            
                        case 'verbose'
                            if ~(isscalar(cv))
                                error('linear_svm:train:invalidarg', ...
                                    'verbose should be a logical scalar.');
                            end
                            verbose = logical(cv);
                            
                        case 'solver'
                            if ~isa(cv, 'function_handle')
                                error('linear_svm:train:invalidarg', ...
                                    'solver should be a function handle.');
                            end
                            solver = cv;
                            
                        otherwise
                            error('linear_svm:train:invalidarg', ...
                                'Unsupported option name %s', cn);
                    end                    
                end
            end
            
            % main
            
            if verbose
                fprintf('Start training SVM on %d samples ...\n', n);
            end
            
            
            % construct problem
            if verbose
                fprintf('\tconstructing QP problem ...\n');
            end
            prb = linear_svm_prob(X, y, c);
            
            % solve the problem
            if verbose
                fprintf('\tsolving QP problem ...\n');
            end
            sol = solver(prb);
            
            % make model
            if verbose
                fprintf('\tmaking SVM model ...\n');
            end
            obj = linear_svm.from_sol(X, y, sol);                            
                        
            if verbose
                fprintf('SVM training completed.\n');
            end                        
        end        
    end
    
    
end


