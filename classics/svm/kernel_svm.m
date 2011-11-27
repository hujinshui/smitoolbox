classdef kernel_svm
    % The class that represents a kernel support vector machine
    %
    
    % Created by Dahua Lin, on Apr 7, 2011
    %
    
    properties(GetAccess='public', SetAccess='private')
        
        dim;        % the input space dimension
        ns;         % the number of support vectors
        Xs;         % support vectors [d x ns]
        ys;         % the label of support vectors [1 x ns]
                
        alpha;      % the dual coefficients of the support vectors [1 x ns]
        beta;       % the primal coefficients of the support vectors [1 x ns]
        b;          % the offset value
        
        kerf;       % the kernel function handle        
    end
    
    
    methods
        
        function obj = kernel_svm(Xs, ys, alpha, beta, b, kf)
            % Construct a kernel SVM model
            %
            %   obj = kernel_svm(Xs, ys, a, b, kf)
            %       constructs a kernel SVM model.
            %   
            %       Input arguments:
            %       - Xs:   the support vectors [d x ns]
            %       - ys:   the label of support vectors [1 x ns or ns x 1]
            %       - alpha:  the alpha vector [1 x ns or ns x 1]
            %       - beta:   the beta vector [1 x ns or ns x 1]
            %       - b:    the offset value
            %       - kf:   the kernel function handle
            %
            
            % verify inputs
            
            if ~(isfloat(Xs) && ndims(Xs) == 2 && isreal(Xs))
                error('kernel_svm:invalidarg', ...
                    'Xs should be a 2D real matrix.');
            end
            if ~isa(Xs, 'double'); Xs = double(Xs); end
            [d, n] = size(Xs);
            
            if ~(isfloat(ys) && isvector(ys) && isreal(ys) && length(ys) == n)
                error('kernel_svm:invalidarg', ...
                    'ys should be a real vector of length ns.');
            end
            if ~isa(ys, 'double'); ys = double(ys); end
            if size(ys, 1) > 1; ys = ys.'; end
            
            if ~(isfloat(alpha) && isreal(alpha) && ndims(alpha) == 2 ...
                    && isvector(alpha) && length(alpha) == n)
                error('kernel_svm:invalidarg', ...
                    'alpha should be a real vector of length n.');
            end
            if ~isa(alpha, 'double'); alpha = double(alpha); end
            if size(alpha, 1) > 1; alpha = alpha.'; end
            
            if ~(isfloat(beta) && isreal(beta) && ndims(beta) == 2 ...
                    && isvector(beta) && length(beta) == n)
                error('kernel_svm:invalidarg', ...
                    'beta should be a real vector of length n.');
            end
            if ~isa(beta, 'double'); beta = double(beta); end
            if size(beta, 1) > 1; beta = beta.'; end
            
            if ~(isfloat(b) && isreal(b) && isscalar(b))
                error('kernel_svm:invalidarg', ...
                    'b should be a real scalar.');
            end
            b = double(b);
            
            if ~isa(kf, 'function_handle')
                error('kernel_svm:invalidarg', ...
                    'kf should be a function handle.');
            end
            
            % make object
            
            obj.dim = d;
            obj.ns = n;
            obj.Xs = Xs;
            obj.ys = ys;
            obj.alpha = alpha;
            obj.beta = beta;
            obj.b = b;
            obj.kerf = kf;                        
        end
        
        
        function r = predict(obj, X)
            % Compute the predictor value on new samples
            %
            %   r = obj.predict(X);
            %
            
            kf = obj.kerf;
            
            Kx = kf(obj.Xs, X);  % Kx: ns x n
            r = obj.beta * Kx + obj.b;
        end
                
    end
    
    
    %% Static functions for model construction/training
    
    
    methods(Static)
        
        function obj = from_psol(X, y, beta, b, kf)
            % Construct a kernel SVM from primal solution
            %
            %   obj = kernel_svm.from_sol(X, y, beta, b, kf);
            %
            %   Input arguments:
            %       - X:        all training samples [d x n]
            %       - y:        the labels of all training samples 
            %       - beta:     the solved beta coefficients
            %       - b:        the offset value
            %       - kf:       the kernel function handle
            %
            
            % verify inputs
            
            if ~(isfloat(X) && ndims(X) == 2 && isreal(X))
                error('kernel_svm:from_sol:invalidarg', ...
                    'X should be a 2D real matrix.');
            end
            n = size(X, 2);
            
            if ~(isfloat(y) && isvector(y) && isreal(y) && length(y) == n)
                error('kernel_svm:invalidarg', ...
                    'ys should be a real vector of length ns.');
            end
            if size(y, 1) > 1; y = y.'; end;   
            
            if ~(isfloat(beta) && isreal(beta) && isvector(beta) && length(beta) == n)
                error('kernel_svm:invalidarg', ...
                    'beta should be a real vector of size d x 1.');
            end
            if size(beta, 1) > 1; beta = beta.'; end
            
            if ~(isfloat(b) && isscalar(b) && isreal(b))
                error('kernel_svm:invalidarg', ...
                    'b should be a real scalar.');
            end
            
            if ~isa(kf, 'function_handle')
                error('kernel_svm:invalidarg', ...
                    'kf should be a function handle.');
            end            
            
            si = find(beta ~= 0);
            Xs_ = X(:, si);
            ys_ = y(:, si);           
            bv = beta(:, si);
            av = bv ./ ys_;
            
            obj = kernel_svm(Xs_, ys_, av, bv, b, kf);            
        end
        
        
        function obj = from_dsol(X, y, K, alpha, kf, c)
            % Construct a kernel SVM from dual solution
            %
            %   obj = kernel_svm.from_sol(X, y, K, a, kf, c);
            %
            %   Input arguments:
            %       - X:        all training samples [d x n]
            %       - y:        the labels of all training samples
            %       - K:        the kernel matrix of X [n x n]
            %       - alpha:    the solved alpha vector [n x 1]
            %       - kf:       the kernel function handle
            %       - c:        the upper bound of alpha values
            %
            
            % verify inputs
            
            if ~(isfloat(X) && ndims(X) == 2 && isreal(X))
                error('kernel_svm:from_sol:invalidarg', ...
                    'X should be a 2D real matrix.');
            end
            n = size(X, 2);
            
            if ~(isfloat(y) && isvector(y) && isreal(y) && length(y) == n)
                error('kernel_svm:invalidarg', ...
                    'ys should be a real vector of length ns.');
            end
            if size(y, 1) > 1; y = y.'; end;
            
            if ~(isfloat(K) && isequal(size(K), [n n]) && isreal(K))
                error('kernel_svm:from_sol:invalidarg', ...
                    'K should be a real matrix of size n x n.');
            end
            
            if ~(isfloat(alpha) && isreal(alpha) && isvector(alpha) && length(alpha) == n)
                error('kernel_svm:invalidarg', ...
                    'alpha should be a real vector of length n.');
            end
            
            if ~isa(kf, 'function_handle')
                error('kernel_svm:invalidarg', ...
                    'kf should be a function handle.');
            end
            
            if ~(isfloat(c) && isscalar(c) && isreal(c) && c > 0)
                error('kernel_svm:invalidarg', ...
                    'c should be a real scalar with c > 0.');
            end
               
            % main
            
            % select support vectors
            
            et = 1e-10 * max(abs(alpha));           
            si = find(alpha > et);
            
            Xs_ = X(:, si);
            ys_ = y(si);
            a_ = alpha(si);
            if size(a_, 1) > 1
                a_ = a_.';
            end
            
            % compute b
            
            sj = find(alpha > et & alpha < c - et);
            
            Kb = K(si, sj);
            yb = y(sj);
            bs = 1 ./ yb - (a_ .* ys_) * Kb; 
            b_ = mean(bs);
            
            obj = kernel_svm(Xs_, ys_, a_, a_ .* ys_, b_, kf);            
        end
        
        
        function obj = train(X, y, kf, c, varargin)
            % Train a kernel SVM on data
            %
            %   obj = kernel_svm.train(X, y, kf, c);
            %
            %       trains a kernel SVM model on input data.
            %
            %       Input arguments:
            %       - X:    the matrix of input samples [d x n]
            %       - y:    the labels of the input samples [vector of
            %               length n], the value of y(i) is -1 or 1.
            %       - kf:   the kernel function handle
            %
            %               K = kf(X1, X2);
            %
            %               Suppose X1 and X2 respectively have n1 and n2
            %               columns, then the size of K is n1 x n2, with
            %               K(i,j) be the kernel value of X1(:,i) and
            %               X2(:,j).
            %
            %       - c:    the upper bound of the alpha values.
            %
            %       One can input additional options in form of name/value
            %       pairs:
            %       - 'kermat':  pre-computed kernel matrix on X.
            %                    If omitted, the function uses kf to
            %                    compute it.
            %
            %       - 'verbose': whether to show the progress of training.
            %                    (default = false)
            %
            %       - 'solver':  the function handle to solve the
            %                    qp_problem. (default = '@mstd_solve')
            %
            
            % verify input arguments
            
            if ~(isfloat(X) && ndims(X) == 2 && isreal(X))
                error('kernel_svm:train:invalidarg', ...
                    'X should be a 2D real matrix.');
            end
            n = size(X, 2);
            
            if ~(isfloat(y) && isvector(y) && isreal(y) && length(y) == n)
                error('kernel_svm:train:invalidarg', ...
                    'y should be a real vector of length n.');
            end
            
            if ~isa(kf, 'function_handle')
                error('kernel_svm:train:invalidarg', ...
                    'kf should be a function handle.');
            end
            
            if ~(isfloat(c) && isscalar(c) && isreal(c) && c > 0)
                error('kernel_svm:train:invalidarg', ...
                    'c should be a real positive scalar.');
            end
            
            % parse options
            
            K =[];
            verbose = false;
            solver = @mstd_solve;
            
            if ~isempty(varargin)
                
                onames = varargin(1:2:end);
                ovals = varargin(2:2:end);
                
                if ~(iscellstr(onames) && numel(onames) == numel(ovals))
                    error('kernel_svm:train:invalidarg', ...
                        'The option list is invalid.');
                end
                
                for i = 1 : numel(onames)
                    cn = onames{i};
                    cv = ovals{i};
                    
                    switch lower(cn)
                        case 'kermat'
                            if ~(isfloat(cv) && isreal(cv) && isequal(size(cv), [n n]))
                                error('kernel_svm:train:invalidarg', ...
                                    'kermat should be an n x n real matrix.');
                            end
                            K = cv;
                            
                        case 'verbose'
                            if ~(isscalar(cv))
                                error('kernel_svm:train:invalidarg', ...
                                    'verbose should be a logical scalar.');
                            end
                            verbose = logical(cv);
                            
                        case 'solver'
                            if ~isa(cv, 'function_handle')
                                error('kernel_svm:train:invalidarg', ...
                                    'solver should be a function handle.');
                            end
                            solver = cv;
                            
                        otherwise
                            error('kernel_svm:train:invalidarg', ...
                                'Unsupported option name %s', cn);
                    end                    
                end
            end
                        
            % main            
            
            if verbose
                fprintf('Start training SVM on %d samples ...\n', n);
            end
            
            % prepare kernel matrix                        
            if isempty(K)
                if verbose
                    fprintf('\tcomputing kernel matrix ...\n');
                end
                K = kf(X, X);
                K = (K + K') * 0.5;
            end
            
            % construct problem
            if verbose
                fprintf('\tconstructing QP problem ...\n');
            end
            P = kernel_svm_prob(K, y, c);
            
            % solve the problem
            if verbose
                fprintf('\tsolving QP problem ...\n');
            end
            a = solver(P);
            
            % make model
            if verbose
                fprintf('\tmaking SVM model ...\n');
            end
            obj = kernel_svm.from_dsol(X, y, K, a, kf, c);                            
            
            
            if verbose
                fprintf('SVM training completed.\n');
            end
        end
        
        
    end
    
    
end