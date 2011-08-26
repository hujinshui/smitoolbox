classdef linear_gauss_genmodel < gen_model
    % The class to implement a Linear Gaussian generative model 
    %
    % The model is formulated as follows
    %
    %   x ~ N(B * mu, Cx) or N(mu, Cx) when B = 1
    %   y ~ N(A * x,  Cy) or N(x,  Cy) when A = 1
    %
    %   parameters:  
    %       - x:        [xdim x 1 double vector]
    %
    %   hyper-params:  
    %       - mu:       [udim x 1 double vector]
    %
    %   design-params:  
    %       - Cx:       [pdmat struct with d = xdim and n = 1]
    %       - Cy:       [pdmat struct with d = ydim and n = 1]
    %       - A:        [ydim x xdim double matrix or a scalar]
    %       - B:        [xdim x udim double matrix or a scalar]          
    %
    
    % Created by Dahua Lin, on Aug 26, 2011
    %
    
    %% properties
    
    properties(GetAccess='public', SetAccess='private')
        
        udim;       % the dimension of mu's space
        xdim;       % the dimension of parameter (x)'s space 
        ydim;       % the dimension of product (y)'s space
                
        Cx;         % the covariance in generating x
        Cy;         % the covariance in generating y
        A;          % the x --> y transform matrix or scalar
        B;          % the mu --> x transform matrix or scalar
        
        Gx;         % a Gaussian distribution: N(0, Cx)
        Gy;         % a Gaussian distribution: N(0, Cy)
    end
    
    %% Construction
    
    methods
        
        function model = linear_gauss_genmodel(Cx, Cy, A, B)
            % Constructs a linear Gaussian generative model
            %
            %   model = linear_gauss_genmodel(Cx, Cy);
            %   model = linear_gauss_genmodel(Cx, Cy, A);
            %   model = linear_gauss_genmodel(Cx, Cy, A, B);
            %
            %       constructs a linear Gaussian generative model given
            %       the design parameters.
            %
            %       Note A or B can be omitted, in which case, it is
            %       assumed to be identity (i.e. scalar 1).
            %
            
            if ~is_pdmat(Cx)
                error('linear_gauss_genmodel:invalidarg', ...
                    'Cx should be a pdmat struct.');
            end
            
            if ~is_pdmat(Cy)
                error('linear_gauss_genmodel:invalidarg', ...
                    'Cy should be a pdmat struct.');
            end
                        
            xd = Cx.dim;
            yd = Cy.dim;
            
            if ~isa(Cx.v, 'double'); Cx.v = double(Cx.v); end
            if ~isa(Cy.v, 'double'); Cy.v = double(Cy.v); end
            
            if nargin < 3
                A = 1;                
            else
                if ~( isfloat(A) && ...
                        (isscalar(A) || isequal(size(A), [yd, xd])) )
                    
                    error('linear_gauss_genmodel:invalidarg', ...
                        ['A should be either ', ...
                        'a numeric matrix of size ydim x dim, ', ...
                        'or a scalar.']);
                end
                
                if ~isa(A, 'double'); A = double(A); end
            end
            
            if isscalar(A) && (xd ~= yd)
                error('linear_gauss_genmodel:invalidarg', ...
                    'A cannot be a scalar when xdim ~= ydim.');
            end            
            
            if nargin < 4
                B = 1;
                ud = xd;
            else
                if ~( isfloat(B) && ...
                        (isscalar(B) || (ndims(B)==2 && size(B,1)==xd)) )
                    
                    error('linear_gauss_genmodel:invalidarg', ...
                        ['B should be either ', ...
                        'a numeric matrix of size ydim x dim, ', ...
                        'or a scalar.']);
                end                
                ud = size(B, 2);
                if ~isa(B, 'double'); B = double(B); end
            end
            
            Gx_ = gaussd.from_mp(0, Cx, 'ip');
            Gy_ = gaussd.from_mp(0, Cy, 'ip');
            
            model.udim = ud;
            model.xdim = xd;
            model.ydim = yd;
            model.Cx = Cx;
            model.Cy = Cy;
            model.A = A;
            model.B = B;
            model.Gx = Gx_;
            model.Gy = Gy_;
        end
        
    end
    
    %% Evaluation methods
    
    methods
            
        function lpri = logpri(model, mu, X)
            % Evaluates the log-prior of given parameters
            %
            %   lpri = model.logpri(mu, X);
            %
            %       Evaluates the log-prior of the parameters in X 
            %       with respect to the hyper-parameter mu.      
            %
            
            % verify inputs
            
            ud = model.udim;
            xd = model.xdim;
            
            if ~(isfloat(mu) && isequal(size(mu), [ud, 1]))
                error('linear_gauss_genmodel:invalidarg', ...
                    'mu should be a numeric vector of size udim x 1.');
            end
            if ~isa(mu, 'double'); mu = double(mu); end
            
            if ~(isfloat(X) && ndims(X) == 2 && size(X,1) == xd)
                error('linear_gauss_genmodel:invalidarg', ...
                    'X should be a numeric matrix with size(X,1) == xdim.');
            end
            if ~isa(X, 'double'); X = double(X); end
            
            % compute
            
            B_ = model.B;
            if isequal(B_, 1)
                bu = mu;
            else
                bu = B_ * mu;
            end
            
            if all(bu == 0)
                Xmu = X;
            else
                Xmu = bsxfun(@minus, X, bu);
            end                
            
            lpri = model.Gx.logpdf(Xmu);            
        end
        
        function Llik = loglik(model, X, Y)
            % Evaluates the generative log-likelihood 
            %
            %   Llik = model.loglik(X, Y);
            %
            %       Evaluates the log-likelihood of Y with respect
            %       to the model with parameters given by X, in a
            %       pairwise manner.
            %
            
            % verify inputs
            
            xd = model.xdim;
            yd = model.ydim;
            
            if ~(isfloat(X) && ndims(X) == 1 && size(X,1) == xd)
                error('linear_gauss_genmodel:invalidarg', ...
                    'X should be a numeric matrix with size(X,1) == xdim.');
            end
            if ~isa(X, 'double'); X = double(X); end
            
            if ~(isfloat(Y) && ndims(Y) == 2 && size(Y,1) == yd)
                error('linear_gauss_genmodel:invalidarg', ...
                    'Y should be a numeric matrix with size(Y,1) == ydim.');
            end
            if ~isa(Y, 'double'); Y = double(Y); end
            
            % compute
            
            A_ = model.A;
            if isequal(A_, 1)
                AX = X;
            else
                AX = A_ * X;
            end
            
            m = size(AX, 1);
            
            if m == 1
                gy0 = model.Gy;
                Llik = gy0.logpdf(bsxfun(@minus, Y, AX));
            else
                gy = gaussd.from_mp(AX, model.Cy, 'ip');
                Llik = gy.logpdf(Y);
            end            
        end
        
    end
    
    
    %% Sampling & Inference methods
    
    methods
    
        function Y = generate(model, X, n, g)
            % Generates product samples given parameters and labels
            %
            %   Y = model.generate(x, n);
            %
            %       generates n sample from the model with parameter x.
            %
            %   Y = model.generate(X, n, g);
            %
            %       generates samples from multi-models with parameters 
            %       X, based on the grouping given in g, a cell array of 
            %       index vectors.
            %
            %       In particular, the samples with indices in g{k} should
            %       be generated by the k-th model (X(:,k)), and n is the 
            %       total number of samples.
            %
            
            % verify inputs
            
            xd = model.xdim;            
            if ~(isfloat(X) && ndims(X) == 1 && size(X,1) == xd)
                error('linear_gauss_genmodel:invalidarg', ...
                    'X should be a numeric matrix with size(X,1) == xdim.');
            end
            if ~isa(X, 'double'); X = double(X); end
            
            if ~(isnumeric(n) && isscalar(n) && n == fix(n) && n >= 0)
                error('linear_gauss_genmodel:invalidarg', ...
                    'n should be a non-negative integer scalar.');
            end
            
            if nargin < 4 || isempty(g)
                if size(X, 2) ~= 1
                    error('linear_gauss_genmodel:invalidarg', ...
                        'X can have one column when no grouping is given.');
                end
                xgrp = 0;
            else
                if ~(iscell(g) && numel(g) == size(X, 2))
                    error('linear_gauss_genmodel:invalidarg', ...
                        ['g should be a cell array of index vectors', ...
                        'with K cells (K = size(X,2)).']);
                end
                xgrp = 1;
            end
               
            % compute
            
            A_ = model.A;
            if isequal(A_, 1)
                AX = X;
            else
                AX = A_ * X;
            end
                        
            Y = model.Gy.sample(n);
            
            if ~xgrp                
                Y = bsxfun(@plus, Y, AX);                
            else
                for k = 1 : numel(g)
                    gk = g{k};
                    Y(:, gk) = bsxfun(@plus, Y(:, gk), AX(:,k));
                end
            end
            
        end
        
        
        function X = map_estimate_params(model, Y, Z, Mu, umap)
            % Performs Maximum-a-posteriori (MAP) estimation of parameters
            %
            %   x = model.map_estimate_params(Y);
            %   x = model.map_estimate_params(Y, [], mu);
            %
            %       Performs MAP estimation of the parameter based on 
            %       a given set of samples, each with weight 1. 
            %
            %       Here, mu is the prior mean of x, and when omitted, 
            %       MLE instead of MAP is performed. mu can be input in
            %       two ways, a udim x 1 vector, or a cell containing
            %       it (for conformance with gen_model interface).
            %
            %   X = model.map_estimate_params(Y, W, alpha);
            %
            %       Performs MAP estimation based on weighted samples.
            %
            %       W should be a K x n matrix, where W(k, i) is the
            %       contribution weight of the i-th sample to the k-th 
            %       model.
            %
            %       In output, X is an xdim x K matrix, with X(:,k) 
            %       estimated based on W(k,:).
            %
            %   X = model.map_estimate_params(Y, g, alpha);
            %
            %       Performs MAP estimation based on grouped samples.
            %       Here, g is a cell array of index vectors, and X(:,k)
            %       should be estimated based on Y(:, g{k});
            %
            %   X = model.map_estimate_params(Y, .., Mu, umap);
            %
            %       Performs MAP estimation with prior-map (umap).
            %       
            %       Here, Mu, in form of a matrix of size udim x K, 
            %       or a cell array with each cell containing a vector,
            %       gives K distinct prior mean vectors. 
            %
            %       In estimation, the i-th parameter X(:,i) should be
            %       estimated based on the umap(i)-th prior mean in Mu.
            %
            
            
        end
        
        
        function X = sample_params(Y, Z, mu, umap)
            % Samples from the posterior distribution of parameters
            %
            %   x = model.sample_params(Y, [], mu);
            %
            %       samples a parameter x, given a set of data Y, and
            %       the prior mean mu.
            %
            %   X = model.sample_params(Y, W, mu);
            %       
            %       samples a parameter (or a set of parameters) based 
            %       on the weighted observations given in Y. 
            %
            %       W should be a K x n matrix, where W(k, i) is the
            %       contribution weight of the i-th sample to the k-th 
            %       model.
            %
            %       In output, X is an xdim x K matrix, with X(:,k) 
            %       estimated based on W(k,:).
            %
            %   X = model.sample_params(Y, [], mu, n);
            %   X = model.sample_params(Y, w, mu, n);
            %
            %       samples n parameters, given a set of data Y, and the 
            %       prior mean mu.
            %
            %   X = model.sample_params(Y, g, mu);
            %
            %       samples a set of parameters based on grouped 
            %       observations given by Y and g. 
            %
            %       Let g be a cell array with K cells, then X
            %       is an xdim x K matrix, with X(:,k) estimated
            %       based on Y(:, g{k}).
            %
            %   X = model.sample_params(Y, .., Mu, umap);
            %
            %       samples K parameters, given a weighted/grouped set 
            %       of data Y, using a prior-map (umap).
            %
            %       Specifically, X(:,k) should be estimated based on
            %       the prior mean Mu(:, umap(k)).
            %
            %
            
        end
        
        
    end
    
    
end





