classdef glincond
    % The class to represent a Gaussian linear conditional model
    %
    %   The model is formulated as follows:
    %
    %       y | (x, a) ~ N(a' x, sigma_y^2);
    %
    %   Here, a is the coefficient vector to be inferred, (x, y) are
    %   pairs observed.
    %
    
    % Created by Dahua Lin, on Jun 20, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        sig2;   % the observation variance
    end
    
    methods
        
        function mdl = glincond(sig2)
            % Constructs a Gaussian linear conditional model
            %
            %   mdl = glincond(sig2);
            %       Here, sig2 is the observation variance
            %
            
            mdl.sig2 = sig2;            
        end
        
        
        function [dc1, dc2] = posupdate(mdl, X, y, W)
            % Compute the posterior update from (weighted) observations
            %
            %   [dc1, dc2] = mdl.posupdate(X, y);
            %   [dc1, dc2] = mdl.posupdate(X, y, W);
            %       computes the posterior update according to X and y.
            %
            %       Suppose, the prior distribution of a is a Gaussian
            %       distribution with canonical parameters c1 and c2,
            %       then conditioned on X and y, the posterior distribution 
            %       of a remains a Gaussian, with canonical parameters
            %       c1 + dc1 and c2 + dc2.
            %
            %       In the input, X should be a matrix of size d x n, with
            %       x_i given by X(:, i), and y should be a row vector of 
            %       size 1 x n, with y_i given by Y(:, i).
            %       where n is the number of observed samples.
            %       W can be a matrix of size k x n, where each row of W
            %       gives a group of weights.            
            %       If W is omitted or empty, it means there is a group
            %       of weights, where each sample has a weight 1.
            %
            %       In the output, dc1 is a matrix of size d x k, and
            %       dc2 is a symmetric matrix object, with dc2.d == d
            %       and dc2.n == k.
            %
            
            % verify input arguments
            
            if ~(isfloat(X) && ndims(X) == 2)
                error('glincond:posupdate:invalidarg', ...
                    'X should be a numeric matrix.');
            end
            
            if ~(isfloat(y) && ndims(y) == 2 && size(y, 1) == 1)
                error('glincond:posupdate:invalidarg', ...
                    'y should be a numeric row vector.');
            end
            
            if size(X, 2) ~= size(y, 2)
                error('glincond:posupdate:invalidarg', ...
                    'The number of columns must agree for X and y.');
            end
            
            % compute
            
            k = 1 / mdl.sig2;
            
            if nargin < 4 || isempty(W)
                if k == 1    
                    dc1 = X * y.';
                    dc2 = gsymat(X * X');
                else
                    dc1 = X * (y.' * k);
                    dc2 = gsymat(k * (X * X'));
                end
            else
                C1 = bsxfun(@times, W, y).';
                C2 = W;
                if k ~= 1
                    C1 = C1 * k;
                    C2 = C2 * k;
                end
                
                dc1 = X * C1;
                
                m = size(W, 1);
                if m == 1
                    S = bsxfun(@times, X, C2) * X';
                    S = 0.5 * (S + S');
                    dc2 = gsymat(S);
                else        
                    d = size(X, 1);
                    S = zeros(d, d, m, class(X(1) * C2(1))); 
                    for i = 1 : m
                        cS = bsxfun(@times, X, C2(i,:)) * X';
                        S(:,:,i) = 0.5 * (cS + cS');                    
                    end
                    dc2 = gsymat(S);
                end
            end                        
        end
                
    end
    
end

