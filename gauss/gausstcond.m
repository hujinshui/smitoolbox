classdef gausstcond
    % The class to represent a conditional Gaussian model (with transform)
    %
    %   The model is formulated as follows
    %       
    %       y | x ~ N(Ax, Sigma_y) 
    %       
    %   Here, Sigma_y is the conditional covariance matrix.
    %
    
    % Created by Dahua Lin, on Jun 20, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        dx;     % the vector space dimension for x
        dy;     % the vector space dimension for y
        
        A;      % the transform matrix [dy x dx]
        invSy;  % the inverse of Sigma_y
        
        S;      % the gsymat object for A' * invSy * A
        T;      % the matrix for invSy * A
    end
    
    methods
        
        function mdl = gausstcond(A, invSy)
            % Constructs a the Gaussian conditional model
            %
            %   mdl = gausstcond(A, invSy);   
            %       constructs the conditional model.            
            %       Here, A is the transform matrix, and invSy is 
            %       inv(Sigma_y).            
            %
            
            if ~(isobject(invSy) && invSy.n == 1)
                error('gausstcond:invalidarg', ...
                    'invSy should be a symmetric matrix object with single matrix.');
            end
                 
            dy_ = invSy.d;
            
            if ~(isfloat(A) && ndims(A) == 2 && size(A, 1) == dy_)
                error('gausstcond:invalidarg', ...
                    'A should be a numeric matrix with size(A, 1) == dy.');
            end
            
            mdl.dx = size(A, 2);
            mdl.dy = dy_;
            mdl.A = A;
            mdl.invSy = invSy;                                  
            
            Tmat = invSy * A;
            Smat = A' * Tmat;
            mdl.S = gsymat(0.5 * (Smat + Smat'));
            mdl.T = Tmat;
        end
               
        
        function [dc1, dc2] = posupdate(mdl, Y, W)
            % Compute the posterior update from (weighted) observations
            %
            %   [dc1, dc2] = mdl.posupdate(Y);
            %   [dc1, dc2] = mdl.posupdate(Y, W);
            %       computes the posterior update according to Y.
            %
            %       Suppose, the prior distribution of x is a Gaussian
            %       distribution with canonical parameters c1 and c2,
            %       then conditioned on Y, the posterior distribution of
            %       x remains a Gaussian, with canonical parameters
            %       c1 + dc1 and c2 + dc2.
            %
            %       In the input, Y should be a matrix of size d x n,
            %       where n is the number of observed samples.
            %       W can be a matrix of size k x n, where each row of W
            %       gives a group of weights.            
            %       If W is omitted or empty, it means there is a group
            %       of weights, where each sample has a weight 1.
            %
            %       In the output, dc1 is a matrix of size d x k, and
            %       dc2 is a symmetric matrix object, with dc2.d == dx
            %       and dc2.n == k.
            %
            
            dy_ = mdl.dy;
            if ~(isfloat(Y) && ndims(Y) == 2 && size(Y, 1) == dy_)
                error('gausscond:posupdate:invalidarg', ...
                    'Y should be a numeric matrix with size(Y, 1) == dy.');
            end
                     
            if nargin < 3 || isempty(W)
                n = size(Y, 2);
                dc1 = mdl.T' * sum(Y, 2);
                dc2 = n .* mdl.S;
            else
                dc1 = mdl.T' * (Y * W.');
                dc2 = sum(W, 2).' .* mdl.S;
            end            
        end
        
    end
        
end
