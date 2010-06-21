classdef gausscond
    % The class to represent a conditional Gaussian model
    %
    %   The model is formulated as follows
    %       
    %       y | x ~ N(x, Sigma_y) 
    %       
    %   Here, Sigma_y is the conditional covariance matrix.
    %
    
    % Created by Dahua Lin, on Jun 20, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        d;     % the vector space dimension       
        invSy;  % the inverse of Sigma_y
    end
    
    methods
        
        function mdl = gausscond(invSy)
            % Constructs a the Gaussian conditional model
            %
            %   mdl = gausscond(invSy);   
            %       constructs the conditional model.            
            %       Here, invSy is inv(Sigma_y).            
            %
            
            if ~(isobject(invSy) && invSy.n == 1)
                error('gausscond:invalidarg', ...
                    'invSy should be a symmetric matrix object with single matrix.');
            end
                        
            mdl.d = invSy.d;
            mdl.invSy = invSy;                                  
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
            
            d_ = mdl.d;
            if ~(isfloat(Y) && ndims(Y) == 2 && size(Y, 1) == d_)
                error('gausscond:posupdate:invalidarg', ...
                    'Y should be a numeric matrix with size(Y, 1) == d.');
            end
            
            iSy = mdl.invSy;            
            if nargin < 3 || isempty(W)
                n = size(Y, 2);
                dc1 = iSy * sum(Y, 2);
                dc2 = n .* iSy;                                
            else
                dc1 = iSy * (Y * W.');
                dc2 = sum(W, 2).' .* iSy;
            end            
        end
        
    end
        
end

