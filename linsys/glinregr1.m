classdef glinregr1 < gcondmodel
    % The class representing a Gaussian linear regression model on 
    % 1D outputs
    %
    %   In this model
    %       y_i = c' * x_i + e
    %
    %   where e follows a normal law with variance sigma^2.
    %
    
    % Created by Dahua Lin, on April 14, 2010
    %
    
    methods
        
        function obj = glinregr1(cdim)
            % Constructs a Gaussian linear regression model
            %
            %   obj = glinregr1(cdim);
            %       constructs a linear regression model with
            %       input dimension cdim.
            %
            
            if ~(isnumeric(cdim) && isscalar(cdim) && cdim == fix(cdim) && cdim >= 1)
                error('glinregr1:invalidarg', ...
                    'cdim should be a positive integer scalar.');
            end
            
            obj = obj@gcondmodel(cdim);
        end
        
        function [U1, U2] = conj_update(obj, data, W)
            % Compute the conjugate update for Gaussian prior on c
            %
            %   [U1, U2] = conj_update(obj, {X, y, s2});
            %   [U1, U2] = conj_update(obj, {X, y, s2}, W);
            %
            %       computes the conjuate update quantities for
            %       Gaussian prior on the coefficients. 
            %
            %       The second argument gives the observed x-y pairs.
            %       It is a cell array comprised of three items:
            %       - X:  the matrix input vectors, size is d x n
            %             where d is the input dimension, and n
            %             is the number of observations.
            %       - y:  the row vector of output values, size is 1 x n.
            %       - s2: the variance of outputs, which can be either
            %             a scalar or a 1 x n row vector.
            %
            %       W is the matrix of observation weights. Its size
            %       is m x n, where m is the number of models to estimate.
            %       It can be omitted or empty, in which case, all 
            %       samples are with weight 1.
            %
            
            % verify input
            
            [X, y, s2, n] = check_data(obj, data);
            
            if ~isempty(W)
                if ~(isfloat(W) && isreal(W) && ndims(W) == 2 && size(W, 2) == n)
                    error('glinregr1:conj_update:invalidarg', ...
                        'W should be either empty or an m x n real matrix.');
                end                
            end
            
            % main
            
            % combine weights and variances
            if isempty(W)
                W = 1 ./ s2;
                m = 1;
            else
                m = size(W, 1);
                if isscalar(s2)
                    W = W * (1 / s2);
                else
                    if m == 1
                        W = W ./ s2;
                    else
                        W = bsxfun(@times, W, 1 ./ s2);
                    end
                end
            end
            
            % compute conjugate updates
            
            if m == 1
                Wy = W .* y;
            else
                Wy = bsxfun(@times, W, y);
            end
            
            U1 = X * Wy';
            
            d = size(X, 1);            
            if m == 1
                if isscalar(W)
                    U2 = (X * W) * X';                    
                else
                    U2 = bsxfun(@times, X, W) * X';
                end
                U2 = 0.5 * (U2 + U2');
            else                
                U2 = zeros(d, d, m, class(X(1) * W(1)));
                for i = 1 : m
                    cS = bsxfun(@times, X, W(i,:)) * X';
                    U2(:,:,i) = 0.5 * (cS + cS)';
                end
            end
            
            U2 = pdmat_gen(U2);
        end
        
        
        function L = eval_loglik(obj, C, data)
            % Evaluates the log likelihood of data
            %
            %   L = eval_loglik(obj, C, {X, y, s2});
            %       evaluates the log-likelihood of data with
            %       respect to the coefficients given by columns of
            %       C. 
            %
            %       The data is given in form of {X, y, s2}, which
            %       is comprised of the input vectors, output values,
            %       and output variances.
            %
            %       If the input dimension is d, then the size of C
            %       should be d x m, and then L is of size m x n,
            %       where n is the number of samples in data.
            %
            
            [X, y, s2] = check_data(obj, data);
            d = size(X, 1);
            
            if ~(isfloat(C) && isreal(C) && ndims(C) == 2 && size(C,1) == d)
                error('glinregr1:invalidarg', ...
                    'C should be a real matrix of size d x m.');
            end
            m = size(C, 2);
            
            if m == 1
                Z = y - C' * X;
                D = (Z.^2) ./ s2;
            else
                Z = bsxfun(@minus, y, C' * X);
                D = bsxfun(@times, Z.^2, 1 ./ s2);
            end
            
            a = log(2*pi) + log(s2);
            
            if isscalar(a) || m == 1
                L = D + a;
            else
                L = bsxfun(@plus, D, a);
            end
            
            L = (-0.5) * L;            
        end
                        
    end
    
    
    methods(Access='private')
        
        function [X, y, s2, n] = check_data(obj, data)
            % verify the integrity of the data
            
            d = obj.cdim;
            
            if ~(iscell(data) && numel(data) == 3)
                error('glinregr1:invaliddata', ...
                    'data should be a cell array with three cells.');
            end
            
            X = data{1};
            y = data{2};
            s2 = data{3};
            
            if ~(isfloat(X) && isreal(X) && ndims(X) == 2 && size(X,1) == d)
                error('glinregr1:invaliddata', ...
                    'X should be a real matrix with size(X,1) == obj.cdim.');
            end
            
            if ~(isfloat(y) && isreal(y) && ndims(y) == 2 && size(y,1) == 1)
                error('glinregr1:invaliddata', ...
                    'y should be a real row vector.');
            end
            
            if ~(isfloat(s2) && isreal(s2) && ndims(s2) == 2)
                error('glinregr1:invaliddata', ...
                    's2 should be a real row vector.');
            end
            
            n = size(X, 2);
            
            if size(y, 2) ~= n
                error('glinregr1:invaliddata', ...
                    'X and y should have the same number of columns.');
            end
            
            if ~(isscalar(s2) || (size(s2, 1) == 1 && size(s2, 2) == n))
                error('glinregr1:invaliddata', ...
                    's2 should be either a scalar or a 1 x n row vector.');
            end            
        end
        
        
    end
    
end