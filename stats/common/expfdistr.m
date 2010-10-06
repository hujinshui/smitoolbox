classdef expfdistr
    % The base class of all exponential family distribution
    %
    
    % Created by Dahua Lin, on Mar 18, 2010
    % Modified by Dahua Lin, on Apr 13, 2010
    %
    
    properties(GetAccess='public', SetAccess='protected')
        nmodels;    % the number of models contained in the object        
        logpar;     % the value of log-partition function
    end
    
    methods        
        
        function P = prob(obj, X, i)
            % Evaluate the pmf or pdf of the given samples
            %
            %   P = obj.prob(X);
            %   P = obj.prob(X, i);
            %       evaluates the probability mass function or probability
            %       density function of the selected models on the samples
            %       given by X.
            %
            %       i is the index or a vector of indices of the selected 
            %       model(s). If i is omitted, then all models will be 
            %       evaluated.
            %
            %       If there are m selected models, and n samples, then 
            %       in the output, P should be a matrix of size m x n, 
            %       with P(i, j) being the result of evaluating i-th model
            %       on the j-th sample
            %
            
            if nargin <= 2
                LP = logprob(obj, X); 
            else
                LP = logprob(obj, X, i);
            end
            P = exp(LP);
        end        
        
        function LP = logprob(obj, X, i)
            % Evaluate the log pmf or pdf of the given samples            
            %
            %   LP = obj.logprob(X);
            %       evaluates the logarithm of the probability mass
            %       function or probability density function of the
            %       selected models on the samples given by X.
            %
            %       i is the index or a vector of indices of the selected 
            %       model(s). If i is omitted, then all models will be 
            %       evaluated.
            %
            %       If there are m selected models and n samples, then 
            %       in the output LP should be a matrix of size m x n, 
            %       with LP(i, j) being the result of evaluating i-th 
            %       model on the j-th sample            
            %
            %       If there are m selected models and n samples, then 
            %       in the output LP should be a matrix of size m x n, 
            %       with LP(i, j) being the result of evaluating i-th 
            %       model on the j-th sample
            %
            
            lpv = obj.logpar;
            if isempty(lpv)   
                error('expfdistr:logprob:nologpar', ...
                    'The log-partition function is not available.');
            end
            
            if nargin <= 2
                i = [];                
            else
                lpv = lpv(i(:));
            end
            
            Y = compute_clinterm(obj, X, i);            
            LB = compute_logbase(obj, X);
            H = LB - lpv;
            
            if isscalar(H) || size(Y, 1) == 1
                LP = Y + H;
            else
                LP = bsxfun(@plus, Y, H.');
            end                            
        end                                   
    end        
    
    methods(Abstract)
        
        LT = compute_clinterm(obj, X, i);            
        % Evaluate the canonical linear term of the selected models on 
        % given samples                        
           
        L = compute_logbase(obj, X);
        % compute the logarithm of base measure on given samples
        % it can return 0 or a scalar the base measure is uniform at all
        % samples
        
    end
        
end
