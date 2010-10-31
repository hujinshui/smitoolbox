classdef multlabeler_std < labeler
    % The standard multi-class labeler
    %
    
    % Created by Dahua Lin, on Apr 23, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        prior;  % the dirichlet distribution prior
    end
    
    methods
        function obj = multlabeler_std(pri)
            % Constructs a standard multi-class labeler
            %
            %   obj = multlabeler_std();
            %       constructs a standard multi-class labeler without
            %       prior, meaning that the estimation of labeling
            %       parameter would be in a maximum likelihood manner.
            %
            %   obj = multlabeler_std(pri);
            %       constructs a standard multi-class labeler with
            %       a dirichlet prior given by pri.
            %
            %       pri is an object of class dirichletd. 
            %       The estimation of labeling parameter is in 
            %       maximum a posteriori manner.
            %
            
            if nargin >= 1 && ~isempty(pri)
                if ~(isa(pri, 'dirichletd') && pri.nmodels == 1)
                    error('multlabeler_std:invalidarg', ...
                        'pri should be an object of class dirichletd with single model.');
                end
                
                obj.prior = pri;
            end            
        end        
        
        function p = estimate_param(obj, Q)
            % Estimates the labeling parameter from soft assignment
            %
            %   p = estimate_param(obj, Q);
            %
            
            pri = obj.prior;                                    
            
            if ~(isfloat(Q) && isreal(Q) && ndims(Q) == 2)
                error('multlabeler_std:invalidarg', ...
                    'Q should be a real matrix.');
            end
            
            if ~isempty(pri)
                if size(Q, 1) ~= pri.dim
                    error('multlabeler_std:invalidarg', ...
                        'The size of Q is not consistent with the prior.');
                end
            end
            
            p = sum(Q, 2);
                                    
            if ~isempty(pri)
                p = p + (pri.alpha - 1);
            end
            
            p = p / sum(p);
        end
        
        function Pb = accept_param(obj, p)
            % Accept a labeling prior p and form a labeling problem
            %
            %   Pb = accept_param(obj, p);
            %
            
            pri = obj.prior;
            
            if (~isfloat(p) && isreal(p) && ndims(p) == 2 && size(p, 2) == 1)
                error('multlabeler_std:accept_param:invalidarg', ...
                    'p should be a real column vector.');
            end
            
            if ~isempty(pri)
                if size(p, 1) ~= pri.dim
                    error('multlabeler_std:accept_param:invalidarg', ...
                        'The size of p is not consistent with the prior.');
                end
            end
            
            Pb = multlabelpb_std(p);
        end
        
        function L = eval_logpri(obj, p)
            % Estimates the log-prior of the parameter
            %
            %   L = eval_logpri(obj, p);
            %        
            
            pri = obj.prior;
            if isempty(pri)                
                L = 0;
            else
                L = pri.logprob(p);
            end
        end
        
    end
end