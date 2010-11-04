classdef binlabeler_std < labeler
    % The standard binary labeler
    %
    
    % Created by Dahua Lin, on Apr 23, 2010
    %
    
    properties(GetAccess='public', SetAccess='private')
        prior;  % the beta distribution prior
    end
    
    methods
        
        function obj = binlabeler_std(pri)
            % Constructs a standard binary labeler
            %
            %   obj = binlabeler_std();
            %       constructs a standard binary labeler without
            %       prior, meaning that the estimation of labeling
            %       parameter would be in a maximum likelihood manner.
            %
            %   obj = binlabeler_std(pri);
            %       constructs a standard binary labeler with
            %       a beta prior given by pri.
            %
            %       pri is an object of class betad. 
            %       The estimation of labeling parameter is in 
            %       maximum a posteriori manner.
            %
            
            if nargin >= 1 && ~isempty(pri)
                if ~(isa(pri, 'betad') && pri.nmodels == 1)
                    error('binlabeler_std:invalidarg', ...
                        'pri should be an object of class betad with single model.');
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
            
            if ~(isfloat(Q) && isreal(Q) && ndims(Q) == 2 && size(Q, 1) == 1)
                error('binlabeler_std:invalidarg', ...
                    'Q should be a real row vector.');
            end
            
            s = sum(Q, 2);
            n = size(Q, 2);
            
            if isempty(pri)
                p = s / n;
            else
                a = pri.alpha;
                b = pri.beta;
                p = (s + (a - 1)) / (n + (a+b-2));
            end
        end
        
        function Pb = accept_param(obj, p) %#ok<MANU>
            % Accept a labeling prior p and form a labeling problem
            %
            %   Pb = accept_param(obj, p);
            %
            
            if ~(isfloat(p) && isreal(p) && isscalar(p))
                error('binlabeler_std:invalidarg', ...
                    'p should be a real scalar.');
            end
            
            Pb = binlabelpb_std(p);
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