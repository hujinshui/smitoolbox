classdef binlabelpb_std < binlabelpb
    % The class to represent the standard binary labeling problem
    %
    %   Each label is considered as independent, and a bernoulli
    %   distribution is used as the prior of labeling.
    %   
    
    properties(GetAccess='public', SetAccess='private')
        pri;         % the prior probability of labels being 1 [scalar or 1 x n]
        logit = 0;   % the logit = log( pri / (1 - pri) ).
                
        nsamples;   % the number of samples (empty: arbitrary)
    end
    
    
    methods
        function obj = binlabelpb_std(p0)
            % Constructs a standard binary labeling problem
            %
            %   obj = binlabeler_std(p0);
            %       constructs the standard binary labeling problem with 
            %       prior distribution of being labeled 1 given by p0.
            %
            %       p0 can be given in either of the following forms:
            %       - empty or omitted, which means that each class has
            %         equal prior probability.
            %       - a scalar that gives the common prior probability
            %         of being generated from model-1.
            %       - a 1 x n row vector with p0(i) giving the prior
            %         probability of i-th sample being generated from
            %         model-1. In this case, the labeler can only be
            %         used to label a fixed number of samples.            
            %
            
            if nargin >= 1 && ~isempty(p0)
                if isfloat(p0)
                    if isscalar(p0) && p0 >= 0 && p0 <= 1
                        obj.pri = p0;
                    elseif ndims(p0) == 2 && size(p0, 1) == 1
                        obj.pri = p0;
                        obj.nsamples = size(p0, 2);
                    else
                        error('binlabelpb_std:invalidarg', ...
                            'p0 should be either a scalar or a row vector.');
                    end
                    
                    obj.logit = log(p0 ./ (1 - p0));
                else
                    error('binlabelpb_std:invalidarg', ...
                        'p0 should be of class double or single.');                
                end
            else
                obj.pri = 0.5;
            end
        end
        
        
        function q = infer_q(obj, L1, L0)
            % Infer the probabilities of labels given the likelihoods
            %
            %   q = obj.infer_q(L1, L0);
            %
                        
            % verify input
            
            n = obj.nsamples;
            if isempty(n)
                n = size(L1, 2);
            end
            if ~(isfloat(L1) && isreal(L1) && isequal(size(L1), [1 n]))
                error('binlabelpb_std:infer_q:invalidarg', ...
                    'L1 should be an 1 x n real row vector.');
            end
            if ~(isfloat(L0) && isreal(L0) && isequal(size(L0), [1 n]))
                error('binlabelpb_std:infer_q:invalidarg', ...
                    'L0 should be an 1 x n real row vector.');
            end
            
            % compute
            
            lg = obj.logit;
            if ~isequal(lg, 0)
                r = (L0 - L1) - lg;
            else
                r = L0 - L1;
            end
            
            q = 1 ./ (1 + exp(r));            
        end
        
        function v = eval_qscore(obj, q)
            % Evaluates the score of the labeling given by q
            %
            %   v = obj.eval_qscore(q);
            %
            
            % verify input
            
            n = obj.nsamples;
            if isempty(n)
                n = size(q, 2);
            end
            if ~(isfloat(q) && isreal(q) && isequal(size(q), [1 n]))
                error('binlabelpb_std:eval_qscore:invalidarg', ...
                    'q should be an 1 x n real row vector.');
            end
                        
            % compute
            
            p0 = obj.pri;
            if isequal(p0, 0.5)
                v = -log(2) * n;
            else
                v = sum(q .* log(p0) + (1 - q) .* log(1 - p0));
            end         
            
            ent = sum(ddentropy([q; 1-q]));
            v = v + sum(ent);
        end
        
    end
end

