classdef multlabelpb_std < multlabelpb
    % The class to represent the standard multi-class labeling problem
    %
    %   Each label is considered as independent, and a multi-class 
    %   discrete distribution is used as the prior of labeling.
    %   
    
    properties(GetAccess='public', SetAccess='private')
        pri;        % the prior probability of labels being 1 [scalar or 1 x n]

        nclasses;   % the number of classes (empty: arbitrary)    
        nsamples;   % the number of samples (empty: arbitrary)
    end
    
    
    methods
        function obj = multlabelpb_std(p0)
            % Constructs a standard multi-class labeling problem
            %
            %   obj = multlabelpb_std(p0);
            %       constructs the standard multi-class labeling problem 
            %       with prior distribution of the labels given by p0.
            %
            %       p0 can be given in either of the following forms:
            %       - empty or omitted, which means that each class has
            %         equal prior probability.
            %       - a m x 1 column vector that gives the common prior 
            %         probabilities. In particular, p0(k) is the prior
            %         probability of being generated from the k-th model.            
            %       - a m x n matrix, with p0(:,j) being the prior
            %         probabilities for the j-th sample.
            %
            %       It is important to ensure that each column of p0
            %       sums to 1.
            %
            
            if nargin >= 1 && ~isempty(p0)
                if isfloat(p0) && ndims(p0) == 2                    
                    [m, n] = size(p0);
                    obj.pri = p0;
                    obj.nclasses = m;
                    if n > 1
                        obj.nsamples = n;
                    end
                else
                    error('multilabeler_std:invalidarg', ...
                        'p0 should be a real matrix.');                
                end                            
            end
        end
        
        
        function Q = infer_q(obj, L)
            % Infer the probabilities of labels given the likelihoods
            %
            %   Q = obj.infer_q(L);
            %
                        
            % verify input
            
            m = obj.nclasses;
            n = obj.nsamples;
            
            if~(isfloat(L) && isreal(L) && ndims(L) == 2)
                error('multlabelpb_std:infer_q:invalidarg', ...
                    'L should be a real matrix');
            end
            
            if ~isempty(m)
                if size(L, 1) ~= m
                    error('multlabelpb_std:invalidarg', ...
                        'L should contain m rows, where m == obj.nclasses');
                end
            end
            
            if ~isempty(n)
                if size(L, 2) ~= n
                    error('multlabelpb_std:invalidarg', ...
                        'L should contain n columns, where n == obj.nsamples');
                end
            end
                        
            % compute
            
            p0 = obj.pri;
            if ~isempty(p0)
                lpri = log(p0);
                if isequal(size(lpri), size(L))
                    H = lpri + L;
                else                    
                    H = bsxfun(@plus, lpri, L);
                end
            else
                H = L;
            end
            
            H = bsxfun(@minus, H, max(H, [], 1));
            Q = exp(H);
            Q = bsxfun(@times, Q, 1 ./ sum(Q, 1));
        end
        
        
        function v = eval_qscore(obj, Q)
            % Evaluates the score of the labeling given by q
            %
            %   v = obj.eval_qscore(Q);
            %
                        
            % verify input
            
            m = obj.nclasses;
            n = obj.nsamples;
            
            if~(isfloat(Q) && isreal(Q) && ndims(Q) == 2)
                error('multlabelpb_std:eval_qscore:invalidarg', ...
                    'Q should be a real matrix');
            end
            
            if ~isempty(m)
                if size(Q, 1) ~= m
                    error('multlabelpb_std:eval_qscore:invalidarg', ...
                        'Q should contain m rows, where m == obj.nclasses');
                end
            end
            
            if ~isempty(n)
                if size(Q, 2) ~= n
                    error('multlabelpb_std:eval_qscore:invalidarg', ...
                        'Q should contain n columns, where n == obj.nsamples');
                end
            end
            
            % compute
            
            p0 = obj.pri;
            if isempty(p0)
                [m, n] = size(Q);                
                v = -log(m) * n;
            else
                if size(p0, 2) == 1
                    v = sum(log(p0)' * Q);
                else
                    v = sum(sum(log(p0) .* Q, 1));
                end
            end
            
            ent = sum(ddentropy(Q));
            v = v + sum(ent);
        end
        
    end
end

