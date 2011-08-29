classdef rand_labeler < smi_func
    % An SMI function to randomly assign labels 
    %
    %   Suppose we intend to assign N samples, each with a label
    %   in [1, K]
    %   
    %   There are three ways to generate random labeling.
    %
    %   - generate a vector of labels, of size 1 x N (opchar = 'v');
    %   - generate an K x N matrix, such that each column is
    %     an indicator vector (opchar = 'b');
    %   - generate an K x N matrix, such that each column is a 
    %     distribution of labels (opchar = 'p');
    %
    
    % Created by Dahua Lin, on Aug 29, 2011
    %
    
    properties(GetAccess='public', SetAccess='private')
        K;          % the number of distinct labels 
        N;          % the number of samples to assign label to
        opchar;     % the char to indicate how to generate the labels
    end
    
    methods
        
        %% Construction
        
        function obj = rand_labeler(K, N, op)
            % Constructs a random labeler
            %
            %   obj = rand_labeler(K, N, op);
            %
            %   Input arguments:
            %       - K:    the number of distinct labels
            %       - N:    the number of samples to assign labels to
            %       - op:   the char indicating the way it generates,
            %               which can be either of 'v', 'b', or 'p'. 
            %
            
            % verify input
            
            if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 1)
                error('rand_labeler:invalidarg', ...
                    'K should be a positive integer scalar.');
            end
            
            if ~(isnumeric(N) && isscalar(N) && N == fix(N) && N >= 0)
                error('rand_labeler:invalidarg', ...
                    'N should be a non-negative integer scalar.');
            end
            
            if ~(ischar(op) && isscalar(op))
                error('rand_labeler:invalidarg', ...
                    'op should be a char scalar.');
            end
                        
            outlet.name = 'result';
            outlet.type = 'double';
            
            op = lower(op);                        
            switch op
                case 'v'
                    outlet.size = [1, N];
                case {'b', 'p'}
                    outlet.size = [K, N];
                otherwise
                    error('rand_labeler:invalidarg', ...
                        'The op char is invalid.');
            end
            
            
            obj = obj@smi_func([], outlet);
            
            obj.K = K;
            obj.N = N;
            obj.opchar = op;
            
        end
        
        
        %% Implementation
        
        function tf = test_slots(obj, ~, ~) %#ok<MANU>
            tf = true;            
        end
        
        
        function L = evaluate(obj, ~)           
            % Performs random label generation
            
            K_ = obj.K;
            N_ = obj.N;
            
            switch obj.opchar
                case 'v'
                    L = randi(K_, [1, N_]);
                case 'b'
                    I = randi(K_, [1, N_]);
                    L = zeros(K_, N_);
                    L(I + (0:N_-1) * K_) = 1;
                case 'p'
                    L = rand(K_, N_);
                    L = bsxfun(@times, L, 1 ./ sum(L, 1));
            end            
        end
        
                
    end
    
end
