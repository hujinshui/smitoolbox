classdef psampler
    % The base class for all samplers that draw samples from a specific
    % distribution
    %
    
    % Created by Dahua Lin, on Mar 19, 2010
    %
    
    properties(GetAccess='public', SetAccess='protected')
        
        dim;      % the dimensionality of the underlying space  
        nmodels;  % the number of samplers contained in the object
        rstream;  % the associated random number stream
    end
    
    methods(Abstract)
        
        X = draw(obj, i, n);
        % Draw independent identically distributed samples from the
        % underlying distribution
        %
        %   i: the index or index array of selected model(s)
        %      if i is empty, then all models are selected.
        %
        %   n: if n is a scalar, then draw n samples from each selected
        %      model consequtively, otherwise, n can be a 1 x m vector
        %      (m is the number of selected models), and it draws
        %      n(j) samples from the j-th selected model.
        %      
    end
    
    methods(Access='protected')
        
        function [i, n] = check_draw_args(obj, i, n)
            % a convenient function to check the validity of the arguments
            % input to the draw method, and convert them to full form
            %
            
            m = obj.nmodels;
            
            if isempty(i)
                i = 1 : m;
            else                
                assert(isnumeric(i) && isvector(i) && all(i == fix(i) & i >= 1 & i <= m), ...
                    'psampler:draw:invalidarg', ...
                    'The model index(indices) are specified incorrectly.');                
            end
            
            if size(i, 1) > 1
                i = i.';
            end
            m = numel(i);
            
            assert((isscalar(n) || isequal(size(n), [1 m])) && all(n == fix(n) & n >= 0), ...
                'psampler:draw:invalidarg', ...
                'The number(s) of samples are specified incorrectly.');
            
            if isscalar(n) && m > 1
                n = n * ones(1, m);
            end                        
        end
        
    end
    
end