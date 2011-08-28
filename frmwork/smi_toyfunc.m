classdef smi_toyfunc < smi_func
    % This is a toy sampling function class
    %
    % The purposes of this class:
    %   - for debugging and/or testing a smi framework and gibbs sampling
    %   - serve as an example to illustrate how to write a sample step
    %    class
    %
    
    % Created by Dahua Lin, Aug 24, 2011
    %
        
    %% methods
    
    methods
        
        function obj = smi_toyfunc(ins, outs)
            % Construct a smi_toyfunc object
            %
            %   obj = smi_toyfunc(ins, outs);
            %       constructs the object by specifying the input
            %       and output slots.
            %
            %       Both ins and outs are a struct, whose fields
            %       are names of the slots, and the values of the
            %       corresponding fields are the dimensions
            %       of vector that can be used for the slot.
            %
            
            if ~(isstruct(ins) && isscalar(ins))
                error('smi_toyfunc:invalidarg', ...
                    'ins should be a struct scalar.');
            end
            
            if ~(isstruct(outs) && isscalar(outs))
                error('smi_toyfunc:invalidarg', ...
                    'outs should be a struct scalar.');
            end
            
            in_names = fieldnames(ins);
            out_names = fieldnames(outs);
            
            n_in = numel(in_names);
            n_out = numel(out_names);
            
            inlets = struct('name', in_names, 'type', 'double', 'size', []);
            outlets = struct('name', out_names, 'type', 'double', 'size', []);
            
            for i = 1 : n_in
                inlets(i).size = ins.(in_names{i});
            end
            
            for i = 1 : n_out
                outlets(i).size = outs.(out_names{i});
            end    
            
            obj = obj@smi_func(inlets, outlets);
        end
        
    
        function tf = test_slots(obj, inflags, outflags)
            % Test whether a particular connection config is valid.
            %
            %   tf = obj.test_slots(inflags, outflags);
            %
            
            n_in = numel(obj.input_slots);
            n_out = numel(obj.output_slots);
            
            if ~(islogical(inflags) && isequal(size(inflags), [1, n_in]))
                error('smi_toyfunc:invalidarg', ...
                    'inflags should be a logical vector of size 1 x n_in.');
            end
            
            if ~(islogical(outflags) && isequal(size(outflags), [1, n_out]))
                error('smi_toyfunc:invalidarg', ...
                    'outflags should be a logical vector of size 1 x n_out.');
            end
            
            tf = true; 
        end
        
                
        function varargout = evaluate(obj, outflags, varargin)
            % Run the sampling step
            %
            %   [out1, out2, ...] = obj.run(outflags, in1, in2, ...);
            %
            %       Runs the sampling step to produce required outputs.
            %
            %       outflags specifies which output is wanted. 
            %
            %       If outflags(i) is true, it means that the caller 
            %       wants the output via the i-th output slot. 
            %
            %       If outflags(i) is false, it means the caller does
            %       not want the corresponding output, the function
            %       may still yield a valid output for outi, or simply
            %       set it to empty.
            %
            
            n_in = numel(obj.input_slots);
            n_out = numel(obj.output_slots);
            varargout = cell(1, n_out);
            
            for i = 1 : n_in
                d = obj.input_slots(i).size;
                cv = varargin{i};
                if ~isempty(cv) && ~isequal(size(cv), [d, 1])
                    error('smi_toyfunc:invalidarg', 'Invalid input argument');
                end
            end            
            
            for i = 1 : n_out                  
                if outflags(i) 
                    d = obj.output_slots(i).size;
                    varargout{i} = rand(d, 1);                    
                end                
            end            
        end
        
    end
    
end


