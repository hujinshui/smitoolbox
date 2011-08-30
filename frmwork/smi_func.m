classdef smi_func
    % A base class for all smi framework functional module
    %
    
    % Created by Dahua Lin, on Aug 24, 2011
    %
    
    properties(GetAccess='public', SetAccess='private')
        input_slots;    % the list of input slots
        output_slots;   % the list of output slots
    end
        
    methods
        function obj = smi_func(in_slots, out_slots)
            % Constructs an smi_func base object
            %
            %   obj = smi_func(in_slots, out_slots);
            %       constructs an smi_func base object, given 
            %       the specification of its input slots and
            %       output slots.
            %
            %       Both in_slots and out_slots are struct array,
            %       with each element containing basic information
            %       of a slot.
            %
            %       Specifically, each element has the following fields:
            %       
            %       - name:     the name of the slot
            %       - type:     the type(class) of variable associated
            %                   with this slot            
            %       - size:     the size of the variable associated with
            %                   this slot, which should be a single number 
            %                   for a column-vector, a pair of numbers for 
            %                   a matrix, etc.
            %
            
            if ~isempty(in_slots)
                if ~(isstruct(in_slots) && ...
                        all(isfield(in_slots, {'name', 'type', 'size'})))
                    error('smi_func:invalidarg', ...
                        'The input slot specification is invalid.');
                end
            else
                in_slots = [];
            end
            
            if ~isempty(out_slots)
                if ~(isstruct(out_slots) && ...
                        all(isfield(out_slots, {'name', 'type', 'size'})))
                    error('smi_func:invalidarg', ...
                        'The output slot specification is invalid.');
                end
            else
                out_slots = [];
            end
            
            obj.input_slots = in_slots;
            obj.output_slots = out_slots;
        end
        
    end
        
    
    methods(Abstract)
        
        tf = test_slots(obj, inflags, outflags);
        % Test whether a given input/output pattern is acceptable
        %
        %   tf = obj.test_slots(inflags, outflags);
        %
        %       Suppose, the number of inputs and that of outputs are
        %       respectively n and m. Then inflags and outflags should be
        %       respectively 1 x n and 1 x m logical vectors, with each
        %       entry indicating whether the corresponding slot is used.
        %
                
        varargout = evaluate(obj, outflags, varargin);
        % Evaluate the function
        %
        %    ... = obj.evaluate(outflags, ...);
        %
        %   Here, outflags is an 1 x num_output_slots logical vector.
        %   outflags(i) indicates whether the i-th output is wanted.
        %
    end
    
end
