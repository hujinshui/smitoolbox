classdef xgs_func
    % A base class for all xgs functions
    %
    
    % Created by Dahua Lin, on Aug 24, 2011
    %
    
    properties(Abstract)
        num_input_slots;    % the number of input slots
        num_output_slots;   % the number of output slots
    end
    
    
    methods(Abstract)
        name = get_slot_name(obj, dir, i);
        % Gets the name of a slot
        %
        %   obj.get_slot_name('in', i);
        %       returns the name of the i-th input slot
        %
        %   obj.get_slot_name('out', i);
        %       returns the name of the i-th output slot
        %
        
        info = get_slot_info(obj, name);
        % Retrieves the information of a slot given its name
        % 
        %   info = obj.get_slot_info(name);
        %
        %       The returned info is a struct comprised of at least 
        %       the following fields:
        %
        %       - in_id:    the index of the named slot in input slots
        %                   ([] if the slot is not for input)
        %       - out_id:   the index of the named slot in output slots
        %                   ([] if the slot is not for output)
        %       - type:     the type(class) of variable to be transfered
        %                   via this slot
        %       - size:     the size of the variable to be transfered via
        %                   this slot, which should be a single number for
        %                   a column-vector, a pair of numbers for a matrix,
        %                   etc. 
        
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
