classdef toy_spstep
    % This is a toy sampling step class
    %
    % The purposes of this class:
    %   - for debugging and/or testing a xgs model framework
    %   - serve as an example to illustrate how to write a sample step
    %    class
    %
    
    % Created by Dahua Lin, Aug 24, 2011
    %
    
    %% properties
    
    properties(GetAccess='public', SetAccess='private')
        inlet_names;    % the cell array of names of input slots
        outlet_names;   % the cell array of names of output slots
        
        in_dims;        % the vector dimensions for input slots
        out_dims;       % the vector dimensions for output slots
    end
    
    properties(Dependent)
        num_input_slots;
        num_output_slots;        
    end    
    
    methods        
        function v = get.num_input_slots(obj)
            v = numel(obj.inlet_names);
        end        
        function v = get.num_output_slots(obj)
            v = numel(obj.outlet_names);
        end        
    end
    
    
    %% methods
    
    methods
        
        function obj = toy_spstep(ins, outs)
            % Construct a toy_spstep object
            %
            %   obj = toy_spstep(ins, outs);
            %       constructs the object by specifying the input
            %       and output slots.
            %
            %       Both ins and outs are a struct, whose fields
            %       are names of the slots, and the values of the
            %       corresponding fields are the dimensions
            %       of vector that can be used for the slot.
            %
            
            if ~(isstruct(ins) && isscalar(ins))
                error('toy_spstep:invalidarg', ...
                    'ins should be a struct scalar.');
            end
            
            if ~(isstruct(outs) && isscalar(outs))
                error('toy_spstep:invalidarg', ...
                    'outs should be a struct scalar.');
            end
            
            obj.inlet_names = fieldnames(ins);
            obj.outlet_names = fieldnames(outs);
            
            ps = struct2cell(ins);
            obj.in_dims = [ps{:}];
            
            qs = struct2cell(outs);
            obj.out_dims = [qs{:}];
        end
        
        
        function info = get_slot_info(obj, name)
            % Get the information of a slot
            %
            %   info = obj.get_slot_info(name);
            %       gets the information of a slot with the input name.
            %
            %       info is a struct with the following fields:
            %       - in_id:    The index of the slot in the inputs.
            %                   ([] if the named slot is not for input)
            %       - out_id:   The index of the slot in the outputs.
            %                   ([] if the named slot is not for output)
            %
            %       The function raises an error if the given name
            %       is not the same of any slot.
            %
            
            if ~ischar(name)
                error('toy_spstep:invalidarg', ...
                    'The slot name must be a char string.');
            end
            
            [tf_i, id_i] = ismember(name, obj.inlet_names);
            [tf_o, id_o] = ismember(name, obj.outlet_names);
            
            if ~tf_i && ~tf_o
                error('toy_spstep:invalidname', ...
                    'The input name %s is not a slot name.', name);
            end
            
            if ~tf_i; id_i = []; end
            if ~tf_o; id_o = []; end
            
            info.in_id = id_i;
            info.out_id = id_o;
        end
        
        
        function tf = test_slots(obj, inflags, outflags)
            % Test whether a particular connection config is valid.
            %
            %   tf = obj.test_slots(inflags, outflags);
            %
            
            n_in = obj.num_input_slots;
            n_out = obj.num_output_slots;
            
            if ~(islogical(inflags) && isequal(size(inflags), [1, n_in]))
                error('toy_spstep:invalidarg', ...
                    'outflags should be a logical vector of size 1 x n_out.');
            end
            
            if ~(islogical(outflags) && isequal(size(outflags), [1, n_out]))
                error('toy_spstep:invalidarg', ...
                    'outflags should be a logical vector of size 1 x n_out.');
            end
            
            tf = true; 
        end
        
        
        function tf = verify_args(obj, outflags, varargin)
            % Verify whether the input arguments are valid
            %
            %   tf = obj.verify_args(outflags, ...);
            %       If the input arguments are valid and sufficient
            %       for the desired output patterns, it returns true;
            %       otherwise, false.
            %
            
            n_out = obj.num_output_slots;
            if ~(islogical(outflags) && isequal(size(outflags), [1, n_out]))
                error('toy_spstep:invalidarg', ...
                    'outflags should be a logical vector of size 1 x n_out.');
            end
            
            n_in = obj.num_input_slots;            
            idims = obj.in_dims;
            
            if numel(varargin) ~= n_in
                error('toy_spstep:invalidarg', ...
                    'The number of inputs are not valid.');
            end
            
            tf = true;
            for i = 1 : n_in
                v = varargin{i};
                tf = isempty(v) || ...
                    (isfloat(v) && ndims(v) == 2 && ...
                    size(v,1) == idims(i) && size(v,2) == 1);
                if ~tf
                    return;
                end
            end            
        end
        
        
        function varargout = run(obj, outflags, varargin)
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
            
            n_out = obj.num_output_slots;
            varargout = cell(1, n_out);
            
            odims = obj.out_dims;
            for i = 1 : n                
                if outflags(i) 
                    varargout{i} = rand(odims(i), 1);                    
                end                
            end            
        end
        
    end
    
end


