classdef toy_spfunc < xgs_func
    % This is a toy sampling function class
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
        
        function obj = toy_spfunc(ins, outs)
            % Construct a toy_spfunc object
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
                error('toy_spfunc:invalidarg', ...
                    'ins should be a struct scalar.');
            end
            
            if ~(isstruct(outs) && isscalar(outs))
                error('toy_spfunc:invalidarg', ...
                    'outs should be a struct scalar.');
            end
            
            obj.inlet_names = fieldnames(ins);
            obj.outlet_names = fieldnames(outs);
            
            ps = struct2cell(ins);
            obj.in_dims = [ps{:}];
            
            qs = struct2cell(outs);
            obj.out_dims = [qs{:}];
        end
        
        
        function name = get_slot_name(obj, dir, i)
            % Gets the name of a slot
            %
            %   obj.get_slot_name('in', i);
            %       returns the name of the i-th input slot
            %
            %   obj.get_slot_name('out', i);
            %       returns the name of the i-th output slot
            %
            
            if strcmpi(dir, 'in')
                name = obj.inlet_names{i};
            else
                name = obj.outlet_names{i};
            end                
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
                error('toy_spfunc:invalidarg', ...
                    'The slot name must be a char string.');
            end
                                                
            [tf_i, id_i] = ismember(name, obj.inlet_names);
            
            if tf_i
                info.in_id = id_i;
                info.out_id = [];
                info.type = 'double';
                info.size = obj.in_dims(id_i);
                return;
            end
            
            [tf_o, id_o] = ismember(name, obj.outlet_names);
            
            if tf_o
                info.in_id = [];
                info.out_id = id_o;
                info.type = 'double';
                info.size = obj.out_dims(id_o);                
            else
                error('toy_spfunc:invalidname', ...
                    'The input name %s is not a slot name.', name);
            end
        end
        
        
        function tf = test_slots(obj, inflags, outflags)
            % Test whether a particular connection config is valid.
            %
            %   tf = obj.test_slots(inflags, outflags);
            %
            
            n_in = obj.num_input_slots;
            n_out = obj.num_output_slots;
            
            if ~(islogical(inflags) && isequal(size(inflags), [1, n_in]))
                error('toy_spfunc:invalidarg', ...
                    'outflags should be a logical vector of size 1 x n_out.');
            end
            
            if ~(islogical(outflags) && isequal(size(outflags), [1, n_out]))
                error('toy_spfunc:invalidarg', ...
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
            
            n_in = obj.num_input_slots;
            n_out = obj.num_output_slots;
            varargout = cell(1, n_out);
            
            idims = obj.in_dims;
            for i = 1 : n_in
                cv = varargin{i};
                if ~isempty(cv) && ~isequal(size(cv), [idims(i), 1])
                    error('toy_spfunc:invalidarg', 'Invalid input argument');
                end
            end            
            
            odims = obj.out_dims;
            for i = 1 : n_out                
                if outflags(i) 
                    varargout{i} = rand(odims(i), 1);                    
                end                
            end            
        end
        
    end
    
end


