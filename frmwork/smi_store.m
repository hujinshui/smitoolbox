classdef smi_store < smi_func
    % A trivial smi function that just stores a fix value
    %
    % Each smi_store instance captures a fixed value, which can be
    % of any type and any size. 
    %
    % It has one output slot, named value, which simply outputs the
    % value that it contains.
    %
    
    % Created by Dahua Lin, on Aug 26, 2010
    %
    
    
    properties(GetAccess='public', SetAccess='private')
        value;
        vsize;
    end
    
    
    methods
        function obj = smi_store(value)
            % Construct an smi_store instance
            %
            %   obj = smi_store(value);
            %       constructs an smi_store object with the value to 
            %       be captured.            
            %
            
            obj = obj@smi_func(0, 1);
            
            obj.value = value;
            vs = size(value);
            if numel(vs) == 2 && vs(2) == 1
                vs = vs(1);
            end
            obj.vsize = vs;
        end        
        
        function [dir, id] = get_slot_id(obj, name) %#ok<MANU>
            if strcmp(name, 'value')
                dir = 'out';
                id = 1;
            else
                error('smi_store:invalidarg', ...
                    'Invalid port name for smi_store.');
            end
        end
        
        function info = get_slot_info(obj, dir, i)
            
            if strcmp(dir, 'out') && i == 1
                info.name = 'value';
                info.type = class(obj.value);
                info.size = obj.vsize;
            else
                error('smi_store:invalidarg', ...
                    'smi_store has only one output slot.');
            end
                
        end 
        
        function tf = test_slots(obj, inflags, outflags) %#ok<MANU,INUSL>
            tf = 1;            
        end
                
        function v = evaluate(obj, outflags) %#ok<INUSD>
            v = obj.value;
        end

    end
end