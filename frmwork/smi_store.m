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
    end
    
    
    methods
        function obj = smi_store(value)
            % Construct an smi_store instance
            %
            %   obj = smi_store(value);
            %       constructs an smi_store object with the value to 
            %       be captured.            
            %
            
            outlet.name = 'value';
            outlet.type = class(value);
                                    
            vs = size(value);
            if numel(vs) == 2 && vs(2) == 1
                vs = vs(1);
            end
            outlet.size = vs;
            
            obj = obj@smi_func([], outlet);
            
            obj.value = value;

        end        
                
        function tf = test_slots(obj, ~, ~) %#ok<MANU>
            tf = 1;            
        end
                
        function v = evaluate(obj, outflags) %#ok<INUSD>
            v = obj.value;
        end

    end
end