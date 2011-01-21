function S = gurobi_params(varargin)
% Make a Gurobi parameter struct
%
%   S = gurobi_params('name1', value1, 'name2', value2, ...);
%
%       creates a Gurobi parameter struct to control the solver.
%       The parameters are input in a name value list form.
%
%       Please refer to the reference manual of Gurobi for a list of
%       all available parameters.
%
%       In addition, we add the additional parameters as follows:
%
%       - 'Display':   integer 0 or 1, to control whether to display
%                      procedural information. (default = 0)
%
%       Note that the parameter names are case-sensitive.
%

% Created by Dahua Lin, on Jan 21, 2011
%


%% Prepare the type map

tmap = containers.Map('KeyType', 'char', 'ValueType', 'char');

% Termination

tmap('BarIterLimit') = 'int';
tmap('Cutoff') = 'dbl';
tmap('IterationLimit') = 'dbl';
tmap('NodeLimit') = 'dbl';
tmap('SolutionLimit') = 'int';
tmap('TimeLimit') = 'dbl';

% Tolerance

tmap('FeasibilityTol') = 'dbl';
tmap('IntFeasTol') = 'dbl';
tmap('MarkowitzTol') = 'dbl';
tmap('MIPGap') = 'dbl';
tmap('MIPGapAbs') = 'dbl';
tmap('OptimalityTol') = 'dbl';
tmap('PSDTol') = 'dbl';

% Simplex

tmap('Method') = 'int';
tmap('PerturbValue') = 'dbl';
tmap('ObjScale') = 'dbl';
tmap('ScaleFlag') = 'int';
tmap('SimplexPricing') = 'int';
tmap('Quad') = 'int';
tmap('NormAdjust') = 'int';

% Barrier

tmap('BarConvTol') = 'dbl';
tmap('BarCorrectors') = 'int';
tmap('BarOrder') = 'int';
tmap('Crossover') = 'int';
tmap('CrossoverBasis') = 'int';

% MIP

tmap('Heuristics') = 'dbl';
tmap('ImproveStartGap') = 'dbl';
tmap('ImproveStartTime') = 'dbl';
tmap('NodefileDir') = 'str';
tmap('NodefileStart') = 'dbl';
tmap('NodeMethod') = 'int';
tmap('PumpPasses') = 'int';
tmap('RINS') = 'int';
tmap('SubMIPNodes') = 'int';
tmap('Symmetry') = 'int';
tmap('VarBranch') = 'int';
tmap('MIPFocus') = 'int';
tmap('SolutionNumber') = 'int';

% MIP cuts

tmap('Cuts') = 'int';

tmap('CliqueCuts') = 'int';
tmap('CoverCuts') = 'int';
tmap('FlowCoverCuts') = 'int';
tmap('FlowPathCuts') = 'int';
tmap('GUBCoverCuts') = 'int';
tmap('ImpliedCuts') = 'int';
tmap('MIPSepCuts') = 'int';
tmap('MIRCuts') = 'int';
tmap('ModKCuts') = 'int';
tmap('ZeroHalfCuts') = 'int';
tmap('NetworkCuts') = 'int';
tmap('SubMIPCuts') = 'int';

tmap('CutAggPasses') = 'int';
tmap('GomoryPasses') = 'int';

% Others

tmap('Aggregate') = 'int';
tmap('AggFill') = 'int';
tmap('DisplayInterval') = 'int';
tmap('IISMethod') = 'int';
tmap('LogFile') = 'int';
tmap('OutputFlag') = 'int';
tmap('PreCrush') = 'int';
tmap('PreDepRow') = 'int';
tmap('PreDual') = 'int';
tmap('PreMIQPMethod') = 'int';
tmap('PrePasses') = 'int';
tmap('Presolve') = 'int';
tmap('ResultFile') = 'str';
tmap('Threads') = 'int';

% Additionals

tmap('Display') = 'int';


%% make the parameter struct

S = [];
if ~isempty(varargin)
    names = varargin(1:2:end);
    vals = varargin(2:2:end);
    
    if ~(iscellstr(names) && length(names) == length(vals))
        error('gurobi_params:invalidarg', ...
            'The name/value list is invalid.');
    end
    
    for i = 1 : length(names)
        
        name = names{i};
        v = vals{i};
        
        if ~isKey(tmap, name)
            error('gurobi_params:invalidarg', ...
                'Unsupported parameter name %s', name);
        end
        
        ty = tmap(name);
        
        switch ty
            case 'int'
                if ~(isscalar(v) && isnumeric(v) && v == fix(v))
                    error('gurobi_params:invalidarg', ...
                        'The value of parameter %s should be an integer', name);
                end
                v = int32(v);
                
            case 'dbl'
                if ~(isscalar(v) && isnumeric(v) && isreal(v))
                    error('gurobi_params:invalidarg', ...
                        'The value of parameter %s should be a real scalar', name);
                end
                v = double(v);
                
            case 'str'
                if ~(ischar(v) && ndims(v) == 2 && size(v,1) == 1)
                    error('gurobi_params:invalidarg', ...
                        'The value of parameter %s should be a string', name);
                end                                
        end
        
        S.(name) = v;                
    end
end


