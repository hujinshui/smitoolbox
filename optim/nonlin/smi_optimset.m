function options = smi_optimset(op0, varargin)
% Set optimization options (for SMI optimization functions)
%
%   options = smi_optimset(funcname);
%
%       Gets the default optimization setting for a smi optimization
%       function (such as 'gdfmin', 'bfgsfmin', and 'newtonfmin').
%
%   options = smi_optimset(funcname, 'name1', value1, 'name2', value2, ...);
%   options = smi_optimset(options0, 'name1', value1, 'name2', value2, ...);
%
%       Sets the optimization setting. It starts with the default setting
%       of a particular function or the supplied initial setting, and
%       updates it using the values in the ensuing name/value pairs.
%

%   Created by Dahua Lin, on Jan 5, 2011
%

%% Initialize options

if ischar(op0) || isa(op0, 'function_handle')
    options = feval(op0, 'options');
else
    if ~(isstruct(op0) && numel(op0) == 1)
        error('smi_optimset:invalidarg', ...
            'The first argument to smi_optimset is invalid.');
    end
    options = op0;
end

    
%% Update

if ~isempty(varargin)
    
    ns = varargin(1:2:end);
    vs = varargin(2:2:end);
    
    if numel(ns) ~= numel(vs) || ~iscellstr(ns)
        error('smi_optimset:invalidarg', ...
            'The input option list is invalid.');
    end
    
    for i = 1 : numel(ns)
        
        nam = ns{i};
        v = vs{i};
        
        switch lower(nam)
            case 'maxiter'
                if ~(isnumeric(v) && isscalar(v) && v >= 1)
                    error('smi_optimset:invalidarg', ...
                        'MaxIter should be a positive integer scalar.');
                end
                options.MaxIter = v;
                
            case 'tolfun'
                if ~(isfloat(v) && isscalar(v) && v > 0)
                    error('smi_optimset:invalidarg', ...
                        'TolFun should be a positive scalar.');
                end
                options.TolFun = v;
                
            case 'tolx'
                if ~(isfloat(v) && isscalar(v) && v > 0)
                    error('smi_optimset:invalidarg', ...
                        'TolX should be a positive scalar.');
                end
                options.TolX = v;
                
            case 'display'
                mon = optim_mon(v);
                options.Monitor = mon;
                
            case 'monitor'
                if ~isobject(v)
                    error('smi_optimset:invalidarg', ...
                        'Monitor should be an object.');
                end
                options.Monitor = v;
                
            case 'initinvhess'
                if ~(isfloat(v) && ndims(v) == 2 && size(v,1) == size(v,2))
                    error('smi_optimset:invalidarg', ...
                        'InitInvHess should be a numeric square matrix.');
                end
                options.InitInvHess = v;
                
            case 'directnewton'
                if ~((isnumeric(v) || islogical(v)) && isscalar(v))
                    error('smi_optimset:invalidarg', ...
                        'DirectNewton should be a numeric/logical scalar.');
                end
                options.DirectNewton = logical(v);
                
            otherwise
                error('smi_optimset:invalidarg', ...
                    'Unsupported optimization option %s', nam);
        end
                
    end
        
end

