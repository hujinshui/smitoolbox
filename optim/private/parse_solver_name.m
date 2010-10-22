function [pname, sname] = parse_solver_name(name, fname)
% parse the solver name into package name and solver name
%

i = find(name == '.', 1);

if isempty(i)
    error([fname ':invalidarg'], ...
        'The solver name %s is invalid', name);
end

pname = name(1:i-1);
sname = name(i+1:end);

