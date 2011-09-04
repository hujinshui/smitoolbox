function [Aeq, beq] = check_lin_eqc(d, Aeq, beq, funname)
% Check the validity of linear equality constraints
%
%   [Aeq, beq] = check_lin_eqc(d, Aeq, beq, funname);
%

if ~(isfloat(Aeq) && ndims(Aeq) == 2)
    error([funname ':invalidconstr'], ...
        'The equality coefficient matrix Aeq should be a numeric matrix.');
end

[m, n] = size(Aeq);
if n ~= d
    error([funname ':invalidconstr'], ...
        'The requirement size(Aeq, 2) == d is not satisfied.');
end

if ~(isfloat(beq) && ndims(beq) == 2 && size(beq, 2) == 1)
    error([funname ':invalidconstr'], ...
        'The right-hand-side of equality beq should be a numeric column vector.');
end
if size(beq,1) ~= m
    error([funname ':invalidconstr'], ...
        'The right-hand-side of equality is inconsistent with Aeq.');
end
