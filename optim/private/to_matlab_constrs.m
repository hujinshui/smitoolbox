function [A, b, Aeq, beq] = to_matlab_constrs(prob)
% convert generic constraint representation to matlab constraints
%

A0 = prob.A;
bl = prob.bl;
bu = prob.bu;

if isempty(A0)
    A = [];
    b = [];
    
else
    if isempty(bl)
        A = A0;
        b = bu;
    elseif isempty(bu)
        A = -A0;
        b = -bl;
    else
        fl = find(isfinite(bl));
        fu = find(isfinite(bu));
        
        A = [A0(fu, :); -A0(fl, :)];
        b = [bu(fu); -bl(fl)];
    end
end

Aeq = prob.Aeq;
beq = prob.beq;
