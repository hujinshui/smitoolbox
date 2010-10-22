function [A, bl, bu] = to_tomlab_constrs(prob)
% convert generic constraint representation to tomlab constraints
%
% Note that mosek also a form like tomlab's.
%

A = prob.A;
bl = prob.bl;
bu = prob.bu;

if ~isempty(A);
    if isempty(bl)
        bl = -inf(size(A,1), 1);
    end
    if isempty(bu)
        bu = inf(size(A,1), 1);
    end
end

if ~isempty(prob.Aeq)
    A = [A; prob.Aeq];
    bl = [bl; prob.beq];
    bu = [bu; prob.beq];
end





