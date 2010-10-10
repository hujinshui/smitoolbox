function [x, info] = lp_solve_ex1(solver, printLevel)
% An example showing how to use lp_solve to solve LP problem
%
%   lp_solve_ex1(solver, printLevel);
%

% Created by Dahua Lin, on Sep[ 25, 2010
%


% construct problem

c = [4 5 6];

A = [1 1 0; 1 -1 0; 7 12 0];
bl = [11 -inf 35]';
bu = [inf 5 inf]';

Aeq = [-1 -1 1];
beq = 0;

prb= lp_problem(c, A, {bl, bu}, Aeq, beq, 0, []);

% solve

[x, info] = lp_solve(prb, solver, [], printLevel);

% Note the ground-truth answer is [8, 3, 11]';
