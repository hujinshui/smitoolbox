function S = lp_problem(c, A, b, Aeq, beq, l, u)
% Constructs a struct to represent a generic Linear programming problem
%
%   A Linear programming problem is formulated as:
%
%       minimize f(x) = c' * x;
%
%       s.t.    bl <= A * x <= bu
%               Aeq * x = beq
%               l <= x <= u
%
%   S = lp_problem(c, A, b);
%   S = lp_problem(c, A, b, Aeq, beq); 
%   S = lp_problem(c, A, b, Aeq, beq, l, u);
%       constructs the problem.
%
%       Input arguments:
%       - c:    the objective coefficients
%       - A:    the coefficient matrix of the inequalities
%       - b:    the bounds of the inequalities, which can be in form of
%               bu or {bl, bu}.
%       - Aeq:  the coefficient matrix of the equalities
%       - beq:  the right hand side values of the equalities
%       - l:    the lower bound of the solution entries
%       - u:    the upper bound of the solution entries
%
%       Note that when some contraints are not used, the corresponding
%       inputs can be omitted or set to empty.
%
%       The output is a struct, which has a field type = 'lp', and the
%       following fields: d, c, A, bl, bu, Aeq, beq, l, u.
%
%   Remarks
%   -------
%       - x0 is the initial guess of the solution. If can be empty or
%         omitted, when it is available. Note that x0 serves just a hint,
%         whether it is used depends on particular solver chosen to solve
%         the problem.
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 25, 2010
%

%% main

error(nargchk(3, 8, nargin));

% c
if ~(isnumeric(c) && isvector(c) && isreal(c) && ~isempty(c))
    error('lp_problem:invalidarg', 'c should be a non-empty real vector.');         
end
if ~isa(c, 'double')
    c = double(c);
end
if size(c, 2) > 2
    c = c.';
end
d = length(c);

% constraints
[A, bl, bu] = check_lin_iec(d, A, b, 'lp_problem');

if nargin >= 4
    [Aeq, beq] = check_lin_eqc(d, Aeq, beq, 'lp_problem');
else
    Aeq = [];
    beq = [];
end

% bounds

if nargin >= 6
    if nargin < 7
        u = [];
    end
    [l, u] = check_bnds(d, l, u, 'lp_problem');
else
    l = [];
    u = [];
end

    
% output

S = struct( ...
    'type', 'lp', ...
    'd', d, ...
    'c', c, ...
    'A', A, 'bl', bl, 'bu', bu, ...
    'Aeq', Aeq, 'beq', beq, ...
    'l', l, 'u', u);


