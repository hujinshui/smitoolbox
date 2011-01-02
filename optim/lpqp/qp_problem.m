function S = qp_problem(H, f, A, b, Aeq, beq, l, u)
% Constructs a struct to represent a generic Quadratic programming problem
%
%   A Linear programming problem is formulated as:
%
%       minimize f(x) = (1/2) * x' * H * x + f' * x;
%
%       s.t.    bl <= A * x <= bu
%               Aeq * x = beq
%               l <= x <= u
%
%   S = qp_problem(H, f, A, b);
%   S = qp_problem(H, f, A, b, Aeq, beq); 
%   S = qp_problem(H, f, A, b, Aeq, beq, l, u);
%       constructs the problem.
%
%       Input arguments:
%       - H:    the quadratic coefficient matrix
%       - f:    the linear coefficient vector
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
%       - Created by Dahua Lin, on Jan 2, 2010
%

%% main

error(nargchk(4, 8, nargin));

% H 
d = size(H,1);
if ~(isfloat(H) && isreal(H) && ndims(H) == 2 && ~isempty(H) && d==size(H,2))
    error('qp_problem:invalidarg', 'H should be a non-empty square matrix.');
end

% f
if ~(isnumeric(f) && isvector(f) && isreal(f) && ~isempty(f) && d==numel(f))
    error('qp_problem:invalidarg', 'f should be a non-empty real vector.');         
end

if ~isa(f, 'double'); f = double(f); end
if size(f, 2) > 2; f = f.'; end

% constraints
if ~isempty(A)
    [A, bl, bu] = check_lin_iec(d, A, b, 'qp_problem');
else
    A = [];
    bl = [];
    bu = [];
end

if nargin >= 4
    [Aeq, beq] = check_lin_eqc(d, Aeq, beq, 'qp_problem');
else
    Aeq = [];
    beq = [];
end

% bounds

if nargin >= 6
    if nargin < 7
        u = [];
    end
    [l, u] = check_bnds(d, l, u, 'qp_problem');
else
    l = [];
    u = [];
end

    
% output

S = struct( ...
    'type', 'qp', ...
    'd', d, ...
    'H', H, ...
    'f', f, ...
    'A', A, 'bl', bl, 'bu', bu, ...
    'Aeq', Aeq, 'beq', beq, ...
    'l', l, 'u', u);


