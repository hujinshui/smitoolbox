function sol = dpmm_demo()
% A program to demonstrate the use of DPMM
%
%   dpmm_demo;
%

% Created by Dahua Lin, on Sep 21, 2011
%

%% Prepare data

d = 2;
Cx = pdmat('s', d, 1);
Cu = pdmat('s', d, 5^2);

centers = [-1 0; 1 0; 0 1]' * 5;
K0 = size(centers, 2);

n = 1000;

Xs = cell(1, K0);
for k = 1 : K0
    Xs{k} = gsample(centers(:,k), Cx, n);
end
X = [Xs{:}];

%% Construct model

gm = gaussgm(Cx, Cu);
model = gauss_npmodel(gm, zeros(d, 1), X);

%% Construct solution

alpha = 1;
sol = dpmm_solution(model, alpha);

T = 10;
for t = 1 : T
    sol.update_labels();
    sol.update_atoms();
end


%% Visualize

% figure;
% title('DPMM (Gauss) Demo');
% plot(X(1,:), X(2,:), '.');
% axis equal;
% 
% A = sol.get_atoms();
% A = [A{:}];
% cnts = sol.get_atom_counts();
% 
% A = A(:, cnts > n / 2);
% 
% hold on;
% plot(A(1,:), A(2,:), 'r+', 'MarkerSize', 20, 'LineWidth', 2);
% hold off;

