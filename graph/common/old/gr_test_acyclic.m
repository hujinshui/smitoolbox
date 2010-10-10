function [tf, tord] = gr_test_acyclic(G)
% Test whether a graph is acylic and compute topological order if yes
%
%   tf = gr_test_acyclic(G);
%       tests whether the graph G (either an affinity matrix or
%       a mgraph struct) is an acylic graph. 
%
%   [tf, tord] = gr_test_acylic(G);
%       additionally returns the topological order of nodes.
%

%% verify input

G = mgraph(G);

%% main

if nargout <= 1
    tf = gr_test_acyclic_cimp(G);
else
    [tf, tord] = gr_test_acyclic_cimp(G);
end

