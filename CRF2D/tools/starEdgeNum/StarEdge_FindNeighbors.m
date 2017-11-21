function [nbrs,edges] = Star_FindNeighbors(V,E,node)
% Given [V,E] from a Star edge representation, returns:
%
%   nbrs: The node numbers of the neighbors of 'node';
%   edges: The edge numbers of the edges between 'node' and its neighbors
%       (in the same order)

edges = V(node):V(node+1)-1;
nbrs = E(edges);