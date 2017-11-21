function [E] = StarEdge_EdgeNum(V,E,n1,n2)
% Returns the edge number of the edge between nodes n1 and n2
% using the [V,E] from a star edge representation,
%
% If the node does not exist, returns an empty matrix
%
% Although simple, this function is innefficient and will be slow
%   if n1 has a lot of edges:
%   - avoid calling it if possible (use StarEdge_FindNeighbors instead)
%   - it could be made more efficient with an auxiliary data structure

[nbrs,edges] = StarEdge_FindNeighbors(V,E,n1);
E = edges(nbrs == n2);
