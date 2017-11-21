function [n1,n2] = StarEdge_GetNodes(V,E,edge)
% Returns the node numbers associated with 'edge'

n1 = 1;
while(V(n1+1) <= edge)
   n1 = n1 + 1; 
end

n2 = E(edge);