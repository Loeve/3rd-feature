function [V,E,nedges] = StarEdge_MakeEdgeNums_Lattice2(nr,nc)
% Returns the [V,E] vectors representing 2D lattice
% edges in the 'star' edge representation
%
% V ranges from 1:nnodes+1, and points to the starting point of 
%   each node's neighbors in E.  The last element points past the end of E
%   to make referencing into E with V easier.
%
% The nodes are assumed to be labeled in COLUMN MAJOR order
% ie. a 3 by 2 lattice would have node labelings: [1 4 7;2 5 8;3 6 9]
%
% E ranges from 1:nedges, and contains the node numbers of 
%   of the neighbors of V
%
% This function only supports connectivity of 4
%   (r=1 von Neumann neighbors)
% 
% If you are going to call this function many times, consider 
% vectorizing or mexing

node = 1;
edge = 1;
for j = 1:nc
    for i = 1:nr
        V(node) = edge;
        
        % Left neighbor
        if legal(i,j-1,nr,nc)
            E(edge) = node-nr;
            edge = edge+1;
        end
        
        % Top neighbor
        if legal(i-1,j,nr,nc)
            E(edge) = node-1;
            edge = edge+1;
        end
        
        % Bottom neigbor
        if legal(i+1,j,nr,nc)
            E(edge) = node+1;
            edge = edge+1;
        end
        
        % Right neighbor
        if legal(i,j+1,nr,nc)
            E(edge) = node+nr;
            edge = edge+1;
        end
        
        
        % Increment node counter
        node = node + 1;
    end
end
V(node)=edge;
nedges = edge-1;

function [bool] = legal(i,j,nrows,ncols)
    if i >= 1 && j >= 1 && i <= nrows && j <=ncols
       bool = 1;
    else
        bool = 0;
    end
