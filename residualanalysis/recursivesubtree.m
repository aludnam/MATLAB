function [ixZ_vec, endleaves_vec] = recursivesubtree(Z, ixZ, ixZ_vec, endleaves_vec)
% [ixZ_vec, endleaves_vec] = recursivesubtree(Z, ixZ, ixZ_vec, endleaves_vec)
% recursively finds all subnodes and leaves from the node with index ixZ in
% matrix Z
% initialized as:
% [ixZ_vec, endleaves_vec] = recursivesubtree(Z, ixZ,[],[])
% Z output from 'linkage' function 
% ixZ -  index of the node in matrix Z (eg ixZ=(find(Z(:,3)==distance_value)))
% Z = linkage(ccds,'complete');
% [H,T,perm] = dendrogram(Z,0, 'colorthreshold',0.3);
% then subtree of node with ixZ = 50 in dendrogram
% [H,T,perm] = dendrogram(Z,0);
% can be enhanced by eg
% set(H(recursivesubtree(Z,50,[],[])),'LineWidth',2)
ixZ_vec(end+1) = ixZ;
value = Z(ixZ, 1:2) - (length(Z)+1);
if value(1) > 0
    ixZ = value(1);
    ixZ_vec(end+1) = ixZ;
    [ixZ_vec, endleaves_vec] = recursivesubtree(Z, ixZ, ixZ_vec, endleaves_vec);
else 
    endleaves_vec(end+1) = Z(ixZ, 1);
end
if value(2) > 0
    ixZ = value(2);
    ixZ_vec(end+1) = ixZ;
    [ixZ_vec, endleaves_vec] = recursivesubtree(Z, ixZ, ixZ_vec, endleaves_vec);
else 
    endleaves_vec(end+1) = Z(ixZ, 2);
end
