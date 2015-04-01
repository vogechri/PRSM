% returns the gradients of the scalar field u
% these are defined on the dual elements of the grid
% which is here (linear basis elements) triangles
% note that there are roughly twice as many triangles than vertices 
function [dx, dy, d_xy, dixy] = getEdgeWeightsForSegmentation (u, varWindow, maxDiff)

minVal    = 0.05;
% simple version:

oldWay =1;
if oldWay

g = edge_Tensor(u);

dx = g(:,:,1);dx = dx(:,1:end-1);
dy = g(:,:,4);dy = dy(1:end-1,:);

d_xy = (g(:,:,1) + g(:,:,4) + g(:,:,2) + g(:,:,3))./2;d_xy = d_xy(1:end-1,1:end-1);
dixy = (g(:,:,1) + g(:,:,4) - g(:,:,2) - g(:,:,3))./2;dixy = dixy(1:end-1,1:end-1);

dx = max( minVal, dx );
dy = max( minVal, dy );

dixy = max( minVal, dixy );
d_xy = max( minVal, d_xy );

else

dx   = edge_Tensor_Dir(u, 0.5, 0,   [1,0]);
dy   = edge_Tensor_Dir(u, 0.0, 0.5, [0,1]);
d_xy = edge_Tensor_Dir(u, 0.5, 0.5, [1./sqrt(2),1./sqrt(2)]);
dixy = edge_Tensor_Dir(u, 0.5, 0.5, [1./sqrt(2),-1./sqrt(2)]);

dx = max( minVal, dx );
dy = max( minVal, dy );

dixy = max( minVal, dixy );
d_xy = max( minVal, d_xy );

end

return;
