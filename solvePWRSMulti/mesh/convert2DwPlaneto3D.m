function [ p3d ] = convert2DwPlaneto3D( N_,  p2d_, K_linv)
% N_ is an input plane, represented as a normal scaled by depth, p2d_ is a
% list of pixel locations, K_linv is the inverse camera intrinsics.
% Intersects the ray represented by the 2D pixel locations, normalized by
% the camera intrinsics, and the plane represented by the input normal to
% find corresponding metric 3D points, these are returned in p3d
%
% based on reconstruc3DFlowHom

p2d_(:,3) = 1;           % Homogenize
p2d_ = (K_linv * p2d_');  % normalize out intrinsics
d = abs(1./(p2d_'*N_));  % compute depth from ray/plane intersection
p3d = p2d_' .* d;        % use depth to metric-ize

end

