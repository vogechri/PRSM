%% This function computes the real world 3d positions and the 3d flow. 
%% The latter is fairly easy with the flow as input, but the positions are 
%% computed using the depth and the camera matrix of the reference camera.
%% Note that the matrix R_l has to be the camera matrix adjusted to the 
%% current image size
%%
%% R_l       : the camera matrix (first 3x3 part)
%% u         : the image (pixel) coordinates
%% lambda_uv : the scalars turning the image coordiantes (u) into 3d
%%             coordinates (basically the focal distance)
%% d : the depth
%% wX : the flow in x direction in world units
%% wY : the flow in x direction in world units
%% wZ : the flow in x direction in world units
%%
%% pos3d  : the 3d positions 
%% flow3d : the flow in 3d 
%%
function [pos3d, flow3d] = ...
  Compute_3D_Position(R_l, u, lambda_uv, d, wx, wy, wz )

[M N K] = size(u);

%Viewport = [ (N)/2 0 (N+1)/2 ; 0 (M)/2 (M+1)/2 ; 0  0  1];
%iViewport = inv (Viewport);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flow3d = zeros(M,N,K);
pos3d  = zeros(M,N,K);

% OR:
Q = inv (R_l);
%Q = Q * iViewport;

qu = zeros(M,N,K);
qu(:,:,1) = lambda_uv .* (Q(1,1) * u(:,:,1) + Q(1,2) * u(:,:,2) + Q(1,3) * u(:,:,3));
qu(:,:,2) = lambda_uv .* (Q(2,1) * u(:,:,1) + Q(2,2) * u(:,:,2) + Q(2,3) * u(:,:,3));
qu(:,:,3) = lambda_uv .* (Q(3,1) * u(:,:,1) + Q(3,2) * u(:,:,2) + Q(3,3) * u(:,:,3));

pos3d(:,:,1) = d .* qu(:,:,1);
pos3d(:,:,2) = d .* qu(:,:,2);
pos3d(:,:,3) = d .* qu(:,:,3);

flow3d(:,:,1) = wx;
flow3d(:,:,2) = wy;
flow3d(:,:,3) = wz;

end