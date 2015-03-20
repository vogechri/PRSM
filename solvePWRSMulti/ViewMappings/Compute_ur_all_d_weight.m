% compute the position in the right image at timestep t (in image coords.)
% and the weight for the data term equations for the depth wD. Further an
% occlusion maop is generated (oMap).
%
% see 3D scene flow proposal
%
% ur    : pixel posistion in right image.
% dur_d : derivative of pixel posisiotn by depth
% oMap the occlusion map for this image, parameterized by the reference
% camera
%
% R_l : the reference camera matrix
% R_r : the camera matrix of the current image
% epi : the epipole in the cutrent image (of the reference camera)
% u   : the pixel coordinates in the reference view
% lambda_uv : the scalars to transfer image into 3d coordiantes for the
% reference camera
% d: the current depth estimate
%
function [ur, dur_d, oMap, urDepth] = ...
  Compute_ur_all_d_weight(R_l, R_r, epi, u, lambda_uv, d, doOMap, SegImg )

if ~exist('SegImg','var')
  SegImg = 0;
end

[M N K] = size(u);

qu = zeros(M,N,K);
zz = zeros(M,N,K);
ur = zeros(M,N,2);
Q = R_r * inv (R_l);

qu(:,:,1) = lambda_uv .*(Q(1,1) * u(:,:,1) + Q(1,2) * u(:,:,2) + Q(1,3) * u(:,:,3));
qu(:,:,2) = lambda_uv .*(Q(2,1) * u(:,:,1) + Q(2,2) * u(:,:,2) + Q(2,3) * u(:,:,3));
qu(:,:,3) = lambda_uv .*(Q(3,1) * u(:,:,1) + Q(3,2) * u(:,:,2) + Q(3,3) * u(:,:,3));

zz(:,:,1) = epi(1) + d .* qu(:,:,1);
zz(:,:,2) = epi(2) + d .* qu(:,:,2);
zz(:,:,3) = epi(3) + d .* qu(:,:,3);
%255 128-> 1 129: 0.015351 -2.537391 2.256235 w-comp all OK ?
ur (:,:,1) = zz(:,:,1) ./ zz(:,:,3);
ur (:,:,2) = zz(:,:,2) ./ zz(:,:,3);

nn = zz(:,:,3).^2;

dur_d (:,:,1) = qu(:,:,1) * epi(3) - qu(:,:,3) * epi(1);
dur_d (:,:,2) = qu(:,:,2) * epi(3) - qu(:,:,3) * epi(2);

dur_d (:,:,1) = dur_d (:,:,1) ./ nn;
dur_d (:,:,2) = dur_d (:,:,2) ./ nn;

urDepth = zz(:,:,3)./lambda_uv;

%%%%%%%%%%%%%%%%%%%%%%%%
% compute occlusion Map:
%%%%%%%%%%%%%%%%%%%%%%%%
oMap = ones(M,N);
if doOMap
%  oMap  = buildOccMap ( ur, zz(:,:,3) );
  oMap = buildOccMapSubPix ( ur, urDepth, SegImg );
end