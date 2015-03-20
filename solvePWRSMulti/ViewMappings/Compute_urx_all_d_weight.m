% compute the position in the right image at timestep t+1 (in image coords.)
% and the weight for the data term equations for the depth wD. Further an
% occlusion maop is generated (oMap). 
%
% see 3D scene flow proposal
%
% ur    : pixel posistion in right image.
% durx_dd : derivative of pixel position by depth
% durx_dw1 : derivative of pixel position by xDirection of the flow
% durx_dw2 : derivative of pixel position by yDirection of the flow
% durx_dw3 : derivative of pixel position by zDirection of the flow
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
% wx: the current wx estimate
% wy: the current wy estimate
% wz: the current wz estimate
% doOccMap: compute occlusion map : yes or no
%
function [urx, durx_dd, durx_dw1, durx_dw2, durx_dw3, oMap, urxDepth] = ...
  Compute_urx_all_d_weight(R_l, R_r, u, epi, lambda_uv, d, wx, wy, wz, doOMap, ...
  SegImg )

if ~exist('SegImg','var')
  SegImg = 0;
end

[M, N, K] = size(u);

qu = zeros(M,N,K);
rw = zeros(M,N,K);
zz = zeros(M,N,K);

urx = zeros(M,N,2);

Q = R_r * inv (R_l);

qu(:,:,1) = lambda_uv .* (Q(1,1) * u(:,:,1) + Q(1,2) * u(:,:,2) + Q(1,3) * u(:,:,3));
qu(:,:,2) = lambda_uv .* (Q(2,1) * u(:,:,1) + Q(2,2) * u(:,:,2) + Q(2,3) * u(:,:,3));
qu(:,:,3) = lambda_uv .* (Q(3,1) * u(:,:,1) + Q(3,2) * u(:,:,2) + Q(3,3) * u(:,:,3));

rw(:,:,1) = R_r(1,1) * wx + R_r(1,2) * wy + R_r(1,3) * wz;
rw(:,:,2) = R_r(2,1) * wx + R_r(2,2) * wy + R_r(2,3) * wz;
rw(:,:,3) = R_r(3,1) * wx + R_r(3,2) * wy + R_r(3,3) * wz;

zz(:,:,1) = epi(1) + d .* qu(:,:,1) + rw(:,:,1);
zz(:,:,2) = epi(2) + d .* qu(:,:,2) + rw(:,:,2);
zz(:,:,3) = epi(3) + d .* qu(:,:,3) + rw(:,:,3);

nn = zz(:,:,3).^2;

urx (:,:,1) = zz(:,:,1) ./ zz(:,:,3);
urx (:,:,2) = zz(:,:,2) ./ zz(:,:,3);

% derivatives all denom * (d/d_f nom)- nom * (d/d_f denom) / denom^2
durx_dd (:,:,1) = qu(:,:,1) .* zz(:,:,3) - qu(:,:,3) .* zz(:,:,1);
durx_dd (:,:,2) = qu(:,:,2) .* zz(:,:,3) - qu(:,:,3) .* zz(:,:,2);

durx_dw1 (:,:,1) = R_r(1,1) .* zz(:,:,3) - R_r(3,1) .* zz(:,:,1);
durx_dw1 (:,:,2) = R_r(2,1) .* zz(:,:,3) - R_r(3,1) .* zz(:,:,2);

durx_dw2 (:,:,1) = R_r(1,2) .* zz(:,:,3) - R_r(3,2) .* zz(:,:,1);
durx_dw2 (:,:,2) = R_r(2,2) .* zz(:,:,3) - R_r(3,2) .* zz(:,:,2);

durx_dw3 (:,:,1) = R_r(1,3) .* zz(:,:,3) - R_r(3,3) .* zz(:,:,1);
durx_dw3 (:,:,2) = R_r(2,3) .* zz(:,:,3) - R_r(3,3) .* zz(:,:,2);

durx_dd (:,:,1) = durx_dd(:,:,1) ./ nn;
durx_dd (:,:,2) = durx_dd(:,:,2) ./ nn;

durx_dw1 (:,:,1) = durx_dw1(:,:,1) ./ nn;
durx_dw1 (:,:,2) = durx_dw1(:,:,2) ./ nn;

durx_dw2 (:,:,1) = durx_dw2(:,:,1) ./ nn;
durx_dw2 (:,:,2) = durx_dw2(:,:,2) ./ nn;

durx_dw3 (:,:,1) = durx_dw3(:,:,1) ./ nn;
durx_dw3 (:,:,2) = durx_dw3(:,:,2) ./ nn;

urxDepth = zz(:,:,3)./lambda_uv;

%%%%%%%%%%%%%%%%%%%%%%%%
% compute occlusion Map:
oMap = ones(M,N);
if doOMap
%  oMap = buildOccMap ( urx, zz(:,:,3) );
  oMap = buildOccMapSubPix ( urx, urxDepth, SegImg );
end
% also occluded / out of image == invisible
oMap (d + wz < 0) = 0;
