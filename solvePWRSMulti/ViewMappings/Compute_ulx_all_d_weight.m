% compute the position in the reference image at timestep t+1 (in image coords.)
% and the weight for the data term equations for the depth wD. Further an
% occlusion maop is generated (oMap).
%
% see 3D scene flow proposal
%
% ur    : pixel posistion in right image.
% dulx_d : derivative of pixel posisiotn by depth
% dulx_w1 : derivative of pixel posisiotn by xDirection of the flow
% dulx_w2 : derivative of pixel posisiotn by yDirection of the flow
% dulx_w3 : derivative of pixel posisiotn by zDirection of the flow
% oMap the occlusion map for this image, parameterized by the reference
% camera
%
% R_l : the reference camera matrix
% u   : the pixel coordinates in the reference view
% lambda_uv : the scalars to transfer image into 3d coordiantes for the
% reference camera
% d: the current depth estimate
% wx: the current wx estimate
% wy: the current wy estimate
% wz: the current wz estimate
% doOccMap: compute occlusion map : yes or no
%
function [ulx, dulx_d, dulx_w1, dulx_w2, dulx_w3, oMap, ulxDepth] = ...
  Compute_ulx_all_d_weight(R_l, u, lambda_uv, d, wx, wy, wz, doOMap, ...
  SegImg )

if ~exist('SegImg','var')
  SegImg = 0;
end

[M, N, K] = size(u);

zz  = zeros(M,N,K);
rw  = zeros(M,N,K);
ulx = zeros(M,N,2);

dulx_d = zeros(M,N,2);
dulx_w1 = zeros(M,N,2);
dulx_w2 = zeros(M,N,2);
dulx_w3 = zeros(M,N,2);

uluv(:,:,1) = u(:,:,1) .* lambda_uv;
uluv(:,:,2) = u(:,:,2) .* lambda_uv;
uluv(:,:,3) = u(:,:,3) .* lambda_uv;

rw(:,:,1) = R_l(1,1) * wx + R_l(1,2) * wy + R_l(1,3) * wz;
rw(:,:,2) = R_l(2,1) * wx + R_l(2,2) * wy + R_l(2,3) * wz;
rw(:,:,3) = R_l(3,1) * wx + R_l(3,2) * wy + R_l(3,3) * wz;

zz(:,:,1) = d .* uluv(:,:,1) + rw(:,:,1);
zz(:,:,2) = d .* uluv(:,:,2) + rw(:,:,2);
zz(:,:,3) = d .* uluv(:,:,3) + rw(:,:,3);

nn = zz(:,:,3).^2;

ulx (:,:,1) = zz(:,:,1) ./ zz(:,:,3);
ulx (:,:,2) = zz(:,:,2) ./ zz(:,:,3);

% derivatives all denom * (d/d_f nom)- nom * (d/d_f denom) / denom^2
dulx_d(:,:,1) = zz(:,:,3) .* uluv(:,:,1) - zz(:,:,1) .* uluv(:,:,3);
dulx_d(:,:,2) = zz(:,:,3) .* uluv(:,:,2) - zz(:,:,2) .* uluv(:,:,3);

dulx_w1(:,:,1) = zz(:,:,3) * R_l(1,1) - zz(:,:,1) * R_l(3,1);
dulx_w1(:,:,2) = zz(:,:,3) * R_l(2,1) - zz(:,:,2) * R_l(3,1);

dulx_w2(:,:,1) = zz(:,:,3) * R_l(1,2) - zz(:,:,1) * R_l(3,2);
dulx_w2(:,:,2) = zz(:,:,3) * R_l(2,2) - zz(:,:,2) * R_l(3,2);

dulx_w3(:,:,1) = zz(:,:,3) * R_l(1,3) - zz(:,:,1) * R_l(3,3);
dulx_w3(:,:,2) = zz(:,:,3) * R_l(2,3) - zz(:,:,2) * R_l(3,3);

dulx_d(:,:,1)  = dulx_d(:,:,1)  ./ nn;
dulx_w1(:,:,1) = dulx_w1(:,:,1) ./ nn;
dulx_w2(:,:,1) = dulx_w2(:,:,1) ./ nn;
dulx_w3(:,:,1) = dulx_w3(:,:,1) ./ nn;

dulx_d(:,:,2)  = dulx_d(:,:,2)  ./ nn;
dulx_w1(:,:,2) = dulx_w1(:,:,2) ./ nn;
dulx_w2(:,:,2) = dulx_w2(:,:,2) ./ nn;
dulx_w3(:,:,2) = dulx_w3(:,:,2) ./ nn;

ulxDepth = zz(:,:,3)./lambda_uv;

oMap = ones(M,N);
if doOMap
%  oMap = buildOccMap ( ulx, zz(:,:,3) );
  oMap = buildOccMapSubPix ( ulx, ulxDepth, SegImg );
end
% also occluded / out of image == invisible
oMap (d + wz < 0) = 0;
