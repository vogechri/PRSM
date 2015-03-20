%%% compute position in all other images and derivatives wrt. 3d motion and
%%% depth -- way too much for this implementation, we only need the pixel
%%% positions 
function [ur, dur_d, ...
  ulx, dulx_d, dulx_w1, dulx_w2, dulx_w3, ...
  urx, durx_dd, durx_dw1, durx_dw2, durx_dw3, ...
  validR, validRT, validLT ] = ...
  computeImageProperties(R_l, R_r, u, epi, lambda_uv, d0_, wx0_, wy0_, wz0_, ...
  doLinInterp2, doOccMap)

[M, N] = size (d0_);

[ur, dur_d, oMap] = Compute_ur_all_d_weight(R_l, R_r, epi, u, lambda_uv, d0_, doOccMap);

% boundary handling % promote to epi-constraints as well ???
if doLinInterp2
%  m = (ur(:,:,1) > N) | (ur(:,:,1) < 1) | (ur(:,:,2) > M) | (ur(:,:,2) < 1) | (oMap(:,:) == 0);
  m = (oMap(:,:) == 0);
else
  m = (ur(:,:,1) > N+ 0.5) | (ur(:,:,1) < 1- 0.5) | (ur(:,:,2) > M+ 0.5) | (ur(:,:,2) < 1- 0.5) | (oMap(:,:) == 0);
  %m = (ur(:,:,1) > N) | (ur(:,:,1) <= 1) | (ur(:,:,2) > M) | (ur(:,:,2) <= 1) | (oMap(:,:) == 0);
end
validR = repmat(m, [1,1,5]);

% right at t+1
[urx, durx_dd, durx_dw1, durx_dw2, durx_dw3, oMap] = ...
  Compute_urx_all_d_weight(R_l, R_r, u, epi, lambda_uv, d0_, wx0_, wy0_, wz0_, doOccMap);

% boundary handling % promote to epi-constraints as well ???
if doLinInterp2
%  m = (urx(:,:,1) > N) | (urx(:,:,1) < 1) | (urx(:,:,2) > M) | (urx(:,:,2) < 1) | (oMap(:,:) == 0);
  m = (oMap(:,:) == 0);
else
  m = (urx(:,:,1) > N+ 0.5) | (urx(:,:,1) < 1- 0.5) | (urx(:,:,2) > M+ 0.5) | (urx(:,:,2) < 1- 0.5) | (oMap(:,:) == 0);
%  m = (urx(:,:,1) > N) | (urx(:,:,1) <= 1) | (urx(:,:,2) > M) | (urx(:,:,2) <= 1) | (oMap(:,:) == 0);
end
validRT = repmat(m, [1,1,5]);

% left at t+1
[ulx, dulx_d, dulx_w1, dulx_w2, dulx_w3, oMap] = ...
  Compute_ulx_all_d_weight(R_l, u, lambda_uv, d0_, wx0_, wy0_, wz0_, doOccMap );
% boundary handling % promote to epi-constraints as well ???

if doLinInterp2
  m = (oMap(:,:) == 0);
else
  m = (ulx(:,:,1) > N + 0.5) | (ulx(:,:,1) < 1 - 0.5) | (ulx(:,:,2) > M + 0.5) | (ulx(:,:,2) < 1 - 0.5) | (oMap(:,:) == 0);
  %m = (ulx(:,:,1) > N) | (ulx(:,:,1) <= 1) | (ulx(:,:,2) > M) | (ulx(:,:,2) <= 1) | (oMap(:,:) == 0);
end
validLT = repmat(m, [1,1,5]);
