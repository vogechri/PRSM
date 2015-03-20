function [disp, uDisp, vDisp, dDisp] = convert3Dto2D(ref, cam, d, wx, wy, wz, j)

% which cam should be used
if ~exist('j', 'var')
  j = 1;
end

[u, lambda_uv] = Compute_lambdaUV(ref.R, size(d,2), size(d,1));

[~, ~, ulx_ ] = computeImageProperties( ref.R, cam(j).R, u, cam(j).epi, ...
  lambda_uv, d, wx, wy, wz, 0, 0);

%%% what is correct here ???, Hog Ball needs d+wz, check with transforming
%%% back and forth, DEPENDS on the interpretation in 3d space what is right
%%% last row of K: - - - -> + else - wz
%dDisp = depth2disp(cat(3, ulx_, ones(size(d))), d-wz, cam(j).epi, ref.R, cam(1).R, lambda_uv, cam(j).q2dx);
dDisp = depth2disp(cat(3, ulx_, ones(size(d))), d+wz, cam(j).epi, ref.R, cam(1).R, lambda_uv, cam(j).q2dx);
disp  = depth2disp(u, d, cam(j).epi, ref.R, cam(j).R, lambda_uv, cam(j).q2dx);
uDisp = ulx_(:,:,1) - u(:,:,1);
vDisp = ulx_(:,:,2) - u(:,:,2);
dDisp = dDisp - disp;

end