function [d, wx, wy, wz, pos3d, pos3d_t1] = convert2Dto3D(ref, cam, d_, wx_, wy_, wz_, j)

if ~exist('j', 'var')
  j = 1;
end

%maxDisp = -0.01;% right if rectified only!
maxDisp = max(size(d_)); % stupid value

[u, lambda_uv] = Compute_lambdaUV(ref.R, size(d_,2), size(d_,1));

ut1 = u; ut1(:,:,1) = ut1(:,:,1) + wx_; ut1(:,:,2) = ut1(:,:,2) + wy_;
% depth independant of position in image if rectified
depthT1 = disp2depth( min(d_+wz_, maxDisp), ut1, cam(j).epi, ref.R, cam(j).R, lambda_uv, cam(j).q2dx); 
depthT  = disp2depth( min(d_, maxDisp), u, cam(j).epi, ref.R, cam(j).R, lambda_uv, cam(j).q2dx);

% d_t1 must be interpolated is d + depth(disp (ut1))
pos3d_t1 = Compute_3D_Position(ref.R, ut1, lambda_uv, depthT1, wx_, wy_, wz_ );
pos3d    = Compute_3D_Position(ref.R, u, lambda_uv, depthT, wx_, wy_, wz_ );

% probe:
%{
Q = ref.R;
qu(:,:,1) = (Q(1,1) * pos3d_t1(:,:,1) + Q(1,2) * pos3d_t1(:,:,2) + Q(1,3) * pos3d_t1(:,:,3));
qu(:,:,2) = (Q(2,1) * pos3d_t1(:,:,1) + Q(2,2) * pos3d_t1(:,:,2) + Q(2,3) * pos3d_t1(:,:,3));
qu(:,:,3) = (Q(3,1) * pos3d_t1(:,:,1) + Q(3,2) * pos3d_t1(:,:,2) + Q(3,3) * pos3d_t1(:,:,3));
qu(:,:,1) = qu(:,:,1) ./ qu(:,:,3);
qu(:,:,2) = qu(:,:,2) ./ qu(:,:,3);
% FAIL
max(abs(qu(:,:,1:2)-ulx_(:,:,1:2)))
%}

d = depthT;
wx = pos3d_t1(:,:,1) - pos3d(:,:,1);
wy = pos3d_t1(:,:,2) - pos3d(:,:,2);
wz = pos3d_t1(:,:,3) - pos3d(:,:,3); % really ?? when and where and why

end