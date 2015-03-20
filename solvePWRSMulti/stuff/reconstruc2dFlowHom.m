%%% input Normals, Rotations and translation, Sgementation, lines on which
%%% the points lie (K_l ^-1 * pixel)
function flow =  reconstruc2dFlowHom( ref, cam, N_res, Rt_res, Seg, ver )

if ~exist('ver','var')
  ver = 1;% rotation per center
end

[M, N, ~] = size(ref.I(1).I);
u        = ones(M,N,3);
u(:,:,1) = repmat( [1:N],  M, 1 );
u(:,:,2) = repmat( [1:M]', 1, N );

flow = zeros(M,N,4);

p2d_ = permute( reshape ( inv(cam(1).Kl) * reshape(permute (u , [3,1,2]), 3, N*M), [3,M,N]), [2,3,1]);
if ver 
  center = findPlaneCenter( Seg, p2d_, N_res);
else
  center = zeros(size(N_res(1:3,:)));
end
%cam = Homography_cam_pixel       (ref, cam, 1, N_res, Rt_res, Seg, 1 );

cam = Homography_cam_pixel_Center(ref, cam, 1, N_res, Rt_res, Seg, 1, center);
cam = Homography_cam_pixel_Center(ref, cam, 1, N_res, Rt_res, Seg, 2, center);
ref = Homography_ref_pixel_Center(ref, cam, 2, N_res, Rt_res, Seg, center);

flow(:,:,1) = cam.I(1).u(:,:,1) - u(:,:,1);
flow(:,:,2) = ref.I(2).u(:,:,1) - u(:,:,1);
flow(:,:,3) = ref.I(2).u(:,:,2) - u(:,:,2);
flow(:,:,4) = cam.I(2).u(:,:,1) - ref.I(2).u(:,:,1);
