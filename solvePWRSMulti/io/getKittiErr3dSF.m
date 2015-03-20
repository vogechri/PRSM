%%% kitti error measures 
function [occErr, noccErr, epes] = getKittiErr3dSF ( Seg,ref,cam, N_, RT_, RC )

global doKittiErrors;
if doKittiErrors == 0
 epes.epe = 0;epes.epe_noc = 0;epes.epeD = 0;epes.epe_nocD = 0;
 occErr.err2 = 1;occErr.err3 = 1;occErr.err4 = 2;occErr.err5 = 3;
 occErr.err2f = 1;occErr.err3f = 2;occErr.err4f = 3;occErr.err5f = 4;
 noccErr.err2 = 4;noccErr.err3 = 5;noccErr.err4 = 1;noccErr.err5 = 3;
 noccErr.err2f = 1;noccErr.err3f = 2;noccErr.err4f = 3;noccErr.err5f = 4;
 return;
end

if ~exist('RC', 'var') 
  RC = 1;
end

[m,n] = size(cam.I(1).I);
u  = ones(m,n,3);
u(:,:,1) = repmat( [1:n],  m, 1 );
u(:,:,2) = repmat( [1:m]', 1, n );
p2d_ = permute( reshape ( inv(cam(1).Kl) * reshape(permute (u , [3,1,2]), 3, n*m), [3,m,n]), [2,3,1]);

flow_ =  reconstruc3DFlowHom( N_(1:3,:), RT_, Seg, p2d_, RC );

[occErr, noccErr, epes] = getKittiErr3dSF_flow ( ref, cam, flow_ );
