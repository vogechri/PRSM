function [occErr, noccErr, epes] = getKittiErr3dSF_flow ( ref, cam, flow_ )

[disp, uDisp, vDisp, dDisp] = convert3Dto2D(ref, cam, flow_(:,:,1), flow_(:,:,2), flow_(:,:,3), flow_(:,:,4));
[occErr, noccErr, epes] = getKittiErrSF ( disp, uDisp, vDisp );

% figure(1), imagesc(disp)
% figure(3), imagesc(uDisp)
% figure(2), imagesc(vDisp)