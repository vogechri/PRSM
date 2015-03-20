%%% apply egomotion proposals for ICCV13 version
function [N_res, Rt_res] = egomotionStep(ref, cam, N_res, Rt_res, Seg, p, ...
  dD, maxMot, p2d_, centers2D, oracle, weight_x, weight_y, weight_xy, weight_ixy)

dt   = p.dt;
ts   = p.ts;
ds   = 0.03;%p.ds*0.5;
dj   = p.dj;
tj   = p.tj;
oob  = p.oob;

flow_4 =  reconstruc3DFlowHom( N_res, Rt_res, Seg, p2d_, 1 );
[N_lin, Rt_lin, Rt_glob] = egoMotion(ref, cam, Seg, flow_4, N_res(1:3,:));
fprintf('Egomotion L2\n');

[newIds, nn, rr] = QPBO_Jitter(ref.I(1).I, cam(1).I(1).I, ref.I(2).I, cam(1).I(2).I, permute(p2d_, [3,1,2]), cam(1).Kl, Seg.Ids, N_lin(1:3,:), dt, ds, ts, Seg.Edges, Seg.Weights, cam(1).Tr, Rt_lin, centers2D, cam(1).Rr, dj, tj, dD, oob, Seg.Img, int32(oracle.oobs), int32(oracle.oobsFlow), int32(oracle.oobsFlow2), maxMot, ones(size(weight_x')),  ones(size(weight_y')),  ones(size(weight_xy')), ones(size(weight_ixy')), Seg.Img', cam(1).Kr, 10);

N_ego = nn(:,newIds+1); Rt_ego = rr(:,:,newIds+1);
getKittiErr3dSF ( Seg,ref,cam, N_ego, Rt_ego );

fprintf('Egomotion Normals\n');
N_scale = egoMotionNormals(ref, cam, Seg, N_lin, Rt_lin, N_res, Rt_res);

[newIds, nn, rr] = QPBO_Jitter(ref.I(1).I, cam(1).I(1).I, ref.I(2).I, cam(1).I(2).I, permute(p2d_, [3,1,2]), cam(1).Kl, Seg.Ids, N_scale(1:3,:), dt, ds, ts, Seg.Edges, Seg.Weights, cam(1).Tr, Rt_lin, centers2D, cam(1).Rr, dj, tj, dD, oob, Seg.Img, int32(oracle.oobs), int32(oracle.oobsFlow), int32(oracle.oobsFlow2), maxMot, ones(size(weight_x')),  ones(size(weight_y')),  ones(size(weight_xy')),  ones(size(weight_ixy')), Seg.Img', cam(1).Kr, 10);

N_scale = nn(:,newIds+1); Rt_scale = rr(:,:,newIds+1);
% does the normal scaling work ?
getKittiErr3dSF ( Seg,ref,cam, N_scale, Rt_scale );

fprintf('Egomotion Normal\n');
N_scaleN = egoMotionNormalsN(ref, cam, Seg, N_lin, N_res, Rt_res, p2d_, Rt_glob);

[newIds, nn, rr] = QPBO_Jitter(ref.I(1).I, cam(1).I(1).I, ref.I(2).I, cam(1).I(2).I, permute(p2d_, [3,1,2]), cam(1).Kl, Seg.Ids, N_scaleN(1:3,:), dt, ds, ts, Seg.Edges, Seg.Weights, cam(1).Tr, Rt_lin, centers2D, cam(1).Rr, dj, tj, dD, oob, Seg.Img, int32(oracle.oobs), int32(oracle.oobsFlow), int32(oracle.oobsFlow2), maxMot, ones(size(weight_x')),  ones(size(weight_y')),  ones(size(weight_xy')),  ones(size(weight_ixy')), Seg.Img', cam(1).Kr, 10);

N_scaleN = nn(:,newIds+1); Rt_scaleN = rr(:,:,newIds+1);
getKittiErr3dSF ( Seg, ref, cam, N_scaleN, Rt_scaleN );

RT_prop= cat(3, Rt_res, Rt_ego, Rt_scale, Rt_scaleN);
N_prop = cat(2,  N_res,  N_ego,  N_scale,  N_scaleN );
% ego onlyOB: better, mostly in occluded areas

%newIds = QPBO_FuseOcc(ref.I(1).I, cam(1).I(1).I, ref.I(2).I, cam(1).I(2).I, permute(p2d_, [3,1,2]), cam(1).Kr, Seg.Ids, N_prop(1:3,:), dt, ds, ts, Seg.Edges, Seg.Weights, cam(1).Tr, RT_prop, centers2D, cam(1).Rr, dj, tj, dD, oob, Seg.Img, int32(oracle.oobs), int32(oracle.oobsFlow), int32(oracle.oobsFlow2), maxMot, ones(size(weight_x')),  ones(size(weight_y')),  ones(size(weight_xy')),  ones(size(weight_ixy')), Seg.Img');
newIds = QPBO_Fuse(ref.I(1).I, cam(1).I(1).I, ref.I(2).I, cam(1).I(2).I, permute(p2d_, [3,1,2]), cam(1).Kr, Seg.Ids, N_prop(1:3,:), dt, ds, ts, Seg.Edges, Seg.Weights, cam(1).Tr, RT_prop, centers2D, cam(1).Rr, dj, tj, dD, oob, Seg.Img, int32(oracle.oobs), int32(oracle.oobsFlow), int32(oracle.oobsFlow2), maxMot, ones(size(weight_x')),  ones(size(weight_y')),  ones(size(weight_xy')),  ones(size(weight_ixy')), Seg.Img');

N_res  = N_prop(:,newIds+1);
Rt_res = RT_prop(:,:,newIds+1);
[oErr, noErr] = getKittiErr3dSF ( Seg,ref,cam, N_res, Rt_res );
%%%%%%%%%%%% EGO MOTION %%%%%%%%%%%%%%%%%%%
