%%% PWRS from ICCV'13 publication, in a minimal setup, used for reduction
%%% of proposals set to speed VC but also to exploit special algorithms
%%% only implemented for this version -- NOTE: these would also work for VC
%%% The reason whhy this is missing for VC is just limited time
%%% especially local replacement (Jit) helps here a lot
function [ oobOracle, initIds, N_lin, Rt_lin ] = ...
  pwrs_init ( ref, cam, Seg, N_prop, RT_prop, ew, dt, pwrs_ds, ts, ...
              dj, tj, dD, oob, maxMot, par, oracle )

  initIds =0;
  [M,N] = size(Seg.Img);
  u  = ones(M,N,3);u(:,:,1) = repmat( [1:N],  M, 1 );u(:,:,2) = repmat( [1:M]', 1, N );
  p2d_ = permute( reshape ( inv(cam(1).Kl) * reshape(permute (u , [3,1,2]), 3, N*M), [3,M,N]), [2,3,1]);
  centers2D = cat(2, Seg.Centers, ones (size(Seg.Centers,1),1))';

  % 1. helper to build a prior about outside pixel
  flowR_2d = oracle.flowR;
  warp_fr(:,:,1) = Interpol_mex( flowR_2d(:,:,1), oracle.stereo(:,:,1)+u(:,:,1), u(:,:,2) );
  warp_fr(:,:,2) = Interpol_mex( flowR_2d(:,:,2), oracle.stereo(:,:,1)+u(:,:,1), u(:,:,2) );
  warp_fr(:,:,1) = oracle.stereo(:,:,1)+warp_fr(:,:,1);
  oobOracle.oobs = oracle.stereo(:,:,1)+u(:,:,1)<1;
  oobOracle.oobsFlow = oracle.flowL(:,:,1)+u(:,:,1)<1 | oracle.flowL(:,:,2)+u(:,:,2)<1 | oracle.flowL(:,:,1)+u(:,:,1)>N | oracle.flowL(:,:,2)+u(:,:,2)>M;
  oobOracle.oobsFlow2 = warp_fr(:,:,1)+u(:,:,1)<1 | warp_fr(:,:,2)+u(:,:,2)<1 | warp_fr(:,:,1)+u(:,:,1)>N | warp_fr(:,:,2)+u(:,:,2)>M;

if par.doSeg ==1  
  newIds = QPBO_Seg(ref.I(1).I, cam(1).I(1).I, ref.I(2).I, cam(1).I(2).I, permute(p2d_, [3,1,2]), cam(1).Kr, Seg.Ids, N_prop(1:3,:), dt, pwrs_ds, ts, Seg.Edges, Seg.Weights, cam(1).Tr, RT_prop, centers2D, cam(1).Rr, dj, tj, dD, oob, Seg.Img, int32(oobOracle.oobs), int32(oobOracle.oobsFlow), int32(oobOracle.oobsFlow2), maxMot, ones(size(ew.x(:,:,1)')), ones(size(ew.y(:,:,1)')), ones(size(ew.xy(:,:,1)')), ones(size(ew.yx(:,:,1)')), Seg.Img'); 

  initIds = 0:size(newIds)-1;  
  nn  = N_prop(:,newIds+1);
  rr  = RT_prop(:,:,newIds+1);
end

if par.doJit ==1
  [newIds, nnj, rrj] = QPBO_Jitter(ref.I(1).I, cam(1).I(1).I, ref.I(2).I, cam(1).I(2).I, permute(p2d_, [3,1,2]), cam(1).Kr, Seg.Ids, nn(1:3,:), dt, pwrs_ds, ts, Seg.Edges, Seg.Weights, cam(1).Tr, rr, centers2D, cam(1).Rr, dj, tj, dD, oob, Seg.Img, int32(oobOracle.oobs), int32(oobOracle.oobsFlow), int32(oobOracle.oobsFlow2), maxMot, ones(size(ew.x(:,:,1)')), ones(size(ew.y(:,:,1)')), ones(size(ew.xy(:,:,1)')), ones(size(ew.yx(:,:,1)')), Seg.Img', cam(1).Kr, 160);
  
  nn=nnj(:,newIds+1); rr=rrj(:,:,newIds+1);
  getKittiErr3dSF ( Seg, ref,cam(1), nn, rr );
end

if par.doEgo ==1
  [nn, rr] = egomotionStep(ref, cam(1), nn, rr, Seg, par, dD, maxMot, p2d_, ...
    centers2D, oobOracle, ew.x(:,:,1), ew.y(:,:,1), ew.xy(:,:,1), ew.yx(:,:,1));
end

  N_res  = nn;
  Rt_res = rr;
  % problem which center?
  N_lin = N_res;Rt_lin = Rt_res;
  centers = findPlaneCenter( Seg, 0, N_lin );
  for i = 1:size(centers, 2)
    Rt_lin(1:3,4,i) =  Rt_lin(1:3,4,i) - Rt_lin(1:3,1:3,i) * centers(:,i) + centers(:,i);
  end