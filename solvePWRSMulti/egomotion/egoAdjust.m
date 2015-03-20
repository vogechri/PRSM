%%% adjust the in between frame camera motion -- compare the chapter
%%% View-consistent multi-frame extension and equaitons therein
%%% in fact it is better to first project the past solution into the
%%% current frame and then (because moving objects are now aligned)
%%% compute the rigid motion to map both proposal sets onto each other this
%%% is done robustly with gradient descent on the lorentzian reprojection
%%% error, because we already heve a good inital solution: the identity
function [Rt_glob] = egoAdjust(ref, cam, Seg, N_t1,R_t1, N_t0,R_t0, Seg_t0, weighting)

if ~exist('weighting','var')
  weighting = 0;
end

smartSigma = 0.05;% super critical parameter
[M, N, ~] = size(ref.I(1).I);
u = Compute_lambdaUV(ref.R, N, M);
ref.I(1).u = u;

% steps:
% 1. forward project t0: normals and rotations
% 2. N_t1 and R_t1 define 2d flow or goal positions
% N_t0 and R_t0 forward define depth (at t1+1) and 2d position
% goal find R|t to bring depth+posisitons towards goal posistions

p2d_ = permute( reshape ( inv(cam(1).Kl) * reshape(permute (u , [3,1,2]), 3, N*M), [3,M,N]), [2,3,1]);
flow2d_t1 =  reconstruc2dFlowHom( ref, cam, N_t1, R_t1, Seg, 0 );
% flow in left images
u_ = u;u_(:,:,1:2) = u_(:,:,1:2) + flow2d_t1(:,:,2:3);
q2d_2 = reshape(permute (u_ , [3,1,2]), 3,  N*M);
% flow in right images
u_ = u;u_(:,:,1:2) = u_(:,:,1:2) + flow2d_t1(:,:,2:3);
u_(:,:,1) = u_(:,:,1) + flow2d_t1(:,:,4); % or minus ?
r2d_3 = reshape(permute (u_ , [3,1,2]), 3,  N*M);

%%% position and depth at the 2nd frame needed, position: ok
%%% depth d+wz?!
centers2D     = cat(2, Seg.Centers, ones (size(Seg.Centers,1),1))';
centers2D_old = cat(2, Seg_t0.Centers, ones (size(Seg_t0.Centers,1),1))';
[N_t0_FW, R_t0_FW] = ProjectCoverMvps(Seg_t0.Edges, cam(1).Kl, N_t0(1:3,:), R_t0, centers2D_old, centers2D, Seg.Img );
flow2d_t0 =  reconstruc2dFlowHom( ref, cam, N_t0_FW, R_t0_FW, Seg, 0 );
u_ = u;u_(:,:,1:2) = u_(:,:,1:2) + flow2d_t0(:,:,2:3);
p2d_1 = inv(cam(1).Kl) * reshape(permute (u_ , [3,1,2]), 3,  N*M);
flow3d_t0 = reconstruc3DFlowHom   ( N_t0_FW, R_t0_FW, Seg, p2d_, 0 ); % indifferent or ..
depth = flow3d_t0(:,:,1)+flow3d_t0(:,:,4);

% at first we need to map the segmentation in the other image - based on
% the normal stuff, compute H, use H, for the closest pixel inside! :), store position
% and reproject. store pair, position, 2d flow (in a cell for the segment)

% take 1/2/3 pixel from segment far away from outliers use corners or reliable stuff - do not use oob pixel though!
Idsk = Seg.Ids;
globIds = [];

for k=1:numel(Seg.Ids)
  potIds = Idsk{k}+1;
  %  potIds ( nonVisibles( potIds ) ) = [];
  
  if numel(potIds) > 15 %~isempty(potIds)
    
    pos = round(Seg.pixCenters(k,:)+1);
    % if locally distinguishable: add it
    border = 10;
    if pos(1)>border && pos(2)>border&& pos(2)<N-border&& pos(1)<M-border
      if sum(sum ( abs( ref.I(1).I(pos(1),pos(2)) - ref.I(1).I(max(1,pos(1)-1):min(pos(1)+1,M), max(1,pos(2)-1):min(pos(2)+1,N))) >1.0/256 )) > 2
        addId = sub2ind(size(Seg.Img), pos(1), pos(2) );
        globIds(end+1) = addId;
      else
        breakHere =1;
      end
    end
  end
end

iDepth = max( 0.0000001, 1./depth(:,:,1) );
iD   = -iDepth(globIds); ps = p2d_1(:,globIds); q1=q2d_2(1:3,globIds); q2 = r2d_3(1:3,globIds);

doFull =0;
if doFull
  [RT_new, l1_err] = algebraicL2_Rot_B(cat(1,[0;0;0;0;0;0]), cam.Kl, cam.Rr, cam.Tr, eye(3), q1, q2, ps, iD );
  [RT_new, l1_err] = algebraicL2_Rot_B(cat(1,[0;0;0], -RT_new(1:3,4)), cam.Kl, cam.Rr, cam.Tr, RT_new(1:3,1:3), q1, q2, ps, iD, l1_err,0);
  [RT_new, l1_err] = algebraicL2_Rot_B(cat(1,[0;0;0], -RT_new(1:3,4)), cam.Kl, cam.Rr, cam.Tr, RT_new(1:3,1:3), q1, q2, ps, iD, l1_err,0);
  
  RT_new = algebraicL2_Rot  (cat(1,[0;0;0], -RT_new(1:3,4)), cam.Kl, cam.Rr, cam.Tr, RT_new(1:3,1:3), q1, q2, ps, iD, l1_err, smartSigma, 1000);
  % maybe sufficient ?
else

  if weighting % not usefull here
    % use depth to weight similarity - rather use inverse :)
    weight_sim = max(1, abs(flow2d_t0(:,:,1) - flow2d_t1(:,:,1)));
    l1_errTest = repmat(1./weight_sim(globIds).^2, 2,1);

    %  pct = prctile(l1_errTest(1:2:end), 90);
    %  l1_errTest (l1_errTest > pct ) = 0.000000000000000001;
    %  plot( 1:numel(l1_errTest(1:2:end)), l1_errTest(1:2:end))
    %  [RT_new, l1_errTest] = algebraicL2_Rot_B(cat(1,[0;0;0;0;0;0]), cam.Kl, cam.Rr, cam.Tr, eye(3), q1, q2, ps, iD );
    
    RT_new = algebraicL2_Rot  ([0;0;0;0;0;0], cam.Kl, cam.Rr, cam.Tr, eye(3), q1, q2, ps, iD, l1_errTest(:), smartSigma, 2000);
  else
    
    [x, fX_] = minimize([0;0;0;0;0;0], @fitHomoR, [2000 0.05], cam.Kl, cam.Rr, cam.Tr, eye(3), q1(1:2,:), ps, iD, 0, smartSigma, q2(1:2,:) );
    
    Rt = cat(1, x(4:6), x(1:3));
    R_cross = [0, -Rt(6), Rt(5); Rt(6), 0, -Rt(4); -Rt(5), Rt(4), 0];
    
    r = Rt(4:6);
    sinA = norm(r);
    alpha = asin(sinA);
    cosA = cos(alpha);
    
    R1 = zeros(3);
    R1(1,1:3:end) = cosA;
    R1(2,2:3:end) = cosA;
    R1(3,3:3:end) = cosA;
    
    R2 = R_cross;
    R3=(1-cosA) * (r*r')./max( eps, sinA^2);
    Rot = R1 + R2 + R3;
    RT_new = [Rot,Rt(1:3)];
    
    visualizeErr =0;
    if visualizeErr % visualization:
      [~, l1_err] = algebraicL2_Rot_B(cat(1,[0;0;0], -x(4:6)), cam.Kl, cam.Rr, cam.Tr, Rot, q1, q2, ps, iD, 0, 1,0);prctile(l1_err, 50)
      errpix = -ones(size(ref.I(1).I));
      %errpix( globIds ) =  l1_err(1:2:end) + l1_err(2:2:end);
      % better vis with
      for i=-2:2
        for j=-2:2
          errpix( globIds+i+j*M ) =  (l1_err(1:2:end) + l1_err(2:2:end)).^0.25;
        end
      end
      figure(8), imagesc(errpix), colorbar;
    end
  end
end

Rt_glob = RT_new;
cat(1, Rt_glob, [0,0,0,0] );
Rt_glob(4,4)=1;
