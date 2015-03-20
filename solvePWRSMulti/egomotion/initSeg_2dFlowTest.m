%%% generate Proposals from the inputs stereo2d_t0, stereo2d_t1, flow_left, flow_right
%%% stereo2d_t1 is actually not used. 
%%% first fits solution with a simple iterated (2 its) least squares (L1 distance)
%%% minimizing the reprojection error in 2d. then refinement with
%%% lorentzian cost
function [N_lin, Rt_lin] = initSeg_2dFlowTest(ref, cam, Seg, stereo2d_t0, stereo2d_t1, flow_left, flow_right, use2Views, maxSegs)

if ~exist('use2Views','var')
  use2Views =1;
end

% todo fix nonlinear steps for rotations - here selected stuff missing?
% just joint ids and descent ?

planeReps = 4;% for plane fit algebraic - 1 enough? no time save anyway -> save 1s for rotations as well
patchCover = 1; % if on reduces proposals nicely

depthFromPlanes = 1;% good for VC ?
algebraicOnly = 1;% off: 8.5 on : 20s
if ~algebraicOnly
  RT_its = 10;
  N_itsa = 20;
  N_itsb = 20;
else
  RT_its = 2;% lorentzian steps default 2
  N_itsa = 4;% lorentzian steps default 4
  N_itsb = 0;
end

% con not bo 0 for technical reasons
stereo2d_t0 = min(-0.000001, stereo2d_t0);
stereo2d_t1 = min(-0.000001, stereo2d_t1);

[M, N, ~] = size(ref.I(1).I);
u(:,:,1) = repmat( [1:N],  M, 1 );
u(:,:,2) = repmat( [1:M]', 1, N );
u(:,:,3) = ones(M,N);

%stereo, p2d_1 is K^-1 * p
p2d_1 = inv(cam(1).Kl) * reshape(permute (u , [3,1,2]), 3, N*M);
u_ = u;u_(:,:,1) = u_(:,:,1) + stereo2d_t0(:,:,1);

% this is the desired location in pixel
p2d_2 = reshape(permute (u_ , [3,1,2]), 3,  N*M);
p2d_ = permute( reshape ( inv(cam(1).Kl) * reshape(permute (u , [3,1,2]), 3, N*M), [3,M,N]), [2,3,1]);

% flow in left images
u_ = u;u_(:,:,1:2) = u_(:,:,1:2) + flow_left(:,:,1:2);
q2d_2 = reshape(permute (u_ , [3,1,2]), 3,  N*M);
oobsFlow = u_(:,:,1)<1 | u_(:,:,2)<1 |u_(:,:,1)>N | u_(:,:,2)>M;

if use2Views
test3(:,:,1) = Interpol_mex( flow_right(:,:,1), stereo2d_t0(:,:,1)+u(:,:,1), u(:,:,2) );
test3(:,:,2) = Interpol_mex( flow_right(:,:,2), stereo2d_t0(:,:,1)+u(:,:,1), u(:,:,2) );

% HOGUET: YES!! - WHY ??
%test3(:,:,1) = Interpol_mex( flow_right(:,:,1), u(:,:,1), u(:,:,2) );
%test3(:,:,2) = Interpol_mex( flow_right(:,:,2), u(:,:,1), u(:,:,2) );

u_(:,:,1) = u(:,:,1) + stereo2d_t0(:,:,1)+test3(:,:,1);
u_(:,:,2) = u(:,:,2) + test3(:,:,2);
% s2d_2 the pixel position of the flow in the right image to be fitted at
s2d_2 = reshape(permute (u_ , [3,1,2]), 3,  N*M);
% if oob do not use this information
oobsFlow2 = u_(:,:,1)<1 | u_(:,:,2)<1 |u_(:,:,1)>N | u_(:,:,2)>M;
end

N_lin  = zeros(3,numel(Seg.Ids));

validPixel = false(M,N);
validPixel(2:2:end,2:2:end) = true;
%validPixel = true(M,N);
Idsk = Seg.Ids;

%%%% new
if patchCover

  depth = convert2Dto3D(ref, cam, stereo2d_t0(:,:,1), flow_left(:,:,1), flow_left(:,:,2), stereo2d_t1(:,:,1)-stereo2d_t0(:,:,1));
  iDepth = max( 0.0000001, 1./depth(:,:,1) );

  if use2Views
    [IdsJoint, newIdsk, iSegIdsK, newIdsk2, iSegIdsK2, validMap1, validMap2] = ...
      algebraicMatrixQ9LR( cam, eye(3), p2d_2, q2d_2, s2d_2, p2d_1, -iDepth, Seg.Ids, Seg, ~oobsFlow, ~oobsFlow2, maxSegs );
  else
    [newIdsk, iSegIdsK] = algebraicMatrixQ9( cam, eye(3), q2d_2, p2d_2, p2d_1, -iDepth, Seg.Ids, Seg, true(M,N), maxSegs );
    IdsJoint = newIdsk;
    newIdsk2 = {};iSegIdsK2 = {};
%    [IdsJoint, newIdsk, iSegIdsK, newIdsk2, iSegIdsK2, validMap1, validMap2] = ...
%      algebraicMatrixQ9LR( cam, eye(3), p2d_2, q2d_2, [], p2d_1, -iDepth, Seg.Ids, Seg, ~oobsFlow, false(size(oobsFlow)), maxSegs );    
  end

  fprintf( 'Reduced to %d normals\n', numel(newIdsk));
  Idsk  = IdsJoint';
end
%%%%

sigmas= zeros(1,numel(Idsk));
for k=1:numel(Idsk)
  
%%%% new  
if patchCover
  id4s = Idsk{k}+1;
else
  % normal case:
  id4s = Idsk{k}(validPixel(Idsk{k}+1))+1; 
end
%%%
  [N_new, ~, sig]=  algebraicNormal_l1_rep(cam, p2d_1(:,id4s), p2d_2(:,id4s), planeReps );
  sigmas(k) = sig;
  N_lin(:,k) = N_new;
end

N_linSave = N_lin;
[N_lin, fX_] = minimizeVektor(N_lin(:,1:numel(Idsk)), @fitHomoAllN, [N_itsa 0.1], Idsk, cam.Kr*cam.Tr, p2d_1, cam.Kr*cam.Rr*p2d_1, p2d_2(1:2,:), 0.5*sigmas,  int32(validPixel));
[N_lin, fX_] = minimizeVektor(N_lin(:,1:numel(Idsk)), @fitHomoAllN, [N_itsb 0.1], Idsk, cam.Kr*cam.Tr, p2d_1, cam.Kr*cam.Rr*p2d_1, p2d_2(1:2,:), 0.05*sigmas, int32(validPixel));
N_lin(isnan(N_lin)) = N_linSave(isnan(N_lin));

%%% new:
if patchCover
N_linS = N_lin;
for i=1:numel(iSegIdsK)
  for k=1:numel(iSegIdsK{i})
    temp = iSegIdsK{i};
    N_lin(:,temp(k)) = N_linS(:,i);
  end
end
%N_lin2 = N_lin;% alloc memory
for i=1:numel(iSegIdsK2)
  for k=1:numel(iSegIdsK2{i})
    temp = iSegIdsK2{i};
    N_lin(:,temp(k)+numel(Seg.Ids)) = N_linS(:,i);
  end
end
end
%%%
if patchCover
fprintf( 'Reduced to %d rigid motions\n', numel(newIdsk));
Idsk  = newIdsk';
Idsk2 = newIdsk2';
end

Rt_start = zeros(3,3,numel(Seg.Ids));
r_start = zeros(6,numel(Seg.Ids));
sigmas   = zeros(1,numel(Seg.Ids));
for k=1:numel(Idsk)

  if use2Views
  idjoint  = cat(1,Idsk{k}+1, Idsk2{k}+1 );
%  idjoint1 = cat(1, true(size(Idsk{k}+1)), false(size(Idsk2{k}+1)) );
%  idjoint2 = cat(1, false(size(Idsk{k}+1)), true(size(Idsk2{k}+1)) );
  idjoint1 = cat(1, logical( validMap1{k}), ~logical( validMap2{k}) );
  idjoint2 = cat(1, ~logical( validMap1{k}), logical( validMap2{k}) );
  id4s = idjoint;
  if depthFromPlanes
    normal = N_linS(1:3,k);
    points = cat( 2, p2d_(idjoint), p2d_(idjoint + N*M) , p2d_(idjoint + 2*N*M))';
    dd = -(1./(points'*normal));
    iDepth(id4s) = max( 0.0000001, 1./dd(:,:,1) );
  end
  if any(idjoint2)
    [RT_new, l1_err] = algebraicL2_Rot_C([0,0,0,0,0,0]', cam.Kl, cam.Kr, cam.Kr*cam.Tr, eye(3), q2d_2(:,id4s), s2d_2(:,id4s), p2d_1(:,id4s), -iDepth(id4s)', idjoint1, idjoint2);
    %      [RT_new, l1_err] = algebraicL2_Rot_C([0,0,0,0,0,0]', cam.Kl, cam.Kr, cam.Kr*cam.Tr, RT_new(1:3,1:3), q2d_2(:,id4s), s2d_2(:,id4s), p2d_1(:,id4s), -iDepth(id4s)', idjoint1, idjoint2, l1_err);
    else % just 1st view:
      [RT_new, l1_err] = algebraicL2_Rot_B([0,0,0,0,0,0]', cam.Kl, cam.Rr, cam.Tr, eye(3), q2d_2(:,id4s), [], p2d_1(:,id4s), -iDepth(id4s)');
      %      [RT_new, l1_err] = algebraicL2_Rot_B([0,0,0,0,0,0]', cam.Kl, cam.Rr, cam.Tr, RT_new(1:3,1:3), q2d_2(:,id4s), [], p2d_1(:,id4s), -iDepth(id4s)', l1_err);
    end
  else % ue only single view
  id4s = Idsk{k}+1;
  if depthFromPlanes
    normal = N_linS(1:3,k);
    points = cat( 2, p2d_(id4s), p2d_(id4s + N*M) , p2d_(id4s + 2*N*M))';
    dd = -(1./(points'*normal));
    iDepth(id4s) = max( 0.0000001, 1./dd(:,:,1) );
  end
  [RT_new, l1_err] = algebraicL2_Rot_B([0,0,0,0,0,0]', cam.Kl, cam.Rr, cam.Tr, eye(3), q2d_2(:,id4s), [], p2d_1(:,id4s), -iDepth(id4s)');
%    [RT_new, l1_err] = algebraicL2_Rot_B([0;0;0;RT_new(1:3,4)], cam.Kl, cam.Rr, cam.Tr, RT_new(1:3,1:3), q2d_2(:,id4s), [], p2d_1(:,id4s), -iDepth(id4s)' , l1_err);
  end
  sigmas(k) = sum(sqrt(l1_err(1:2:end).^2 + l1_err(2:2:end).^2))/numel(l1_err);
  Rt_start(:,:,k) = RT_new(1:3,1:3);
  r_start(4:6,k)  = RT_new(1:3,4);
end

if use2Views
  % lacks some smartness i.e. booleans to tell where parts are valid from q and s
  Idsk  = IdsJoint';
  Rt_lin = fullRotationFit(r_start(:,1:numel(Idsk)), Idsk, cam, Rt_start, q2d_2, p2d_1, -iDepth, 0.5*sigmas, RT_its, s2d_2);
else
  Rt_lin = fullRotationFit(r_start(:,1:numel(Idsk)), Idsk, cam, Rt_start, q2d_2, p2d_1, -iDepth, 0.5*sigmas, RT_its);
end

% new:
if patchCover
Rt_linS = Rt_lin;
for i=1:numel(iSegIdsK)
  for k=1:numel(iSegIdsK{i})
    temp = iSegIdsK{i};
    Rt_lin(:,:,temp(k)) = Rt_linS(:,:,i);
  end
end
%Rt_lin2 = Rt_lin;
for i=1:numel(iSegIdsK2)
  for k=1:numel(iSegIdsK2{i})
    temp = iSegIdsK2{i};
    Rt_lin(:,:,temp(k)+numel(Seg.Ids)) = Rt_linS(:,:,i);
  end
end
end

%Rt_lin(1:3,1:3,:) = Rt_start(1:3,1:3,:);
%Rt_lin(1:3,4,:)   = r_start(4:6,:);

Rt_R = Rt_lin(1:3,1:3,:);
Rt_T = Rt_lin(1:3,4,:);

T_linSave = reshape(r_start(4:6,:), [3,1, size(r_start,2)]);

Rt_R(isnan(Rt_R)) = Rt_start(isnan(Rt_R));
Rt_T(isnan(Rt_T)) = T_linSave(isnan(Rt_T));

Rt_lin(1:3,1:3,:) = Rt_R;
Rt_lin(1:3,4,:)   = Rt_T;

% convert to standard format
centers = findPlaneCenter( Seg, 0, N_lin(:,1:numel(Seg.Ids)) );
if use2Views 
  centers = cat(2, centers, findPlaneCenter( Seg, 0, N_lin(:,numel(Seg.Ids)+1:end ) ));
end

for i = 1:size(centers, 2)
  Rt_lin(1:3,4,i) =  Rt_lin(1:3,4,i) + Rt_lin(1:3,1:3,i) * centers(:,i) - centers(:,i);
end
N_lin = cat(1, N_lin, ones(1,size(N_lin,2)));

%remove nan's! from the start - there should be none - but who knows
bad_id = ceil(find(isnan(N_lin))/4);
bad_id = bad_id(1:3:end);
for i=1:numel(bad_id)
  N_lin( 1:3, bad_id ) = N_lin( 1:3, bad_id-1 );
end
bad_id = ceil(find(isnan(Rt_lin))/16);
bad_id = bad_id(1:3:end);
for i=1:numel(bad_id)
  Rt_lin( 1:3, 4, bad_id ) =  Rt_lin( 1:3, 4, bad_id-1 );
end
