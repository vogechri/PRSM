% L1 distance in 2d from the 3d solution (also possible from 2d solutions as well !)
function [N_lin, Rt_lin, Rt_glob] = egoMotion(ref, cam, Seg, flowInit, N_lin, matches)

smartSigma = 0.05;% super critical parameter ?

[M, N, ~] = size(ref.I(1).I);
[u, lambda_uv] = Compute_lambdaUV(ref.R, N, M);
ref.I(1).u = u;

p2d_1 = inv(cam(1).Kl) * reshape(permute (u , [3,1,2]), 3, N*M);

%%%%%%%%%%% first get the 2d warp pixel posistion for all iamges involved!
x.d_ = flowInit(:,:,1);x.w_.x = flowInit(:,:,2);x.w_.y = flowInit(:,:,3);x.w_.z = flowInit(:,:,4);
% for all now derivatives computed, timesteps and cameras
%[ref, cam ] = computeImagePropertiesCamOcc_new( ref, u, lambda_uv, cam, x, 1);
[ref, cam ] = computeImagePropertiesCamOcc_new( ref, u, lambda_uv, cam, x, 0); % stable

test  = (ref.I(2).u(:,:,1:2)-ref.I(1).u(:,:,1:2));
%test2 = (cam.I(2).u(:,:,1:2)-cam.I(1).u(:,:,1:2));
test3 = (cam.I(2).u(:,:,1:2)-ref.I(1).u(:,:,1:2));
%%%%%%%%%%%%%%%%%%%

% flow in left images
u_ = u;u_(:,:,1:2) = u_(:,:,1:2) + test;
q2d_2 = reshape(permute (u_ , [3,1,2]), 3,  N*M);

% flow in right images
u_ = u;u_(:,:,1:2) = u_(:,:,1:2) + test3;
r2d_3 = reshape(permute (u_ , [3,1,2]), 3,  N*M);

% at first we need to map the segmentation in the other image - based on
% the normal stuff, compute H, use H, for the closest pixel inside! :), store position
% and reproject. store pair, position, 2d flow (in a cell for the segment)

Rt_lin = zeros(4,4,numel(Seg.Ids));

% take 1/2/3 pixel from segment far away from outliers use corners or reliable stuff - do not use oob pixel though!
Idsk = Seg.Ids;

globIds = [];

for k=1:numel(Seg.Ids)
  potIds = Idsk{k}+1;
%  potIds ( nonVisibles( potIds ) ) = [];
  
  if numel(potIds) > 15 %~isempty(potIds)
    
    pos = round(Seg.pixCenters(k,:)+1);

    % if locally distinguishable: add it
    if pos(1)>30 && pos(2)>30 && pos(2)<N-30 && pos(1)<M-30
    if sum(sum ( abs( ref.I(1).I(pos(1),pos(2)) - ref.I(1).I(max(1,pos(1)-1):min(pos(1)+1,M), max(1,pos(2)-1):min(pos(2)+1,N))) >1.5/256 )) > 2
      addId = sub2ind(size(Seg.Img), pos(1), pos(2) );
      globIds(end+1) = addId;
    else
      breakHere =1;
    end
    end
  end  
end

iDepth = max( 0.0000001, 1./flowInit(:,:,1) );

% from sparse matches -- not needed
if exist('matches','var') && numel(matches) > 100
  iD = -1./Interpol_mex( flowInit(:,:,1), matches(1,:), matches(2,:) )';
  ps = cam(1).Kl \ cat(1, matches(1,:), matches(2,:), ones(1,size(matches,2)) );
  q1 = cat(1, matches(3,:), matches(4,:) );
  q2 = cat(1, Interpol_mex( reshape( r2d_3(1,:), M, N), matches(1,:), matches(2,:))', Interpol_mex( reshape( r2d_3(2,:), M, N), matches(1,:), matches(2,:) )'); 

  vis1 = true(size(matches,2),1);
  vis2 = false(size(matches,2),1); % TODO no right image corespondences

  iD = cat(2, iD, -iDepth(globIds)); 
  ps = cat(2, ps, p2d_1(:,globIds)); 
  q1 = cat(2, q1, q2d_2(1:2,globIds));
  q2 = cat(2, q2, r2d_3(1:2,globIds));

  vis1 = cat(1, vis1, true(numel(globIds),1));
  vis2 = cat(1, vis2, true(numel(globIds),1));
  
  [RT_new, l1_err] = algebraicL2_Rot_C(cat(1,[0;0;0;0;0;0]), cam.Kl, cam.Kr, cam.Kr*cam.Tr, eye(3), q1, q2, ps, iD, vis1, vis2);
  [RT_new, l1_err] = algebraicL2_Rot_C(cat(1,[0;0;0], -RT_new(1:3,4)), cam.Kl, cam.Kr, cam.Kr*cam.Tr, RT_new(1:3,1:3), q1, q2, ps, iD, vis1, vis2, l1_err,0);
  [RT_new, l1_err] = algebraicL2_Rot_C(cat(1,[0;0;0], -RT_new(1:3,4)), cam.Kl, cam.Kr, cam.Kr*cam.Tr, RT_new(1:3,1:3), q1, q2, ps, iD, vis1, vis2, l1_err,0);
  
else % dense input 
  iD   = -iDepth(globIds); ps = p2d_1(:,globIds); q1=q2d_2(1:3,globIds); q2 = r2d_3(1:3,globIds);  
  
  [RT_new, l1_err] = algebraicL2_Rot_B(cat(1,[0;0;0;0;0;0]), cam.Kl, cam.Rr, cam.Tr, eye(3), q1, q2, ps, iD );
  [RT_new, l1_err] = algebraicL2_Rot_B(cat(1,[0;0;0], -RT_new(1:3,4)), cam.Kl, cam.Rr, cam.Tr, RT_new(1:3,1:3), q1, q2, ps, iD, l1_err,0);
  [RT_new, l1_err] = algebraicL2_Rot_B(cat(1,[0;0;0], -RT_new(1:3,4)), cam.Kl, cam.Rr, cam.Tr, RT_new(1:3,1:3), q1, q2, ps, iD, l1_err,0);
end
RT_new = algebraicL2_Rot  (cat(1,[0;0;0], -RT_new(1:3,4)), cam.Kl, cam.Rr, cam.Tr, RT_new(1:3,1:3), q1, q2, ps, iD, l1_err, smartSigma, 1000);

for k=1:numel(Seg.Ids)
  Rt_lin(:,:,k) = cat(1, RT_new, [0,0,0,1]);
end

Rt_glob = RT_new;

% convert to standard format
centers = findPlaneCenter( Seg, 0, N_lin );
for i = 1:size(centers, 2)
  Rt_lin(1:3,4,i) =  Rt_lin(1:3,4,i) + Rt_lin(1:3,1:3,i) * centers(:,i) - centers(:,i);
end
N_lin = cat(1, N_lin, ones(1,size(N_lin,2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
