%%% use the per segment solution to generate a proposal set for the per
%%% pixle refinement. here also proposals from the other views are used iff
%%% those are differnet enegough from the proposal in the canonical view
function [C_projX, newIdsT1 ] = expansionProposals (ref, cam, N_linM, ...
  Rt_linM, centers2DM, camId, Seg, projImg1, centers2D, centers2D2, ...
  vcIdsCam, camCase, Rt_cam )

%camCase == 1: left right cam -- so current time different camera
%camCase == 3: from left cam @ t+1 backwards
%camCase == 5: canonial view to the prev frame 
%camCase == 7: canonial view to the next next frme -- not tested since reimplemented 
%camCase == 7: original cose see below
% Rt_cam && camCase == 3: from left cam @ t-1 forwards

newIdsT1 = vcIdsCam;%vcIds{2};

% reproject to find new cneters here: must reproject from last frame - so
Rt_lr = (cat(1, cat(2,cam(camId).Rr,cam(camId).Tr), [0,0,0,1]));

if camCase == 7
  % normal at t+2 for all
  for i = 1:size(Rt_linM, 3)
    %  Rt_linM_inv(:,:,i) =  backPack.Rt_cam * Rt_linM(:,:,i) * Rt_linM(:,:,i); % BEFORE: note that I invert Rt_cam internally in mex
    Rt_linM_inv(:,:,i) =  inv(backPack.Rt_cam) * Rt_linM(:,:,i) * Rt_linM(:,:,i);
  end
end
if camCase ==5
  for i = 1:size(Rt_linM, 3)
    Rt_linM_inv(:,:,i) =  inv(Rt_linM(:,:,i)) * Rt_cam; % from t to t-1
  end
end
if camCase ==1 % case ==1
  for i = 1:size(Rt_linM, 3)
    Rt_linM_inv(:,:,i) = Rt_lr; % from left at t to right at t
  end
end
if camCase == 3 % % case 3
  Rt_linM_inv = Rt_linM;
end

[N_back, ~, ~] = ProjectMvps(Seg.Edges, cam(camId).Kl, N_linM(1:3,:), Rt_linM_inv, centers2DM(:, 1:size(N_linM,2)), centers2DM(:, 1:size(N_linM,2)), zeros(size(Seg.Img)) );
%[N_back, ~, ~] = ProjProposals6(Seg.Edges, cam(camId).Kl, N_linM(1:3,:), Rt_linM_inv, centers2DM(:, 1:size(N_linM,2)), centers2DM(:, 1:size(N_linM,2)), zeros(size(Seg.Img)) );
% was in camCase == 7 - but why ?
%[N_back, ~, ~] = ProjProposals6(Seg.Edges, cam(camId).Kl, N_linM(1:3,:), Rt_linM_inv, centers2DM, centers2DM, zeros(size(Seg.Img)) );

if camCase ~=3 && camCase ~=7
  if ~exist('Rt_cam','var')
    Rt_lr = inv(Rt_lr);
    for i = 1:size(Rt_linM, 3)
      Rt_linM_inv(:,:,i) =  Rt_lr; % form l to r - but no normals yets
    end
  else
    if camCase ==5
      %% frame t-1
      for i = 1:size(Rt_linM, 3)
        Rt_linM_inv(:,:,i) =  inv(Rt_cam) * Rt_linM(:,:,i); % form t-1 to t - but no normals yets
%        Rt_linM_inv(:,:,i) =  inv(Rt_linM(:,:,i)) * Rt_cam; % from t to t-1
      end
    else
      
    end
  end
else
  for i = 1:size(Rt_linM, 3)
    Rt_linM_inv(:,:,i) =  inv(Rt_linM(:,:,i)); % form t+1 to t
  end
end

[~, ~, C_projX] = ProjectMvps(Seg.Edges, cam(camId).Kl, N_back(1:3,newIdsT1+1), Rt_linM_inv(:,:,newIdsT1+1), centers2D2, centers2D, Seg.Img );
%[~, ~, C_projX] = ProjProposals6(Seg.Edges, cam(camId).Kl, N_back(1:3,newIdsT1+1), Rt_linM_inv(:,:,newIdsT1+1), centers2D2, centers2D, Seg.Img );

% new remove redundant stuff:
for i=1:size(C_projX,2)
  pix = [round(C_projX(2,i)+0.5), round(C_projX(1,i)+0.5)];
  if 1<pix(1) && 1<pix(2) && size(projImg1,1)>pix(1) && size(projImg1,2)>pix(2)
    if projImg1(pix(1), pix(2)) == newIdsT1(i)
      newIdsT1(i) = -1; % invalid
      C_projX(:,i)=-10000;% invalid
    end
  end
end

C_projX(C_projX==-10000) =[];
C_projX=reshape(C_projX, 2, numel(C_projX)/2);
newIdsT1(newIdsT1== -1) = [];

C_projX = inv(cam.Kl) * cat( 1, C_projX, ones(1,size(C_projX,2)) );
% TODO right image missing - but too many then ??
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% well this is how it was done before:
%{
% now for t+2:
%%%%%%%%%%%%%%%%%%%%%%%%%
% normal at t+2 for all
for i = 1:size(Rt_linM, 3)
%  Rt_linM_inv(:,:,i) =  backPack.Rt_cam * Rt_linM(:,:,i) * Rt_linM(:,:,i); % BEFORE: note that I invert Rt_cam internally in mex
  Rt_linM_inv(:,:,i) =  inv(backPack.Rt_cam) * Rt_linM(:,:,i) * Rt_linM(:,:,i);
end
[N_back, ~, ~] = ProjectMvps(Seg.Edges, cam(camId).Kl, N_linM(1:3,:), Rt_linM_inv, centers2DM, centers2DM, zeros(size(Seg.Img)) );
%[N_back, ~, ~] = ProjProposals6(Seg.Edges, cam(camId).Kl, N_linM(1:3,:), Rt_linM_inv, centers2DM, centers2DM, zeros(size(Seg.Img)) );
% now project back, keep centers and test for uniqueness:
for i = 1:size(Rt_linM, 3)
   Rt_linM_inv(:,:,i) =  inv(Rt_linM_inv(:,:,i)); % form t+2 to t
end

[~, ~, C_projX] = ProjectMvps(Seg.Edges, cam(camId).Kl, N_back(1:3,newIds7+1), Rt_linM_inv(:,:,newIds7+1), centers2D2, centers2D, Seg.Img );
%[~, ~, C_projX] = ProjProposals6(Seg.Edges, cam(camId).Kl, N_back(1:3,newIds7+1), Rt_linM_inv(:,:,newIds7+1), centers2D2, centers2D, Seg.Img );

newIdsT1 = newIds7;
for i=1:size(C_projX,2)
  pix = [round(C_projX(2,i)+0.5), round(C_projX(1,i)+0.5)];
  if 1<pix(1) && 1<pix(2) && size(projImg1,1)>pix(1) && size(projImg1,2)>pix(2)
    if projImg1(pix(1), pix(2)) == newIdsT1(i)
      newIdsT1(i) = -1; % invalid
      C_projX(:,i)=-10000;% invalid
    end
  end
end

C_projX(C_projX==-10000) =[];
C_projX=reshape(C_projX, 2, numel(C_projX)/2);
newIdsT1(newIdsT1== -1) = [];

C_projX = inv(cam.Kl) * cat( 1, C_projX, ones(1,size(C_projX,2)) );
centers2DM_pic = cat(2, centers2DM_pic, C_projX);
centerIds = int32(cat(1, centerIds, newIdsT1));
%}