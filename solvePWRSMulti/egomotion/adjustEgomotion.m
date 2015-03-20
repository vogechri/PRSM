%%% handles the computation of the adjusted proposals from different time
%%% steps than the current one
function  [Rt_cam, N_proj, Rt_proj, C_proj, N_old_lin, Rt_old_lin, ...
           S_old_img, centers2D_old, ids_prevFrame, Seg] = ...
  adjustEgomotion (Seg, ref, cam, par, plotXtraPrev, N_prop, RT_prop)

% could optionally hand over the label-sets prev frame and current
roughFolder = '/cluster/scratch_xp/public/vogechri/consistent/';
%roughFolder = 'C:/Users/vogechri/Desktop/work/init/';

if ~par.testing
  prevDir = sprintf('%s/ECCV_Rough3Frame_f%02d/', roughFolder, par.subImg-1);
else
  prevDir = sprintf('%s/Eccv_RoughTest_f%02d_newSGM_red1/', roughFolder, par.subImg-1);
end

% -- temporary folder -- check if exists
if exist(par.tempFolder,'dir')  %&& 0
%  fprintf('loading from temporary folder\n');
  cprintf('cyan','loading from temporary folder\n');
  roughFolder = par.tempFolder;
  prevDir = sprintf('%s/', roughFolder );
end

[m,n] = size(Seg.Img);
SegO = Seg;

% old solution
load( sprintf('%sRoughSolution%03d_%02d', prevDir, par.imgNr, par.subImg-1 ) );
% TODO: Use the oldIds information for uniqueness and so
Rt_old=Rt_old(:,:,oldIds+1);N_old=N_old(:,oldIds+1);
S_old_img = S_old.Img;
Seg = SegO;
S_oldOld = S_old;

ref_inv = ref; cam_inv = cam;
ImgI      = cell(1,2*par.frames+2);
if par.subImg>0
  ImgI = cam.Iold;
end

ref_inv.I(2).I = ref_inv.I(1).I;
cam_inv.I(2).I = cam_inv.I(1).I;
ref_inv.I(1).I = ImgI{1};
cam_inv.I(1).I = ImgI{2};

centers2D     = cat(2, Seg.Centers, ones (size(Seg.Centers,1),1))';
centers2D_old = cat(2, S_oldOld.Centers, ones (size(S_oldOld.Centers,1),1))';

% Rt_old (..lin) is already in the origin centered format
Rt_old_lin = Rt_old;
N_old_lin = N_old;

[N_proj, Rt_proj, C_proj] = ProjectMvps(Seg.Edges, cam(1).Kl, N_old_lin(1:3,:), Rt_old_lin, centers2D_old, centers2D, Seg.Img );
[N_t0_lin, Rt_t0_lin]     = ProjectCoverMvps(Seg.Edges, cam(1).Kl, N_old_lin(1:3,:), Rt_old_lin, centers2D_old, centers2D, Seg.Img );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Solution projected no consistency\n'); % here already siginficantly different!
getKittiErr3dSF ( Seg, ref, cam, N_t0_lin, Rt_t0_lin, 0 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotXtraPrev == 1
  plotAnalysis(ref, cam, N_t0_lin, Rt_old_lin, Seg, u, 20, par, sprintf('%03d_projectedPrevSolution', par.imgNr));
end
NRT_uni1 = cat( 2, N_proj',reshape(Rt_proj,16, size(Rt_proj,3))' );
NRT_uni2 = cat( 2, N_t0_lin',reshape(Rt_t0_lin,16, size(Rt_t0_lin,3))' );
[~,ids_prevFrame]=ismember(NRT_uni2, NRT_uni1, 'rows');

% centersProp = findPlaneCenter( Seg, 0, N_prop(1:3,1:numel(Seg.Ids)) );
% for i = 1:size(centersProp, 2)
%   RT_prop(1:3,4,i) =  RT_prop(1:3,4,i) - RT_prop(1:3,1:3,i) * centersProp(:,i) + centersProp(:,i);
% end
[N_prop, RT_prop] = convert_to_nonCentered( Seg, N_prop, RT_prop );

Rt_camTestProp = egoAdjust(ref, cam, Seg, N_prop, RT_prop, N_old_lin, Rt_old_lin, S_oldOld, 0);
Rt_cam = Rt_camTestProp;

% fix cam motion, proposals from last into current representation: 
for i=1:size(Rt_proj,3)
  Rt_proj(:,:,i) = Rt_cam * Rt_proj(:,:,i);
end

% 
fprintf('After fitting previous Frame Solution projected\n');
getKittiErr3dSF ( Seg, ref, cam, N_proj(:,ids_prevFrame), Rt_proj(:,:,ids_prevFrame), 0 ); % works
if plotXtraPrev ==1
  [u, lambda_uv] = Compute_lambdaUV(ref.R, n, m);
  plotAnalysis(ref, cam, N_proj(:,ids_prevFrame), Rt_proj(:,:,ids_prevFrame), Seg, u, 20, par, sprintf('%03d_projectedPrevSolution', par.imgNr));
end

Rt_cam