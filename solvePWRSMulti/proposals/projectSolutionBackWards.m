%%% include a frame from the future: 
%%% this was the old 4 frame implementation -- however the code was
%%% adjusted and this part was left for future implementation 
%%% thus it can only be used as a reference for this job
function backPack = projectSolutionBackWards(par, Seg, ref, cam, Rt_glob, Seg2, dir11)

plotresults  = 0;
plotXtraPrev = 0;

%, N_lin, Rt_lin)
%[ids_prevFrame, N_proj, Rt_proj, C_proj, Rt_cam] 

% roughFolder = 'C:/Users/vogechri/Desktop/work/init/';
% roughFolder = '/cluster/scratch_xl/public/vogechri/consistent/';

% dir11 = sprintf('%s/JanTemp_Rough_f%d/', roughFolder, par.subImg+1);

% from outside:
%   load( sprintf('%sRoughSolution%03d_%02d', dir10, par.imgNr, par.subImg ) );
%   fprintf('Init-Solution, no consistency\n');  
%   getKittiErr3dSF ( S_old, ref, cam, N_old, Rt_old, 0 );
%   Seg = S_old;N_res = N_old; Rt_res = Rt_old;
%   p2d_ = permute( reshape ( inv(cam(1).Kl) * reshape(permute (u , [3,1,2]), 3, n*m), [3,m,n]), [2,3,1]);
%   flow4 = reconstruc3DFlowHom ( N_res, Rt_res, Seg, p2d_, 0 );% actualy take rather from solution
%   [~,~,Rt_glob] = egoMotion(ref, cam, Seg, flow4, N_res(1:3, 1:numel(Seg.Ids)));


%    N_lin = N_res;Rt_lin = Rt_res;


  [m,n] = size(Seg.Img);
  [u, lambda_uv] = Compute_lambdaUV(ref.R, n, m);
  p2d_ = permute( reshape ( inv(cam(1).Kl) * reshape(permute (u , [3,1,2]), 3, n*m), [3,m,n]), [2,3,1]);

  % old solution
  %  load( sprintf('%sRoughSolution%03d_%02d', dir11, par.imgNr, par.subImg+2 ) );
  %challenges:
  if isfield( ref, 'Inew' )
    load( sprintf('%sRoughSolution%03d_%02d', dir11, par.imgNr+1, par.subImg ) );
  else
    load( sprintf('%sRoughSolution%03d_%02d', dir11, par.imgNr, par.subImg+1 ) );
  end
  
  % TODO: Use the oldIds information for uniqueness and so
  Rt_old=Rt_old(:,:,oldIds+1);N_old=N_old(:,oldIds+1);
  S_old_img = S_old.Img;
  S_oldOld = S_old;


  ref_inv = ref; cam_inv = cam;
  ImgI      = cell(1,2*par.frames+2);
  strNumberi = sprintf( '%06d_%02d', par.imgNr, par.subImg+2 );

  if isfield( ref, 'Inew' )
    ImgI{1} = ref.Inew ;
    ImgI{2} = cam.Inew ;
  else  
    if par.testing
      %ImgI{1} = double(imread(sprintf('%s/%s/%s.png', '../../data/', '/data_stereo_flow/testing/image_0', strNumberi )))/255;
      %ImgI{2} = double(imread(sprintf('%s/%s/%s.png', '../../data/', '/data_stereo_flow/testing/image_1', strNumberi )))/255;
      ImgI{1} = double(imread(sprintf('%s/%s/%s.png', '../../../data/', 'data_stereo_flow/data_stereo_flow/testing/image_0', strNumberi )))/255;
      ImgI{2} = double(imread(sprintf('%s/%s/%s.png', '../../../data/', 'data_stereo_flow/data_stereo_flow/testing/image_1', strNumberi )))/255;
    else
      %ImgI{1} = double(imread(sprintf('%s/%s/%s.png', '../../data/', '/data_stereo_flow/training/image_0', strNumberi )))/255;
      %ImgI{2} = double(imread(sprintf('%s/%s/%s.png', '../../data/', '/data_stereo_flow/training/image_1', strNumberi )))/255;      
      ImgI{1} = double(imread(sprintf('%s/%s/%s.png', '../../../data/', 'data_stereo_flow/data_stereo_flow/training/image_0', strNumberi )))/255;
      ImgI{2} = double(imread(sprintf('%s/%s/%s.png', '../../../data/', 'data_stereo_flow/data_stereo_flow/training/image_1', strNumberi )))/255;
    end
  end

  ref_inv.I(1).I = ref_inv.I(2).I;
  cam_inv.I(1).I = cam_inv.I(2).I;
  ref_inv.I(2).I = ImgI{1};
  cam_inv.I(2).I = ImgI{2};  
  %
if plotXtraPrev
%  plotAnalysis(ref_inv, cam_inv, N_old, Rt_old, S_old, u, 20, par, sprintf('AA%03d_NextSolution_%02d', par.imgNr, par.subImg+1));
end

  centers2D = cat(2, Seg.Centers, ones (size(Seg.Centers,1),1))';
  centers2D_old = cat(2, S_old.Centers, ones (size(S_old.Centers,1),1))';

  Rt_old_lin = Rt_old;

  % fitted from old works better backwards! but not forwards - it is a
  % chain forwards with new normals, and new positions motion of segi x
  % is now in segi y
  flow4_old = reconstruc3DFlowHom ( N_old, Rt_old, S_old, p2d_, 0 );
  [~,~,Rt_glob_old] = egoMotion(ref, cam, S_old, flow4_old, N_old(1:3,:));

  % this is a nice solution, defined already in view centred coordinates
%  getKittiErr3dSF ( S_old, ref, cam, N_old, repmat(Rt_glob, [1,1,size(Rt_old,3)]), 0 ); % OK

  Rt_cam = cat(1, Rt_glob, [0,0,0,1])*inv(cat(1, Rt_glob_old,[0,0,0,1]));
  %  Rt_cam(:,4) =0;   Rt_cam(4,4) =1; only rotation sucks

%%%%
  Rt_old_inv = Rt_old;
  for i = 1:size(Rt_old, 3)
%    Rt_old_inv(:,:,i) = inv(Rt_old(:,:,i)) * inv(Rt_cam); % new : ok?
    Rt_old_inv(:,:,i) = inv(Rt_old(:,:,i) * Rt_cam); % same 
    %TODO:
%    Rt_old_inv(:,:,i) = inv(Rt_old(:,:,i)) * Rt_cam; % original one: wrong?
  end

% this is for getting N_proj and C_proj, ids_prevFrame and Rt_proj
   %call by index - Seg or S_old or what ?
  [N_proj, Rt_proj, C_proj] = ProjProposals6(Seg.Edges, cam(1).Kl, N_old(1:3,:), Rt_old_inv, centers2D_old, centers2D, Seg.Img );
  [N_t0_lin, Rt_t0_lin]     = timeProposals_noCenter5(Seg.Edges, cam(1).Kl, N_old(1:3,:), Rt_old_inv, centers2D_old, centers2D, Seg.Img );

  fprintf('Solution projected no consistency\n');
  getKittiErr3dSF ( Seg, ref, cam, N_t0_lin, Rt_old, 0 );


%if  plotXtraPrev ==1  
%  plotAnalysis(ref, cam, N_t0_lin, Rt_old_lin, Seg, u, 20, par, sprintf('%03d_projectedPrevSolution', par.imgNr));
%end
  NRT_uni1 = cat( 2, N_proj',reshape(Rt_proj,16, size(Rt_proj,3))' );
  NRT_uni2 = cat( 2, N_t0_lin',reshape(Rt_t0_lin,16, size(Rt_t0_lin,3))' );
%  NRT_uni1 = round (NRT_uni1*10000)/10000;
%  NRT_uni2 = round (NRT_uni2*10000)/10000;
  [~,ids_prevFrame]=ismember(NRT_uni2, NRT_uni1, 'rows');

%    N_lin=N_old(:,oldIds+1), Rt_lin=Rt_old(:,:,oldIds+1)

  for i=1:size(Rt_proj,3)
    Rt_proj(:,:,i) = Rt_cam * Rt_old(:,:,i);
  end
  % equivalent:
%   for i=1:size(Rt_proj,3)
%     Rt_t0_lin(:,:,i) = inv(Rt_t0_lin(:,:,i));
%   end
%   getKittiErr3dSF ( Seg, ref, cam, N_t0_lin, Rt_t0_lin, 0 ); % works    
%%%%%

  fprintf('Next Frame Solution projected into current\n');
  getKittiErr3dSF ( Seg, ref, cam, N_proj(:,ids_prevFrame), Rt_proj(:,:,ids_prevFrame), 0 ); % works

if plotXtraPrev
%  plotAnalysis(ref, cam, N_proj(:,ids_prevFrame), Rt_proj(:,:,ids_prevFrame), Seg, u, 20, par, sprintf('AA%03d_testProjNext_%02d', par.imgNr, par.subImg));
end

if plotresults==1
% plot with projected solution next frame:
  [N_projFW, Rt_projFW]     = timeProposals_noCenter5(S_old.Edges, cam(1).Kl, N_proj(:,ids_prevFrame), Rt_proj(:,:,ids_prevFrame), centers2D, centers2D_old, S_old.Img );
  for i=1:size(Rt_projFW,3)
    Rt_projFW(:,:,i) = inv(Rt_cam) * Rt_projFW(:,:,i);
  end
%  plotAnalysis(ref_inv, cam_inv, N_projFW, Rt_projFW, S_old, u, 20, par, sprintf('AA%03d_NextFWSolution_%02d', par.imgNr, par.subImg+1));
end

%{
  Rt_t0_linX = Rt_t0_lin;
  for i=1:size(Rt_t0_lin,3)
    Rt_t0_linX(:,:,i) = Rt_cam * Rt_t0_lin(:,:,i);
  end
  fprintf('global fit\n');
  getKittiErr3dSF ( Seg, ref, cam, N_t0_lin, Rt_t0_linX, 0 ); % now just different
  
  % fix cam motion
  for i=1:size(Rt_proj,3)
    Rt_proj(:,:,i) = Rt_cam * Rt_proj(:,:,i);
  end
%}

% sanity CHECK
%   for i = 1:size(Rt_proj, 3)
%     Rt_proj_inv(:,:,i) =  inv(Rt_proj(:,:,i)) * Rt_cam;
%   end
%   [N_proj_back, Rt_proj_back, C_proj_back] = ProjProposals6(Seg.Edges, cam(1).Kl, N_proj(1:3,:), Rt_proj_inv, centers2D, centers2D_old, Seg.Img );
% CHECK:  N_proj_back == N_oldOld; % ok  
  
  
%  N_prop  = cat(2, cat(1, N_t0, ones(1,size(N_t0,2))), N_prop(1:4,:));
%  RT_prop = cat(3, Rt_t0_add, RT_prop);

%{
  % from prev frame, backward map - does it look ok?
  for i = 1:size(centers2D, 2)
    Rt_inv(:,:,i) =  inv(Rt_t0_linX(:,:,i)) * Rt_cam;
  end
if  plotXtra ==1  
  plotAnalysis(ref_inv, cam_inv, N_t0_lin, Rt_inv, Seg, u, 20, par, sprintf('%03d_backWarp_init', par.imgNr));
end
%   for i = 1:size(centers2D, 2)
%     Rt_inv(:,:,i) =  Rt_cam * inv(Rt_t0_linX(:,:,i));
%   end
%   plotAnalysis(ref_inv, cam_inv, N_t0_lin, Rt_inv, Seg, u, 20, par, sprintf('%03d_testWarp_lin', par.imgNr));
  
% now test differnet stuff: with pic 110
  for i = 1:size(centers2D, 2)
%    Rt_inv(:,:,i) =  inv(Rt_cam * Rt_lin(:,:,i)); % sucks
%    Rt_inv(:,:,i) =  inv( Rt_lin(:,:,i) ); % sucks
    Rt_inv(:,:,i) =  Rt_cam * inv( Rt_lin(:,:,i) ); % ok
%    Rt_inv(:,:,i) =  inv(Rt_cam) * inv( Rt_lin(:,:,i) ); % sucks

   %  Rt_cam =  cat(1, Rt_glob, [0,0,0,1]) * inv(cat(1, Rt_glob_old,[0,0,0,1]));  

   % ok: never ending story
%    Rt_inv(:,:,i) =  inv(cat(1, Rt_glob_old,[0,0,0,1])) * inv(inv(cat(1, Rt_glob, [0,0,0,1])) * Rt_lin(:,:,i));% logic dictates this:

  end
Rt_cam
%}

% now get sol for the views involved ! to start algorithm with  !!!
  centers2D2 = cat(2, Seg2.Centers, ones (size(Seg2.Centers,1),1))';
  centers2DNew = cam.Kl * centers2D2;
  newids_inOldSegi=[];
  for i=1:size(centers2DNew,2)
    newids_inOldSegi(end+1) =  S_old_img( round( centers2DNew(2,i) ), round( centers2DNew(1,i) ));
  end
  
  sol5 = int32(newids_inOldSegi);
  % here are the ids for the previous frame
%  getKittiErr3dSF ( Seg, ref, cam, N_proj(:,sol5+1), Rt_proj(:,:,sol5+1), 0 ); % is it the prev solution?  
  % project that into right view at t-1, projecting needs centered view!
  sol6 = Project_nonCentered( permute(p2d_, [3,1,2]), cam(1).Kl, cam(1).Tr, cam(1).Rr, 0, 0, S_old.Ids, N_old(1:3,:), Rt_old, centers2D_old, centers2D2, Seg2.Edges, Seg2.Edges, Seg2.Img, Seg2.Img, int32([0:size(N_old,2)]), cam(1).Kr, 0 );

  if plotresults==1
  [N_projFW, Rt_projFW]     = ProjProposals6(Seg2.Edges, cam(1).Kl, N_proj(1:3,:), Rt_proj, centers2D, centers2D2, Seg2.Img );
  for i=1:size(Rt_projFW,3)
%    Rt_projFW(:,:,i) = Rt_cam * Rt_projFW(:,:,i); % hmm: not ok ?    
    Rt_projFW(:,:,i) = inv(Rt_cam) * Rt_projFW(:,:,i); % hmm: ok ?
  end
    plotAnalysis(ref_inv, cam_inv, N_projFW(:,sol5+1), Rt_projFW(:,:,sol5+1), Seg2, u, 20, par, sprintf('AA%03d_sol5Solution_%02d', par.imgNr, par.subImg+1));
  end

% wait - should this not be the solution at t+2? this is for t+1 isn't it?
  sol3 = int32(sol5);
  sol4 = int32(sol6);
  sol5 = Project_nonCentered( permute(p2d_, [3,1,2]), cam(1).Kl, cam(1).Tr, cam(1).Rr, 0, 0, S_old.Ids, N_old(1:3,:), Rt_old, centers2D_old, centers2D2, Seg2.Edges, Seg2.Edges, Seg2.Img, Seg2.Img, int32([0:size(N_old,2)]), cam(1).Kr, 1 );
  sol6 = Project_nonCentered( permute(p2d_, [3,1,2]), cam(1).Kl, cam(1).Tr, cam(1).Rr, 0, 0, S_old.Ids, N_old(1:3,:), Rt_old, centers2D_old, centers2D2, Seg2.Edges, Seg2.Edges, Seg2.Img, Seg2.Img, int32([0:size(N_old,2)]), cam(1).Kr, 2 );

backPack.ids_prevFrame = ids_prevFrame;
backPack.N_proj        = N_proj;
backPack.Rt_proj       = Rt_proj;
backPack.C_proj        = C_proj;
backPack.Rt_cam        = Rt_cam;%
%backPack.Rt_cam        = inv(Rt_cam);% or .. see sol5 above should be once! inverted ONLY not teice check cpp file gHom.double motion
backPack.sol5          = int32(sol5);
backPack.sol6          = int32(sol6);
backPack.sol3          = int32(sol3);
backPack.sol4          = int32(sol4);
backPack.Il            = ImgI{1};
backPack.Ir            = ImgI{2};
