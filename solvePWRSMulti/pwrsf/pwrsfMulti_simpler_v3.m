%%% the main function
function [flow2d, Energy] = pwrsfMulti_simpler_v3 ( ref, cam, par, Seg, N_prop, RT_prop, oracle )

Energy   = 0;
pwrs_ds = 0.03;%0.06 -- smoothness iccv version 

fromProposals    = ~(par.doSeg || par.doJit || par.doEgo); % either directly from proposals or intermediate pwrs step
% local replacement strategy for VC-SF -- implemented suboptimally but at least its there - for sure better to follow the pwrs (iccv13) implementation 
lrp            = par.locRep & fromProposals; % alternative to per-run pwrsf
reduceRNTProps = 1; % default: ON -- reduce similar proposals, combining them into a single one
doAutoIntern   = 1; % default: ON -- estimate above from data term value in optimization
localReduction = 0; % default: OFF -- tries to preselect proposals by solving a simpler problem first
gys=2;gxs=2;        % expansion region for local replacement and reduction -- NOT TOO LARGE !
%
refineLoop   = par.refine;% run refining based on loop default: ON  if 16x16 grid
endlevel     = 8;% refinement in 2^-1 steps, so startlevel= 16, 8 , .., endlevel
%
useFWprojectedProposals = 1; % default:on, 3frame off? 2 frame on 
%
plotXtra     = 0; % debug and eso-pics
plotXtraPrev = 0; % debug and eso-pics
plotresults  = 0; % show results analytic -- segmentation data cost, etc 
plotFinal    = 0; % non analytic -- show results
plotresultsInit= 0;% plot even more ..
plotProposals  = 0; % plot proposal solution
%
sW               = 1.0; % weight for spatial weight of segmentation energy
doVideoProposals = par.usePrevProps; % without use3Frames: add proposals do 2Frame
useThreeFrames   = par.use3Frames;

fprintf('\n\n\n s%.0f, j%.0f, e%.0f  fromProps:%d \n\n\n', par.doSeg, par.doJit, par.doEgo, fromProposals);

colors = 0;       % segments colors constant 
temp = find(isnan(N_prop) | isinf(N_prop) );
while ~isempty(temp)
  N_prop(temp)=N_prop(temp-4);
  temp = find(isnan(N_prop) | isinf(N_prop));
end
temp = find(isnan(RT_prop) | isinf(RT_prop));
while ~isempty(temp)
  RT_prop(temp)=RT_prop(temp-16);
  temp = find(isnan(RT_prop) | isinf(RT_prop));
end

%%%%%%%%%%%% video proposals: load former solution if available %%%%%%%%%%
if  par.subImg>=1 && (doVideoProposals == 1 || useThreeFrames ==1)
  try
    [Rt_cam, N_proj, Rt_proj, C_proj, N_old_lin, Rt_old_lin, S_old_img, centers2D_old, ids_prevFrame, Seg] = ...
      adjustEgomotion (Seg, ref, cam, par, plotXtraPrev, N_prop, RT_prop);
  catch
    cprintf('red','Failed to load videoProposals\n');
    doVideoProposals=0;
    useThreeFrames =0;
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[M,N] = size(Seg.Img);

dt   = par.dt;
ts   = par.ts;
ds   = par.ds;
dj   = par.dj;
tj   = par.tj;
oob  = par.oob;
% not used anymore: can be used if known from the data (DataDefinitionsVC.h)
dD     = ceil(300 * size(ref.I(1).I, 2)/1242); % max disp
maxMot = ceil(450 * size(ref.I(1).I, 2)/1242); % max flow

% seg:  both 'out of order' 
oobC = 0.8; occC = 0.8; %strong interaction with occlusion !

u  = ones(M,N,3);u(:,:,1) = repmat( [1:N],  M, 1 );u(:,:,2) = repmat( [1:M]', 1, N );
p2d_ = permute( reshape ( inv(cam(1).Kl) * reshape(permute (u , [3,1,2]), 3, N*M), [3,M,N]), [2,3,1]);
centers2D = cat(2, Seg.Centers, ones (size(Seg.Centers,1),1))';
nSegs = numel(Seg.Ids);

% segmentation of all other images : the same grid for simplicity
Seg2 = Seg; % well do i not care about different segmentations do i ..
centers2D2 = cat(2, Seg2.Centers, ones (size(Seg2.Centers,1),1))';

nImg={};
nImg{1} = numel(Seg.Ids);
nImg{2} = numel(Seg2.Ids);
nImg{3} = numel(Seg2.Ids);
nImg{4} = numel(Seg2.Ids);
if useThreeFrames || numel(cam) >= 2
  nImg{5} = numel(Seg2.Ids);
  nImg{6} = numel(Seg2.Ids);
end

autos = getNoCPenalties (ref, cam, par, Seg, Seg2, useThreeFrames);
if isfield(cam,'Iold') % exist('ref_inv','var')
%  autos = getAutos(ref, cam, par, Seg, Seg2, useThreeFrames, doAutos, cam.Iold);
  ew    = getEdgeWeights(ref, cam, sW, useThreeFrames, cam.Iold);
else
%  autos = getAutos(ref, cam, par, Seg, Seg2, useThreeFrames, doAutos );
  ew    = getEdgeWeights(ref, cam, sW, useThreeFrames );
end

for kk=0:size(RT_prop,3)/nSegs-1
  fprintf('Proposal %d:\n', kk);
  getKittiErr3dSF ( Seg,ref,cam(1), N_prop  (:, kk*nSegs+1:nSegs*(kk+1)), RT_prop(:,:, kk*nSegs+1:nSegs*(kk+1)) );
end

if plotProposals % pre init print
  plotAnalysis(ref, cam(1), N_prop (:, 1:nSegs), RT_prop(:,:, 1:nSegs), Seg, u, 10, par, 1, sprintf('%03d_props', par.imgNr));
  if size(N_prop,2) >= 2*nSegs
    plotAnalysis(ref, cam(1), N_prop (:, nSegs+1:nSegs*2), RT_prop(:,:, nSegs+1:nSegs*2), Seg, u, 10, par, 1, sprintf('%03d_props', par.imgNr));
  end
  if ~useThreeFrames || isfield(cam,'Iold') % ~exist('ref_inv','var')
    colors = plotSegOverlayConsistent(ref.I(1).I, Seg, cam(1).I(1).I, Seg, ref.I(2).I, Seg, cam(1).I(2).I, Seg, 100, par, 0.5, sprintf('SegProps%03d', par.imgNr), colors);
  else
    colors = plotSegOverlayConsistent_3F(ref.I(1).I, Seg, cam(1).I(1).I, Seg, ref.I(2).I, Seg, cam(1).I(2).I, Seg, cam.Iold{1}, Seg, cam.Iold{2}, Seg, 100, par, 0.5, sprintf('SegProps%03d', par.imgNr), colors);
  end
end
%
% save( sprintf( '%s/perPixNew%03d_%02d.mat', par.sFolder, par.imgNr, par.subImg), 'p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11', 'p12', 'p13','p14','p15','p16','p17');

if fromProposals == true
  [N_lin_prop, Rt_lin_prop, centers2D_prop] = convert_to_nonCentered( Seg, N_prop, RT_prop );
  N_linM = N_lin_prop;Rt_linM = Rt_lin_prop;
  initIds = 0:numel(Seg.Ids)-1;
  oobOracle.oobs        = int32(zeros(size(Seg.Img)));
  oobOracle.oobsFlow    = int32(zeros(size(Seg.Img)));
  oobOracle.oobsFlow2   = int32(zeros(size(Seg.Img)));
  oobOracle.oobsFlowRRt = int32(zeros(size(Seg.Img)));
else % fromProposals == false
  [ oobOracle, initIds, N_lin, Rt_lin ] = ...
    pwrs_init ( ref, cam, Seg, N_prop, RT_prop, ew, dt, pwrs_ds, ts, ...
                dj, tj, dD, oob, maxMot, par, oracle );
  N_linM = N_lin;Rt_linM = Rt_lin;    
  fprintf('Start\n');getKittiErr3dSF ( Seg, ref, cam(1), N_lin, Rt_lin, 0 );
  % not needed for VC
  oobOracle.oobs        = int32(zeros(size(oobOracle.oobs)));
  oobOracle.oobsFlow    = oobOracle.oobs;
  oobOracle.oobsFlow2   = oobOracle.oobs;
  oobOracle.oobsFlowRRt = oobOracle.oobs;
end
%%% assume N_lin, etc contain prop-set, 
sol = int32(initIds);
%%%%%%%%%%

if  plotXtra ==1
  Rt_inv = Rt_lin;
  for i = 1:size(Rt_lin, 3)
    Rt_inv(:,:,i)       =  (Rt_lin(:,:,i)) \ Rt_cam;
  end
  plotAnalysis(cam.Iold{1}, cam.Iold{2}, N_lin(:,sol+1), Rt_inv(:,:,sol+1), Seg, u, 20, par, 0, sprintf('%03d_backWardsFromStart', par.imgNr));%global
end

if doVideoProposals
  % here I pick sol5/sol6 to be good for the previous frame:
  centers2DNew = cam(1).Kl * centers2D2;
  newids_inOldSegi=ones(1,size(centers2DNew,2));
  for i=1:size(centers2DNew,2)
    newids_inOldSegi(i) =  S_old_img( round( centers2DNew(2,i) ), round( centers2DNew(1,i) ));
  end

  if ~reduceRNTProps
    solY = 0:size(N_old_lin,2)-1;
    solX = 0:size(N_proj,2)-1;
  else
    % reduction affects also proposal ids - see below
    [solY, N_old_lin, Rt_old_lin] = reduce_NRT_pix (ref,cam, N_old_lin, Rt_old_lin);
    [solX, N_proj, Rt_proj]       = reduce_NRT_pix (ref,cam, N_proj, Rt_proj);
  end
  
  sol5 = int32( solX(newids_inOldSegi+1) );

  % here are the ids for the previous frame
  fprintf( '_proj(:,sol5+1): so in prev frame - actually says nothing\n');
  getKittiErr3dSF ( Seg2, ref, cam(1), N_proj(:,sol5+1), Rt_proj(:,:,sol5+1), 0 ); % is it the prev solution?
  % project that into right view at t-1, projecting needs centered view!
  sol6 = ProjectLabels( permute(p2d_, [3,1,2]), cam(1).Kl, cam(1).Tr, cam(1).Rr, 0, 0, Seg2.Ids, N_old_lin(1:3,:), Rt_old_lin, centers2D_old, centers2D2, Seg2.Edges, Seg2.Edges, Seg2.Img, Seg2.Img, sol5, cam(1).Kr, 0 );
  % and here the ids for the current frame:
  fprintf( 'X_proj(:,ids_prevFrame):\n');
  getKittiErr3dSF ( Seg, ref, cam(1), N_proj(:,solX(ids_prevFrame)+1), Rt_proj(:,:,solX(ids_prevFrame)+1), 0 ); % works

%%% fourth frame in front: 
%  backPack = projectSolutionBackWards(par, Seg, ref, cam(1), Rt_glob, Seg2);

%  [weight_x7, weight_y7, weight_xy7, weight_ixy7] = getEdgeWeightsForSegmentation (backPack.Il, 35, 0.1);
%  [weight_x8, weight_y8, weight_xy8, weight_ixy8] = getEdgeWeightsForSegmentation (backPack.Ir, 35, 0.1);

  % try to enrich proposal set here
  C_projX = inv(cam(1).Kl) * cat( 1, C_proj, ones(1,size(C_proj,2)) );
  
  % simple or normal as before
  if ~exist('N_lin_prop', 'var')
    maxcol = min( size(N_proj,1), size(N_lin,1) ); % both 3 ? entries
    N_linM  = cat(2, N_lin(1:maxcol,:), N_proj(1:maxcol,:));
    Rt_linM = cat(3, Rt_lin, Rt_proj);
    centers2DM = cat(2, centers2D, C_projX);
    
    % why ids_prevFrame?
    % well ids_prevFrame is projected forward, so the guys in the 
    % current image
    %
    % without the positions are the true forward projected ones ! 
    % also outside the current frame ! -- could also use both
    if useFWprojectedProposals
      propIds = int32(cat ( 1, sol',  solX(ids_prevFrame)+numel(sol) ));
    else
      propIds = int32(cat ( 1, sol',  solX+numel(sol) ));
    end
    sol5 = sol5 + size(N_lin,2);% sol5 & sol6 relative to N_proj, etc.
    sol6 = sol6 + size(N_lin,2);
  else
    N_linM  = cat(2, N_lin_prop,  N_proj);
    Rt_linM = cat(3, Rt_lin_prop, Rt_proj);
    centers2DM = cat(2, centers2D_prop, C_projX);

    % why ids_prevFrame? see above 
    if useFWprojectedProposals
      propIds = int32(cat ( 1, [0:size(N_lin_prop,2)-1]', solX(ids_prevFrame)+size(N_lin_prop,2) ));
    else
      propIds = int32(cat ( 1, [0:size(N_lin_prop,2)-1]', solX+size(N_lin_prop,2) ));
    end

    sol5 = sol5 + size(N_lin_prop,2);% sol5 & sol6 relative to N_proj, etc.
    sol6 = sol6 + size(N_lin_prop,2);
  end

else % sol5 etc. from projection

  % 3 frames but no proposals from past
  if useThreeFrames && ~exist('sol5', 'var')
    Rt_linM_inv = Rt_linM;
    for i = 1:size(Rt_linM, 3)
      Rt_linM_inv(:,:,i) = Rt_linM(:,:,i) \ Rt_cam; % from t to t-1
    end
    sol5 = ProjectLabels( permute(p2d_, [3,1,2]), cam(1).Kl, cam(1).Tr, cam(1).Rr, 0, 0, Seg.Ids, N_linM(1:3,:), Rt_linM_inv, centers2D, centers2D2, Seg.Edges, Seg2.Edges, Seg.Img, Seg2.Img, int32(initIds), cam(1).Kr, 1 );
    sol6 = ProjectLabels( permute(p2d_, [3,1,2]), cam(1).Kl, cam(1).Tr, cam(1).Rr, 0, 0, Seg.Ids, N_linM(1:3,:), Rt_linM_inv, centers2D, centers2D2, Seg.Edges, Seg2.Edges, Seg.Img, Seg2.Img, int32(initIds), cam(1).Kr, 2 );
  end
  
  if ~exist('N_lin_prop', 'var')
    N_linM  = N_lin;
    Rt_linM = Rt_lin;
    centers2DM = centers2D;
    propIds = int32(sol);    
  else
    N_linM  = N_lin_prop ;
    Rt_linM = Rt_lin_prop;
    centers2DM = centers2D_prop;
    propIds = int32(0:size(N_lin_prop,2)-1);
  end
end

% with four frames:
%centers2DM = cat(2, centers2DM , inv(cam.Kl) * cat( 1, backPack.C_proj, ones(1,size(backPack.C_proj,2)) ));
%backPack.sol6 = int32(backPack.sol6 + size(N_linM,2));
%backPack.sol5 = int32(backPack.sol5 + size(N_linM,2));
%sol3          = int32(backPack.sol3 + size(N_linM,2));% might interfere with solitiopn for 10 frame mst be best
%sol4          = int32(backPack.sol4 + size(N_linM,2));
%N_linM  = cat(2, N_linM, backPack.N_proj);
%Rt_linM = cat(3, Rt_linM, backPack.Rt_proj);

% reduction all solutions all normals:
if ~reduceRNTProps
  sol_all = 0:numel(N_linM,2)-1;
else
  [sol_all, N_linM, Rt_linM] = reduce_NRT_pix (ref,cam, N_linM, Rt_linM);
  sol  = int32( sol_all(sol +1) );
  cprintf('Comments', 'Reduced to %05d proposals \n', size(N_linM,2));
  getKittiErr3dSF ( Seg, ref,cam(1), N_linM(:,sol+1), Rt_linM(:,:,sol+1), 0 );
end

if exist('sol5', 'var')
  sol5 = int32( sol_all(sol5+1) );  sol6 = int32( sol_all(sol6+1) );
end

propIds = int32(sol_all(propIds+1));
initIds = int32(sol_all(initIds+1));

sol2 = ProjectLabels( permute(p2d_, [3,1,2]), cam(1).Kl, cam(1).Tr, cam(1).Rr, 0, 0, Seg.Ids, N_linM(1:3,:), Rt_linM, centers2D, centers2D2, Seg.Edges, Seg2.Edges, Seg.Img, Seg2.Img, int32(initIds), cam(1).Kr, 0 );
sol3 = ProjectLabels( permute(p2d_, [3,1,2]), cam(1).Kl, cam(1).Tr, cam(1).Rr, 0, 0, Seg.Ids, N_linM(1:3,:), Rt_linM, centers2D, centers2D2, Seg.Edges, Seg2.Edges, Seg.Img, Seg2.Img, int32(initIds), cam(1).Kr, 1 );
sol4 = ProjectLabels( permute(p2d_, [3,1,2]), cam(1).Kl, cam(1).Tr, cam(1).Rr, 0, 0, Seg.Ids, N_linM(1:3,:), Rt_linM, centers2D, centers2D2, Seg.Edges, Seg2.Edges, Seg.Img, Seg2.Img, int32(initIds), cam(1).Kr, 2, 0, int32(sol2), 1, int32(sol3) );

  p1=cat( 3, ref.I(1).I, cam(1).I(1).I); % input images t0
  p2=cat( 3, ref.I(2).I, cam(1).I(2).I); % input image t1
  p3=cat( 3, Seg.Img, Seg2.Img );        % per pixel segmentation
  p4=cat( 3, Seg.Img, Seg2.Img );        % per pixel segmentation
  % parameters: data, smoothness, smooth motion, disp-max, mot-max, ..
  p5=cat(1, dt, 0.4*ds, ts, dj, tj, dD, maxMot, oobC, occC, Seg.PatchSize);

  if doAutoIntern==1
    p5=cat(1, p5, doAutoIntern, par.vcPottsSeg, par.vcEpsSeg, par.gx, par.gy )
   else
    p5=cat(1, p5, 0, par.vcPottsSeg, par.vcEpsSeg, par.gx, par.gy )
  end

  p6=int32(cat(1, sol(:), sol2(:) )); % init solution at time 0 forall cams
  p7=int32(cat(1, sol3(:), sol4(:))); % init solution at time 1 forall cams
  p8=cat(3, cam(1).Kl, cam(1).Kr ); % internal calibration cameras
  p9=cat(3, cam(1).Rr ); % external calib: rotation cam 1 to N
  p10=cat(3, cam(1).Tr );% external calib: translation cam 1 to N
  p11=N_linM(1:3,:); % proposal normals and 
  p12=Rt_linM;       % proposal rigid motion
  p13=cat(2,centers2D, centers2D); % centers of segments per cam at time 0
  p14=cat(2,centers2D, centers2D); % centers of segments per cam at time 1
  p15=int32(propIds); % which proposal at which center ? prop 77 at center 5: propIds[77]=5
  p16=centers2DM;     % expansion centers - the position of the proposal, contains center 5
  p17={Seg.Edges, Seg.Edges}; % defines segment neighborhoods at t0 all cams
  p18={Seg.Edges, Seg.Edges}; % defines segment neighborhoods at t0 all cams
  p19 = autos.Seg;%autoScores; %all cams first, then t1 all cams , etc
  p20 = ew.x;%T; all cams, times:  edgeweights in x dir
  p21 = ew.y;%T;
  p22 = ew.xy;%T;
  p23 = ew.yx;%T;
  p24=int32(oobOracle.oobs);
  p25=int32(oobOracle.oobsFlow);
  p26=int32(oobOracle.oobsFlowRRt);
  %     save( sprintf( '%s/perSegNew%03d_%02d.mat', par.sFolder, par.imgNr, par.subImg), 'p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11', 'p12', 'p13','p14','p15','p16','p17','p18','p19','p20','p21','p22','p23','p24','p25', 'p26');

  if numel(cam) >=2
    p1=cat( 3, p1, cam(2).I(1).I);
    p2=cat( 3, p2, cam(2).I(2).I);

    solx2 = ProjectLabels( permute(p2d_, [3,1,2]), cam(2).Kl, cam(2).Tr, cam(2).Rr, 0, 0, Seg.Ids, N_linM(1:3,:), Rt_linM, centers2D, centers2D2, Seg.Edges, Seg2.Edges, Seg.Img, Seg2.Img, int32(initIds), cam(2).Kr, 0 );
    solx4 = ProjectLabels( permute(p2d_, [3,1,2]), cam(2).Kl, cam(2).Tr, cam(2).Rr, 0, 0, Seg.Ids, N_linM(1:3,:), Rt_linM, centers2D, centers2D2, Seg.Edges, Seg2.Edges, Seg.Img, Seg2.Img, int32(initIds), cam(2).Kr, 2, 0, int32(solx2), 1, int32(sol3) );
    p6=int32(cat(1, p6, solx2(:) ));
    p7=int32(cat(1, p7, solx4(:) ));
    
    p8 =cat(3, cam(1).Kl, cam(1).Kr, cam(2).Kr );
    p9 =cat(3, cam(1).Rr, cam(2).Rr );
    p10=cat(3, cam(1).Tr, cam(2).Tr );

    % same segmentation
    p3=cat( 3, p3, Seg2.Img );
    p4=cat( 3, p4, Seg2.Img );
    p13=cat(2,centers2D, centers2D, centers2D);
    p14=cat(2,centers2D, centers2D, centers2D);
    p17={Seg.Edges, Seg.Edges, Seg.Edges};
    p18={Seg.Edges, Seg.Edges, Seg.Edges};
  end
  
  if plotresultsInit==1 % debugging
    p5s = p5;p5s(end-1) = 0;p5s(end) = 0;
    [newIds, allD, dmX, dataS, energyInfo] = Seg_3SM( p1,p2,p3,p4,p5s,p6,p7,p8,p9, p10,p11,p12,p13,p14,p15,p16,p17,p18,p19, p20,p21,p22,p23,p24,p25,p26);
    plotAnalysis(ref, cam(1), N_linM(:,newIds(1 : nImg{1} )+1), Rt_linM(:,:,newIds(1 : nImg{1} )+1), Seg, u, 10, par, 0, sprintf('SegInit%d_level_%.1f', par.imgNr, -1));
    if isfield(cam,'Iold') % exist('cam_inv','var') && exist('ref_inv','var')
      plotSegResult (cam, ref, Seg, par, dmX, allD, dataS, 'init', cam.Iold{1}, cam.Iold{2});%, dataA
    else
      if exist('dataA','var')
        plotSegResult (cam, ref, Seg, par, dmX, allD, dataS, 'init', dataA );
      else
        plotSegResult (cam, ref, Seg, par, dmX, allD, dataS, 'init' );
      end
    end
  end
  
  if useThreeFrames && isfield(cam,'Iold') % exist ('ref_inv', 'var')
    p27 = cat( 3, cam.Iold{1}, cam.Iold{2} );
    p28 = cat( 3, Seg.Img, Seg2.Img );
    p29 = int32(cat(1, sol5(:), sol6(:)));
    p30 = cat(2,centers2D, centers2D);
    p31 = {Seg.Edges, Seg.Edges};
    p32 = Rt_cam;
  
    [newIds, allD, dmX, dataS, energyInfo] = ..., dataA ] = ...
      Seg_3SM( p1,p2,p3,p4,p5,p6,p7,p8,p9, p10,p11,p12,p13,p14,p15,p16,p17,p18,p19, p20,p21,p22,p23,p24,p25,p26,p27,p28,p29, p30,p31,p32);
  else % idea is to just kepp local props instead of all of them here
    if localReduction
      % reduction step, keep local props only - discard rest
      p5s = p5;p5s(end-1) = gxs;p5s(end) = gys;
      [newIds] = Seg_3SM( p1,p2,p3,p4,p5s,p6,p7,p8,p9, p10,p11,p12,p13,p14,p15,p16,p17,p18,p19, p20,p21,p22,p23,p24,p25,p26);
      % now based on this set new
      fprintf('reduced proposal ids to %04d proposals\n', numel(unique(newIds(1:nImg{1}))) );
      getKittiErr3dSF ( Seg, ref,cam(1), N_linM(:,newIds(1:nImg{1})+1), Rt_linM(:,:,newIds(1:nImg{1})+1), 0 );
      
      p15 = int32(newIds(1:nImg{1})); % could also keep all ? 
      p16 = p16(:, 1:nImg{1});
      %    save( sprintf( '%s/simple2F_%03d_%02d.mat', par.sFolder, par.imgNr, par.subImg), 'p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11', 'p12', 'p13','p14','p15','p16','p17','p18','p19','p20','p21','p22','p23','p24','p25', 'p26');
    end
    [newIds, allD, dmX, dataS, energyInfo] = ..., dataA ] = ...
      Seg_3SM( p1,p2,p3,p4,p5,p6,p7,p8,p9, p10,p11,p12,p13,p14,p15,p16,p17,p18,p19, p20,p21,p22,p23,p24,p25,p26);
    % is this energy lower equal or higher if all proposals are evaluated?
  end

if plotresults==1
  if isfield(cam,'Iold') %exist('cam_inv','var') && exist('ref_inv','var')
    plotSegResult (cam, ref, Seg, par, dmX, allD, dataS, 'SegInit', cam.Iold{1}, cam.Iold{2});%dataA,
  else
    if exist('dataA','var')
      plotSegResult (cam, ref, Seg, par, dmX, allD, dataS, 'SegInit', dataA );
    else
      plotSegResult (cam, ref, Seg, par, dmX, allD, dataS, 'SegInit' );
    end
  end
end

% save proposal set from thie video frame:
if par.saveProposals
  if ~exist(par.tempFolder,'dir')
    mkdir(sprintf('%s/', par.tempFolder));
  end
   N_old = N_linM; Rt_old = Rt_linM;S_old = Seg;oldIds = newIds(1:nImg{1});
   save( sprintf( '%s/RoughSolution%03d_%02d.mat', par.tempFolder, par.imgNr, par.subImg), 'S_old', 'N_old', 'Rt_old', 'oldIds');
end

while lrp >0
  % how to call ? : new input, thus 
  vcIds={};tmp=0;
  for i=1:numel(nImg)
    vcIds{i} = newIds(1+tmp : nImg{i}+tmp );
    tmp = tmp + nImg{i};
  end
  
  p5s = p5;p5s(end-1) = gxs;p5s(end) = gys;% not needed over too large areas
  p6=int32(cat(1, vcIds{1}(:), vcIds{2}(:) ));
  p7=int32(cat(1, vcIds{3}(:), vcIds{4}(:)));
  p15 = int32(vcIds{1});

  %  save( sprintf( '%s/replace_%03d_%02d.mat', par.sFolder, par.imgNr, par.subImg), 'p1','p2','p3','p4','p5s','p6','p7','p8','p9','p10','p11', 'p12', 'p13','p14','p15','p16','p17','p18','p19','p20','p21','p22','p23','p24','p25', 'p26');
  getKittiErr3dSF ( Seg, ref, cam, N_linM(:,1+newIds(1:numel(Seg.Ids))), Rt_linM(:,:,1+newIds(1:numel(Seg.Ids))), 0 );

  % Seg_3SM_mcOff_replace
  if useThreeFrames
    p29 = int32(cat(1, vcIds{5}(:), vcIds{6}(:)) );
    [newIds,~,~,~,~,N_linM, Rt_linM, cadd] = Seg_3SM_locreplace( p1,p2,p3,p4,p5s,p6,p7,p8,p9, p10,p11,p12,p13,p14,p15,p16,p17,p18,p19, p20,p21,p22,p23,p24,p25,p26,p27,p28,p29, p30,p31,p32);
  else
    %%% needs to return new N, RT !!!!!!!
    [newIds,~,~,~,~,N_linM, Rt_linM, cadd] = Seg_3SM_locreplace( p1,p2,p3,p4,p5s,p6,p7,p8,p9, p10,p11,p12,p13,p14,p15,p16,p17,p18,p19, p20,p21,p22,p23,p24,p25,p26);    
  end
  centers2DM = cat( 2, centers2DM, cadd);
  if reduceRNTProps
    [sol_all, N_linM, Rt_linM, ia] = reduce_NRT_pix (ref,cam, N_linM, Rt_linM);
    newIds  = int32( sol_all(newIds+1) );
    centers2DM = centers2DM (:,ia);
    cprintf('Comments', 'Reduced to %05d proposals \n', size(N_linM,2));
  end
  
  p11=N_linM(1:3,:); % proposal normals and
  p12=Rt_linM;       % proposal rigid motion
  p16=centers2DM;
  fprintf('LRP it: %d\n', lrp);
  lrp = lrp -1;
end

  %%%%%%%%%% REFINE %%%%%%%
  level = double(int32(Seg.PatchSize/2));
  while refineLoop==1 && level >= endlevel
    vcIds={};tmp=0;
    for i=1:numel(nImg)
      vcIds{i} = newIds(1+tmp : nImg{i}+tmp );
      tmp = tmp + nImg{i};
    end
    if numel(cam)>1
      xtraCam={};
      xtraCam{1} = vcIds{3};
      xtraCam{2} = vcIds{6};
      
      vcIds{3} = vcIds{4};
      vcIds{4} = vcIds{5};
      vcIds(5:6) = [];
    end

    % plot segmentation
    if plotresults==1
      
      projImg1 = Seg.Img;
      for i=1:max(Seg.Img(:))+1 projImg1(Seg.Img == i-1) = vcIds{1}(i);end;
      projImg2 = Seg.Img;
      for i=1:max(Seg.Img(:))+1 projImg2(Seg.Img == i-1) = vcIds{2}(i);end;
      projImg3 = Seg.Img;
      for i=1:max(Seg.Img(:))+1 projImg3(Seg.Img == i-1) = vcIds{3}(i);end;
      projImg4 = Seg.Img;
      for i=1:max(Seg.Img(:))+1 projImg4(Seg.Img == i-1) = vcIds{4}(i);end;
      Seg1=Seg;Seg1.Img = int32(projImg1);Seg2=Seg;Seg2.Img = int32(projImg2);
      Seg3=Seg;Seg3.Img = int32(projImg3);Seg4=Seg;Seg4.Img = int32(projImg4);
      
      if numel(vcIds)==6 %exist('newIds5', 'var')
        projImg5 = Seg.Img;
        for i=1:max(Seg.Img(:))+1 projImg5(Seg.Img == i-1) = vcIds{5}(i);end;
        projImg6 = Seg.Img;
        for i=1:max(Seg.Img(:))+1 projImg6(Seg.Img == i-1) = vcIds{6}(i);end;
        Seg5=Seg;Seg5.Img = int32(projImg5);Seg6=Seg;Seg6.Img = int32(projImg6);
      end

      if exist('xtraCam', 'var')
        projImg2x = Seg.Img;
        for i=1:max(Seg.Img(:))+1 projImg2x(Seg.Img == i-1) = xtraCam{1}(i);end;
        projImg4x = Seg.Img;
        for i=1:max(Seg.Img(:))+1 projImg4x(Seg.Img == i-1) = xtraCam{2}(i);end;
        Seg2x=Seg;Seg2x.Img = int32(projImg2x);Seg4x=Seg;Seg4x.Img = int32(projImg4x);
      end

      if ~useThreeFrames || ~isfield(cam,'Iold') % ~exist('ref_inv','var')
        colors = plotSegOverlayConsistent(ref.I(1).I, Seg1, cam(1).I(1).I, Seg2, ref.I(2).I, Seg3, cam(1).I(2).I, Seg4, 100, par, 0.5, sprintf('SegLev%02d_Cubes%03d', Seg.PatchSize, par.imgNr), colors);
      else
        colors = plotSegOverlayConsistent_3F(ref.I(1).I, Seg1, cam(1).I(1).I, Seg2, ref.I(2).I, Seg3, cam(1).I(2).I, Seg4, cam.Iold{1}, Seg5, cam.Iold{2}, Seg6, 100, par, 0.5, sprintf('SegLev%02d_Cubes%03d', Seg.PatchSize, par.imgNr), colors);
      end

      plotAnalysis(ref, cam(1), N_linM(:,vcIds{1}+1), Rt_linM(:,:,vcIds{1}+1), Seg, u, 10, par, 0, sprintf('SegLev_smo3%d_level_%.1f', par.imgNr, level));
    end

    % scale up by a factor of 2
    SegOld = Seg;%oldIds = int32(vcIds{i});

    Seg = SegmentImageCube( ref.I(1).I, level, par );
    Seg = correctSegCenters( Seg );
    Seg = setWeights_patchSmooth( Seg, cam(1).Kr );
    centers2D    = cat(2, Seg.Centers, ones (size(Seg.Centers,1),1))';
    Seg2 = Seg;
    % adjust indices of seg2prop, nothing else:
%     if isfield(cam,'Iold') % exist('ref_inv','var')
%       autos = getAutos(ref, cam, par, Seg, Seg, useThreeFrames, doAutos, cam.Iold);
%     else
%       autos = getAutos(ref, cam, par, Seg, Seg, useThreeFrames, doAutos );
%     end
    autos = getNoCPenalties (ref, cam, par, Seg, Seg2, useThreeFrames);
    p19 = autos.Seg;%autoScores;

    % map old solution to new/finer segmentation
    sol = cell (numel(vcIds),1);
    for k=1:numel(vcIds)
      for i = 1:size(Seg.pixCenters, 1)
        pc  = round( Seg.pixCenters (i,:) );
        sol{k}(end+1)  = vcIds{k} ( SegOld.Img( pc(1),pc(2)) +1);
      end
    end
    getKittiErr3dSF ( Seg, ref,cam, N_linM(:,sol{1}+1), Rt_linM(:,:,sol{1}+1), 0 ); % old

    centers2D2 = centers2D;    
    nImg={};
    nImg{1} = numel(Seg.Ids);
    nImg{2} = numel(Seg2.Ids);
    nImg{3} = numel(Seg2.Ids);
    nImg{4} = numel(Seg2.Ids);
    if useThreeFrames || numel(cam) >= 2
      nImg{5} = numel(Seg2.Ids);
      nImg{6} = numel(Seg2.Ids);
    end
    
    p3=cat( 3, Seg.Img, Seg.Img );
    p4=cat( 3, Seg.Img, Seg.Img );
    p5=cat(1,dt, 0.4*ds, ts, dj, tj, dD, maxMot, oobC, occC, Seg.PatchSize);
    if doAutoIntern==1
      p5=cat(1, p5, doAutoIntern, par.vcPottsSeg, par.vcEpsSeg, par.gx, par.gy )
     else
      p5=cat(1, p5, 0, par.vcPottsSeg, par.vcEpsSeg, par.gx, par.gy )
    end

    p6=int32(cat(1, sol{1}(:), sol{2}(:) ));
    p7=int32(cat(1, sol{3}(:), sol{4}(:)));

    p13=cat(2,centers2D, centers2D);
    p14=cat(2,centers2D, centers2D);

    % test: no gain
testNew =0;
if testNew
    projImg1 = Seg.Img;
    for i=1:max(Seg.Img(:))+1 projImg1(Seg.Img == i-1) = sol{1}(i);end;
    [C_projX, newIdsT1 ] = expansionProposals (ref, cam(1), N_linM, Rt_linM, centers2DM, 1, Seg, projImg1, centers2D, centers2D, sol{2}, 1 );
    centerIds      = int32(sol{1});
    centers2DM_pic = cat(2, centers2D, C_projX);
    centerIds      = int32(cat(2, centerIds, newIdsT1));
    % t+1 backwards:
    [C_projX, newIdsT1 ] = expansionProposals (ref, cam(1), N_linM, Rt_linM, centers2DM, 1, Seg, projImg1, centers2D, centers2D2, sol{3}, 3);
    centers2DM_pic = cat(2, centers2DM_pic, C_projX);
    centerIds      = int32(cat(2, centerIds, newIdsT1));
    %
    p15 = int32(centerIds);
    p16 = centers2DM_pic;
else
    p15 = int32(sol{1});
    p16 = centers2D;
end
    p17={Seg.Edges, Seg.Edges};
    p18={Seg.Edges, Seg.Edges};
    
    if useThreeFrames
      p27 = cat( 3, cam.Iold{1}, cam.Iold{2} );
      p28 = cat( 3, Seg.Img, Seg.Img );
      p29 = int32(cat(1, sol{5}(:), sol{6}(:)) );
      p30 = cat(2,centers2D, centers2D );
      p31 = {Seg.Edges, Seg.Edges};
      p32 = Rt_cam;
      
%    save( sprintf( '%s/perSegNew3f_%03d_%02d.mat', par.sFolder, par.imgNr, par.subImg), 'p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11', 'p12', 'p13','p14','p15','p16','p17','p18','p19','p20','p21','p22','p23','p24','p25', 'p26', 'p27','p28','p29','p30','p31','p32');
      
      [newIds, allD, dmX, dataS, energyInfo] = ..., dataA ] = ...
        Seg_3SM( p1,p2,p3,p4,p5,p6,p7,p8,p9, p10,p11,p12,p13,p14,p15,p16,p17,p18,p19, p20,p21,p22,p23,p24,p25,p26,p27,p28,p29, p30,p31,p32);
    else
      [newIds, allD, dmX, dataS, energyInfo] = ..., dataA ] = ...
        Seg_3SM( p1,p2,p3,p4,p5,p6,p7,p8,p9, p10,p11,p12,p13,p14,p15,p16,p17,p18,p19, p20,p21,p22,p23,p24,p25,p26);
    end
    level = double(int32(level/2));
  end
  %%%%%%%%%% END REFINE %%%%%%%  
  
  %%% plot and convert data for per pixel step:
  vcIds={};tmp=0;
  for i=1:numel(nImg)
    vcIds{i} = newIds(1+tmp : nImg{i}+tmp );
    tmp = tmp + nImg{i};
  end
  if numel(cam)>1
    xtraCam={};
    xtraCam{1} = vcIds{3};
    xtraCam{2} = vcIds{6};
    
    vcIds{3} = vcIds{4};
    vcIds{4} = vcIds{5};
    vcIds(5:6) = [];
  end

  getKittiErr3dSF ( Seg, ref,cam(1), N_linM(:,vcIds{1}+1), Rt_linM(:,:,vcIds{1}+1), 0 ); % old

projImg1 = Seg.Img;projImg2 = Seg.Img;projImg3 = Seg.Img;projImg4 = Seg.Img;
if numel(vcIds)==6 %exist('newIds5', 'var')
  projImg5 = Seg.Img;  projImg6 = Seg.Img;
end
if exist('xtraCam', 'var')
  projImg2x = Seg.Img;  projImg4x = Seg.Img;
end

for i=1:max(Seg.Img(:))+1 
  copySegids = Seg.Ids{i}+1;%Seg.Img == i-1;
  projImg1(copySegids) = vcIds{1}(i);
  projImg2(copySegids) = vcIds{2}(i);
  projImg3(copySegids) = vcIds{3}(i);
  projImg4(copySegids) = vcIds{4}(i);
  if numel(vcIds)==6
    projImg5(copySegids) = vcIds{5}(i);
    projImg6(copySegids) = vcIds{6}(i);
  end
  if exist('xtraCam', 'var')
    projImg2x(copySegids) = xtraCam{1}(i);%Seg.Img == i-1;
    projImg4x(copySegids) = xtraCam{2}(i);
  end
end

Seg1=Seg;Seg1.Img = int32(projImg1);Seg2=Seg;Seg2.Img = int32(projImg2);
Seg3=Seg;Seg3.Img = int32(projImg3);Seg4=Seg;Seg4.Img = int32(projImg4);
if numel(vcIds)==6  
  Seg5=Seg;Seg5.Img = int32(projImg5);Seg6=Seg;Seg6.Img = int32(projImg6);
end
if exist('xtraCam', 'var')
  Seg2x=Seg;Seg2x.Img = int32(projImg2x);Seg4x=Seg;Seg4x.Img = int32(projImg4x);
end

% 4 frames: as above
% projImg7 = Seg.Img;
% for i=1:max(Seg.Img(:))+1 projImg7(Seg.Ids{i}+1) = newIds7(i);end;
% projImg8 = Seg.Img;
% for i=1:max(Seg.Img(:))+1 projImg8(Seg.Ids{i}+1) = newIds8(i);end;

if plotresults==1
  selectedBW = zeros(size(Seg.Img));
  for i=1:max(Seg.Img(:))+1 selectedBW(Seg.Img == i-1) = vcIds{1}(i);end;
  f= figure(21);set(f, 'visible','off');imshow(selectedBW,[]), colormap(jet), colorbar;
  export_fig( sprintf('%s/pickedFromLastFrame%03d.png', par.sFolder, par.imgNr), '-m1');
  close(f);
end
N_lin = N_linM;Rt_lin = Rt_linM; % for consistent continuation

if  plotXtra ==1
  plotAnalysis(ref, cam(1), N_lin, Rt_lin, Seg, u, 20, par, 0, sprintf('%03d_Start', par.imgNr));
end

if plotresults==1 && refineLoop==1
    plotAnalysis(ref, cam(1), N_linM(:,vcIds{1}+1), Rt_linM(:,:,vcIds{1}+1), Seg, u, 10, par, 0, sprintf('SegEnd_%03d', par.imgNr));
if isfield(cam,'Iold') % exist('cam_inv','var') && exist('ref_inv','var')
  plotSegResult (cam, ref, Seg, par, dmX, allD, dataS, 'SegEnd', cam.Iold{1}, cam.Iold{2});%dataA, 
else
  if exist('dataA','var')
    plotSegResult (cam, ref, Seg, par, dmX, allD, dataS, 'SegEnd', dataA );
  else
    plotSegResult (cam, ref, Seg, par, dmX, allD, dataS, 'SegEnd' );
  end
end
end

if plotresults==1
  if ~useThreeFrames || ~isfield(cam,'Iold') % exist('ref_inv','var')
    colors = plotSegOverlayConsistent(ref.I(1).I, Seg1, cam(1).I(1).I, Seg2, ref.I(2).I, Seg3, cam(1).I(2).I, Seg4, 100, par, 0.5, sprintf('SegCubes%03d', par.imgNr), colors);
  else
    colors = plotSegOverlayConsistent_3F(ref.I(1).I, Seg1, cam(1).I(1).I, Seg2, ref.I(2).I, Seg3, cam(1).I(2).I, Seg4, cam.Iold{1}, Seg5, cam.Iold{2}, Seg6, 100, par, 0.5, sprintf('SegCubes%03d', par.imgNr), colors);
  end
  %plotSegOverlayConsistent_4F(ref.I(1).I, Seg1, cam.I(1).I, Seg2, ref.I(2).I, Seg3, cam.I(2).I, Seg4, ref_inv.I(2).I, Seg5, cam_inv.I(2).I, Seg6, backPack.Il, Seg7, backPack.Ir, Seg8, 100, par, 0.5, sprintf('Seg250Cubes%03d', par.imgNr));
end

%%% expand regions from other views as well, if distinctive use as proposal
  increaseSet = true;  
  centers2DM_pic = centers2D;
  centerIds      = int32( vcIds{1} );
  if increaseSet
    % from right to left
    [C_projX, newIdsT1 ] = expansionProposals (ref, cam(1), N_linM, Rt_linM, centers2DM, 1, Seg, projImg1, centers2D, centers2D2, vcIds{2}, 1 );
    centers2DM_pic = cat(2, centers2DM_pic, C_projX);
    centerIds      = int32(cat(1, centerIds, newIdsT1));

    % from past forwards
    if numel(vcIds) >5
      [C_projX, newIdsT1 ] = expansionProposals (ref, cam(1), N_linM, Rt_linM, centers2DM, 1, Seg, projImg1, centers2D, centers2D2, vcIds{5}, 5, Rt_cam );
      centers2DM_pic = cat(2, centers2DM_pic, C_projX);
      centerIds      = int32(cat(1, centerIds, newIdsT1));
    end
    
    % t+1 backwards:
    [C_projX, newIdsT1 ] = expansionProposals (ref, cam(1), N_linM, Rt_linM, centers2DM, 1, Seg, projImg1, centers2D, centers2D2, vcIds{3}, 3);
    centers2DM_pic = cat(2, centers2DM_pic, C_projX);
    centerIds      = int32(cat(1, centerIds, newIdsT1));
    
    % t+2 backwards:   % check pwrsfMulti_simpler_v2 for testing & debug t+2 backwards
    if numel(vcIds) >7 && exist('backPack','var')
      [C_projX, newIdsT1 ] = expansionProposals (ref, cam(1), N_linM, Rt_linM, centers2DM, 1, Seg, projImg1, centers2D, centers2D2, vcIds{7}, 7, backPack.Rt_cam );
      centers2DM_pic = cat(2, centers2DM_pic, C_projX);
      centerIds      = int32(cat(1, centerIds, newIdsT1));
    end
  end
%%%

%     p1=cat( 3, ref.I(1).I, cam(1).I(1).I); 
%     p2=cat( 3, ref.I(2).I, cam(1).I(2).I); 
     p3  = cat( 3, projImg1, projImg2 );
     p4  = cat( 3, projImg3, projImg4 );
     p5  = cat(1,dt, 0.4*ds, ts, dj, tj, dD, maxMot, oobC, occC, Seg.PatchSize, Seg.PatchSize+2, par.segWeight, par.vcPottsPix, par.vcEpsPix, doAutoIntern);
     p6  = cat(3, cam(1).Kl, cam(1).Kr );
     p7  = cat(3, cam(1).Rr );    
     p8  = cat(3, cam(1).Tr );
     p9  = N_linM(1:3,:);
     p10 = Rt_linM;
     p11 = centers2DM_pic;
     p12 = int32(centerIds);
     p13 = autos.Pix;
     
     p14 = ew.x;
     p15 = ew.y;
     p16 = ew.xy;
     p17 = ew.yx;
  if useThreeFrames
    p18 = cat( 3, cam.Iold{1}, cam.Iold{2} );
    p19  = cat( 3, projImg5, projImg6 );
    p20 = Rt_cam;
    [allSolutions, allD, dmX, energyInfoPix ] = ... 
      Pix_3SM(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20);
  else
    
    if numel(cam)>=2
      p3  = cat( 3, projImg1, projImg2, projImg2x );
      p4  = cat( 3, projImg3, projImg4, projImg4x );
      p6  = cat(3, cam(1).Kl, cam(1).Kr, cam(2).Kr );
      p7  = cat(3, cam(1).Rr, cam(2).Rr );
      p8  = cat(3, cam(1).Tr, cam(2).Tr );
    end
    [allSolutions, allD, dmX, energyInfoPix ] = ...
      Pix_3SM(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17);
  end
  
  SegNew=Seg;SegNew.Img = int32(allSolutions(:,:,1));SegNew.Ids={};%!! Ids must be empty
  for i=1:max(SegNew.Img(:))+1 SegNew.Ids{i} = find(SegNew.Img == i-1)-1;end;
  flow2d =  reconstruc2dFlowHom( ref, cam(1), N_lin, Rt_lin, SegNew, 0 );
  
  %%% rest: plotting error evaluations, energies, etc.
  getKittiErr3dSF ( SegNew, ref, cam(1), N_lin, Rt_lin, 0 );

  if  plotresults == 1
    Seg1=Seg;Seg1.Img = int32(allSolutions(:,:,1));Seg2=Seg;Seg2.Img = int32(allSolutions(:,:,2));Seg3=Seg;Seg3.Img = int32(allSolutions(:,:,3));Seg4=Seg;Seg4.Img = int32(allSolutions(:,:,4));
    % ?? PLOTTING:
    for i=1:max(Seg1.Img(:))+1
      Seg1.Ids{i} = find(Seg1.Img == i-1)-1;
      Seg2.Ids{i} = find(Seg2.Img == i-1)-1;
      Seg3.Ids{i} = find(Seg3.Img == i-1)-1;
      Seg4.Ids{i} = find(Seg4.Img == i-1)-1;
    end;
    
    if size( allSolutions, 3) > 4
      Seg5=Seg;Seg5.Img = int32(allSolutions(:,:,5));Seg6=Seg;Seg6.Img = int32(allSolutions(:,:,6));
    end
    if size( allSolutions, 3) > 6
      Seg7=Seg;Seg7.Img = int32(allSolutions(:,:,7));Seg8=Seg;Seg8.Img = int32(allSolutions(:,:,8));
    end
    
    if exist('Seg5','var') && isfield(cam,'Iold') % exist('ref_inv','var')
      plotPixelResult( ref, cam, par, dmX, allD, N_lin, Rt_lin, u, ...
        Seg1, Seg2, Seg3, Seg4, Seg5, Seg6, cam.Iold{1}, cam.Iold{2}, colors);
    else
      plotPixelResult( ref, cam, par, dmX, allD, N_lin, Rt_lin, u, Seg1, Seg2, Seg3, Seg4, colors );
    end
  end
  if plotFinal==1
    plotAnalysis(ref, cam(1), N_linM, Rt_linM, SegNew, u, 10, par, 0, sprintf('%03d_perPixel', par.imgNr));
  end
  
  fprintf('Seg\n Energy: %.2f; \n non-sub: %.2f, \n non-sol: %.2f, \n allNodes: %.2f, \n all_Edges: %.2f, \n nRuns: %.2f\n, non-sub pct: %.2f\n, non-sol pct: %.2f\n', ...
    energyInfo(1), energyInfo(2), energyInfo(3), energyInfo(4), energyInfo(5), energyInfo(6), energyInfo(2)/energyInfo(5)*100, energyInfo(3)/energyInfo(5)*100 );
  
  fprintf('Pix\n Energy: %.2f; \n non-sub: %.2f, \n non-sol: %.2f, \n allNodes: %.2f, \n all_Edges: %.2f, \n nRuns: %.2f\n, non-sub pct: %.2f\n, non-sol pct: %.2f\n', ...
    energyInfoPix(1), energyInfoPix(2), energyInfoPix(3), energyInfoPix(4), energyInfoPix(5), energyInfoPix(6), energyInfoPix(2)/energyInfoPix(5)*100, energyInfoPix(3)/energyInfoPix(5)*100 );
  
  energystr = sprintf('Energy: %.2f\nnon-sub: %.2f\nnon-sol: %.2f\nallNodes: %.2f\nall_Edges: %.2f\nnRuns: %.2f\nnon-sub-pct: %.2f\nnon-sol-pct: %.2f\n', ...
    energyInfo(1), energyInfo(2), energyInfo(3), energyInfo(4), energyInfo(5), energyInfo(6), energyInfo(2)/energyInfo(5)*100, energyInfo(3)/energyInfo(5)*100 );
  
  energystr = sprintf('%sEnergy: %.2f\nnon-sub: %.2f\nnon-sol: %.2f\nallNodes: %.2f\nall_Edges: %.2f\nnRuns: %.2f\nnon-sub-pct: %.2f\nnon-sol-pct: %.2f\n', ...
    energystr, energyInfoPix(1), energyInfoPix(2), energyInfoPix(3), energyInfoPix(4), energyInfoPix(5), energyInfoPix(6), energyInfoPix(2)/energyInfoPix(5)*100, energyInfoPix(3)/energyInfoPix(5)*100 );

energystr = sprintf('%sUniqueMVP: %.2f\n', energystr, numel(unique(SegNew.Img)) );

fid = fopen(sprintf('%s/ENERGY_K%03d_%02d_%s.txt', par.sFolder, par.imgNr, par.subImg, date), 'w', 'n');
if fid~=-1
  fwrite(fid, energystr, 'char');
  fclose(fid);
end

return;
