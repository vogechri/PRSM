%%% generate the proposal set by fitting flow and motion parts into the
%%% segments. 
function [ N_prop, RT_prop, oracle ] = generateProposals(par, cam, ref, Seg )

flowstereo2d_SGM = 1;
flowstereo2d     = 0;
N_prop = [];
RT_prop = [];
% lacks parallelism -- which is very possible 2 flows, init_Seg, etc. 
innerIts = 10;warps = 4;pyrscale = 0.9;
dataterm = 2;% 1:CENSUS  2:CSAD 0: SAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flowstereo2d
  % should be parallel:
  stereoT_2d = TGV_flow( 1.25/255, 9.3333, warps, pyrscale, ref.I(1).I, cam.I(1).I, innerIts, 2, 1, 1);stereoT_2d = stereoT_2d(:,:,1);
  if par.computeRflow
    [flowL_2d, flowR_2d] = TGV_flowDouble( 1.25/255, 12.333, warps, pyrscale, ref.I(1).I, ref.I(2).I, cam.I(1).I, cam.I(2).I, innerIts, 2, dataterm, 0);
    % alternatively:
%    flowR_2d   = TGV_flow( 1.25/255, 12.333, warps, pyrscale, cam.I(1).I, cam.I(2).I, innerIts, 2, dataterm, 0);
%    flowL_2d   = TGV_flow( 1.25/255, 12.333, warps, pyrscale, ref.I(1).I, ref.I(2).I, innerIts, 2, dataterm, 0);
  else
    flowL_2d   = TGV_flow( 1.25/255, 12.333, warps, pyrscale, ref.I(1).I, ref.I(2).I, innerIts, 2, dataterm, 0);
  end

  if par.computeRflow
    [N_lin, Rt_lin] = initSeg_2dFlowTest(ref, cam, Seg, stereoT_2d, stereoT_2d, flowL_2d, flowR_2d, 1, par.fitSegs);
  else
    [N_lin, Rt_lin] = initSeg_2dFlowTest(ref, cam, Seg, stereoT_2d, stereoT_2d, flowL_2d, flowL_2d, 0, par.fitSegs);
  end
  
  fprintf('full2d - lpFit\n');
  [occErr, noccErr, epes] = getKittiErr3dSF ( Seg, ref, cam, N_lin, Rt_lin );
  fprintf('full2d\n');
  [occErr, noccErr, epes] = getKittiErrSF ( stereoT_2d(:,:,1), flowL_2d(:,:,1), flowL_2d(:,:,2) );
  
  kittiStr = sprintf('DispPix-occ 2/3/4/5 %.3f & %.3f & %.3f & %.3f\nFlowPix-occ 2/3/4/5 %.3f & %.3f & %.3f & %.3f\n', occErr.err2, occErr.err3, occErr.err4, occErr.err5, occErr.err2f, occErr.err3f, occErr.err4f, occErr.err5f);
  kittiStr = sprintf('%s\nDispPix-noc 2/3/4/5 %.3f & %.3f & %.3f & %.3f\nFlowPix-noc 2/3/4/5 %.3f & %.3f & %.3f & %.3f\n', kittiStr, noccErr.err2, noccErr.err3, noccErr.err4, noccErr.err5, noccErr.err2f, noccErr.err3f, noccErr.err4f, noccErr.err5f);
  kittiStr = sprintf('%s\nDispEPE %.3f & %.3f\nFlowEPE %.3f & %.3f\n', kittiStr, epes.epe_nocD, epes.epeD, epes.epe_noc, epes.epe);
  
  fid = fopen(sprintf('%s/RESULTS_F%03d_%02d_%s.txt', par.sFolder, par.imgNr, par.subImg, date), 'w', 'n');
  if fid~=-1
    fwrite(fid, kittiStr, 'char');
    fclose(fid);
  end
  
  N_prop  = cat( 2, N_prop,  N_lin);
  RT_prop = cat( 3, RT_prop, Rt_lin);
   
  if ~flowstereo2d_SGM && par.generateMoreProposals% 2nd p set with less different mvps
    if par.computeRflow
      [N_lin, Rt_lin] = initSeg_2dFlowTest(ref, cam, Seg, stereoT_2d, stereoT_2d, flowL_2d, flowR_2d, 1, par.fitSegs/2 );
    else
      [N_lin, Rt_lin] = initSeg_2dFlowTest(ref, cam, Seg, stereoT_2d, stereoT_2d, flowL_2d, flowL_2d, 0, par.fitSegs/2 );      
    end
    N_prop  = cat( 2, N_prop,  N_lin);
    RT_prop = cat( 3, RT_prop, Rt_lin);
  end
end %   if 2dflowstereo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flowstereo2d_SGM
  
  if exist('stereoT_2d','var')
    %this or
    %    stereoT_2d_new = stereoT_2d(:,:,1);
    %    stereoT_2d = -getDisparitySGM_proposal(cam, ref, -stereoT_2d_new, 0, 160, 3 );
    % just sgm:
    stereoT_2d = -getDisparitySGM_proposal(cam, ref, 0*ones(size(ref.I(1).I)), 0, 160, 3, 100 );
  else
    stereoT_2d = -getDisparitySGM_proposal(cam, ref, 0*ones(size(ref.I(1).I)), 0, 160, 3, 100 );

    if par.computeRflow
      [flowL_2d, flowR_2d] = TGV_flowDouble( 1.25/255, 12.333, warps, pyrscale, ref.I(1).I, ref.I(2).I, cam.I(1).I, cam.I(2).I, innerIts, 2, dataterm, 0);
    % alternatively: 
%      flowR_2d   = TGV_flow( 1.25/255, 12.333, warps, pyrscale, cam.I(1).I, cam.I(2).I, innerIts, 2, dataterm, 0);
%      flowL_2d   = TGV_flow( 1.25/255, 12.333, warps, pyrscale, ref.I(1).I, ref.I(2).I, innerIts, 2, dataterm, 0);      
    else
      flowL_2d   = TGV_flow( 1.25/255, 12.333, warps, pyrscale, ref.I(1).I, ref.I(2).I, innerIts, 2, dataterm, 0);
    end
  end
  
  % all i need in the end: - and cam.K and .. loading
  % last 1 : use both views to fit
  
  if par.computeRflow
    % delivers more proposals but 'more of the same' proposals, which are
    % reduced later -- in general it is smarter to assign a center and
    % size? of a proposal and grow segments 'near' pixel level but 'throw'
    % stupid stuff
    % alternative multiple sparse inits plus lk-style homography fit
    % e.g. two best fits?
    [N_lin, Rt_lin] = initSeg_2dFlowTest(ref, cam, Seg, stereoT_2d, stereoT_2d, flowL_2d, flowR_2d, 1, par.fitSegs);
  else
    % delivers more proposals but 'more of the same' proposals
    [N_lin, Rt_lin] = initSeg_2dFlowTest(ref, cam, Seg, stereoT_2d, stereoT_2d, flowL_2d, flowL_2d, 0, par.fitSegs);
  end
  
  global doKittiErrors;
  if doKittiErrors
    fprintf('full2d-SGM\n');
    [oErr, noErr] = getKittiErrSF ( stereoT_2d(:,:,1), flowL_2d(:,:,1), flowL_2d(:,:,2) );
    fprintf('full2d-SGM - lpFit\n');
    [oErr, noErr] = getKittiErr3dSF ( Seg, ref, cam, N_lin, Rt_lin );
  end 
  
  N_prop  = cat( 2, N_prop,  N_lin);
  RT_prop = cat( 3, RT_prop, Rt_lin);
  
  if ~flowstereo2d && par.generateMoreProposals % 2nd p set with less different mvps - works also with less
    if par.computeRflow
      [N_lin, Rt_lin] = initSeg_2dFlowTest(ref, cam, Seg, stereoT_2d, stereoT_2d, flowL_2d, flowR_2d, 1, par.fitSegs/2 );
    else
      [N_lin, Rt_lin] = initSeg_2dFlowTest(ref, cam, Seg, stereoT_2d, stereoT_2d, flowL_2d, flowL_2d, 0, par.fitSegs/2 );
    end
    N_prop  = cat( 2, N_prop,  N_lin);
    RT_prop = cat( 3, RT_prop, Rt_lin);
    % next level so nice is 1000, 500 and 200 (with 100 even a bit better but ..)
    % note that duplicates get removed so actually proposal set is
    % 1000+500+200 = 1700 proposals
    if par.computeRflow
      [N_lin, Rt_lin] = initSeg_2dFlowTest(ref, cam, Seg, stereoT_2d, stereoT_2d, flowL_2d, flowR_2d, 1, par.fitSegs/5 );
    else
      [N_lin, Rt_lin] = initSeg_2dFlowTest(ref, cam, Seg, stereoT_2d, stereoT_2d, flowL_2d, flowL_2d, 0, par.fitSegs/5 );
    end
    N_prop  = cat( 2, N_prop,  N_lin);
    RT_prop = cat( 3, RT_prop, Rt_lin);
  end
end %   if 2dflowstereo

%remove nan's!
bad_id = ceil(find(isnan(N_prop))/4);
bad_id = bad_id(1:3:end);
for i=1:numel(bad_id)
  N_prop( 1:3, bad_id ) = N_prop( 1:3, bad_id-1 );
end
bad_id = ceil(find(isnan(RT_prop))/16);
bad_id = bad_id(1:3:end);
for i=1:numel(bad_id)
  RT_prop( 1:3, 4, bad_id ) =  RT_prop( 1:3, 4, bad_id-1 );
end

oracle.stereo = stereoT_2d;
oracle.flowL = flowL_2d;
if exist('flowR_2d','var')
  oracle.flowR = flowR_2d;
else
  oracle.flowR = flowl_2d;
end
end
