function [ N_prop, RT_prop, oracle ] = generateProposals_load(par, cam, ref, Seg, dtaThresh, dispMax)

flowstereo2d     = 1;
flowstereo2d_cen = 1;
N_prop  = [];
RT_prop = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dtaFolder = '../../../data/';% linux
dtaFolder = '../../../../data/variational';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  kittiStr = '';  
  if flowstereo2d_cen

    % last pure 2d try
%  if par.subImg==10
%      stereoT_2d = disp_read_Kitti( sprintf('../../../Geiger/stereoL/Par_pf0.90_w5_i20_l9.33_stt0.00_ce0.15_r2_Dta1KITTI_stereo_%d.png', par.imgNr ) );
%   else
%     stereoT_2d = disp_read_Kitti( sprintf('../../../Geiger/stereoL/Par_pf0.90_w5_i20_l9.33_stt0.00_ce0.15_r2_Dta1KITTI_stereo_%d_%02d.png', par.imgNr, par.subImg ) );
%  end

%  flowL_2d   = flow_read_Kitti( sprintf('../../../Geiger/flowL/Par_pf0.90_w5_i20_l12.33_stt0.00_ce1.25_r2_Dta1KITTI_%03d_%02d.png', par.imgNr, par.subImg));
%  flowR_2d   = flow_read_Kitti( sprintf('../../../Geiger/flowR/Par_pf0.90_w5_i20_l12.33_stt0.00_ce1.25_r2_Dta1KITTI_%03d_%02d.png', par.imgNr, par.subImg));

  if par.subImg==10
    stereoT_2d = disp_read_Kitti( sprintf('%s/disp/%06d_%02d.png',  dtaFolder, par.imgNr, par.subImg));
    flowL_2d   = flow_read_Kitti( sprintf('%s/flowL/%06d_%02d.png', dtaFolder, par.imgNr, par.subImg));
    flowR_2d   = flow_read_Kitti( sprintf('%s/flowR/%06d_%02d.png', dtaFolder, par.imgNr, par.subImg));
  else
    stereoT_2d = disp_read_Kitti( sprintf('../../../Geiger/stereoL/Par_pf0.90_w5_i20_l9.33_stt0.00_ce0.15_r2_Dta1KITTI_stereo_%d_%02d.png', par.imgNr, par.subImg ) );
    flowL_2d   = flow_read_Kitti( sprintf('../../../Geiger/flowL/Par_pf0.90_w5_i20_l12.33_stt0.00_ce1.25_r2_Dta1KITTI_%03d_%02d.png', par.imgNr, par.subImg));
    flowR_2d   = flow_read_Kitti( sprintf('../../../Geiger/flowR/Par_pf0.90_w5_i20_l12.33_stt0.00_ce1.25_r2_Dta1KITTI_%03d_%02d.png', par.imgNr, par.subImg));
  end

  if 0
    % need to do that first to have a nice oob thingy as well:
    DispImg__x1 = getDisparitySGM_proposal(cam, ref, u,stereoT_2d_new, 0, 160, 3 );
    disp_write(DispImg__x1, sprintf('I:/cvogel/dispImagesTraining_videoNew/disparityEstimate%03d_%02d.png', par.imgNr, par.subImg ));
    stereoT_2d = -DispImg__x1;
  end

    % all i need in the end: - and cam.K and .. loading
    % last 1 : use both views to fit
%    [N_lin, Rt_lin] = initSeg_2dFlow(ref, cam, Seg, stereoT_2d, stereoT_2d, flowL_2d, flowR_2d, par, 1);
  [N_lin, Rt_lin] = initSeg_2dFlowTest(ref, cam, Seg, stereoT_2d, stereoT_2d, flowL_2d, flowR_2d, 1, par.fitSegs);

%    [N_lin, Rt_lin] = initSeg_2dFlowTukey(ref, cam, Seg, stereoT_2d, stereoT_2d, flowL_2d, flowR_2d, 0);	
	
%     M = size(cam.I(1).I,1);N = size(cam.I(1).I,2);
%     u  = ones(M,N,3);u(:,:,1) = repmat( [1:N],  M, 1 );u(:,:,2) = repmat( [1:M]', 1, N );
%     plotAnalysis(ref, cam(1), N_lin, Rt_lin, Seg, u, 10, par, 1, sprintf('%03d_genPropD', par.imgNr));
%{   
    fprintf('full2d\n');
    [oErr, noErr] = getKittiErrSF ( stereoT_2d(:,:,1), flowL_2d(:,:,1), flowL_2d(:,:,2) );
    kittiStr = sprintf('%s\nDispPix-occ 2/3/4/5 %.3f & %.3f & %.3f & %.3f\nFlowPix-occ 2/3/4/5 %.3f & %.3f & %.3f & %.3f\n', kittiStr, oErr.err2, oErr.err3, oErr.err4, oErr.err5, oErr.err2f, oErr.err3f, oErr.err4f, oErr.err5f);
    kittiStr = sprintf('%s\nDispPix-noc 2/3/4/5 %.3f & %.3f & %.3f & %.3f\nFlowPix-noc 2/3/4/5 %.3f & %.3f & %.3f & %.3f\n', kittiStr, noErr.err2, noErr.err3, noErr.err4, noErr.err5, noErr.err2f, noErr.err3f, noErr.err4f, noErr.err5f);

    fprintf('full2d - lpFit\n');
    [oErr, noErr] = getKittiErr3dSF ( Seg, ref, cam, N_lin, Rt_lin );
    kittiStr = sprintf('%s\nDispPix-occ 2/3/4/5 %.3f & %.3f & %.3f & %.3f\nFlowPix-occ 2/3/4/5 %.3f & %.3f & %.3f & %.3f\n', kittiStr, oErr.err2, oErr.err3, oErr.err4, oErr.err5, oErr.err2f, oErr.err3f, oErr.err4f, oErr.err5f);
    kittiStr = sprintf('%s\nDispPix-noc 2/3/4/5 %.3f & %.3f & %.3f & %.3f\nFlowPix-noc 2/3/4/5 %.3f & %.3f & %.3f & %.3f\n', kittiStr, noErr.err2, noErr.err3, noErr.err4, noErr.err5, noErr.err2f, noErr.err3f, noErr.err4f, noErr.err5f);
  %}  
    N_prop  = cat( 2, N_prop,  N_lin);
    RT_prop = cat( 3, RT_prop, Rt_lin);
    
   if ~flowstereo2d % 2nd p set with less different mvps - works also with less
     [N_lin, Rt_lin] = initSeg_2dFlowTest(ref, cam, Seg, stereoT_2d, stereoT_2d, flowL_2d, flowR_2d, 1, par.fitSegs/2 );
     N_prop  = cat( 2, N_prop,  N_lin);
     RT_prop = cat( 3, RT_prop, Rt_lin);
   end    
    
  end %   if 2dflowstereo
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if flowstereo2d

    % last pure 2d try
%    stereoT_2d = disp_read_Kitti( sprintf('../../../Geiger/stereo/disparityEstimate%03d_%02d.png', par.imgNr, par.subImg ) );
%    flowL_2d   = flow_read_Kitti( sprintf('../../../Geiger/flowL/Par_pf0.90_w5_i20_l12.33_stt0.00_ce1.25_r2_Dta1KITTI_%03d_10.png', par.imgNr));
%    flowR_2d   = flow_read_Kitti( sprintf('../../../Geiger/flowR/Par_pf0.90_w5_i20_l12.33_stt0.00_ce1.25_r2_Dta1KITTI_%03d_10.png', par.imgNr));

  if exist('stereoT_2d','var')
    stereoT_2d_new = stereoT_2d(:,:,1);
    stereoT_2d = -getDisparitySGM_proposal(cam, ref, -stereoT_2d_new, 0, 160, 3 );
  else  
%  stereoT_2d = disp_read_Kitti( sprintf('%s/stereo/disparityEstimate%03d_%02d.png', dtaFolder, par.imgNr, par.subImg ) );
    stereoT_2d = disp_read_Kitti( sprintf('%s/dispSGM/disparityEstimate%03d_%02d.png', dtaFolder, par.imgNr, par.subImg ) );
    stereoT_2d = -getDisparitySGM_proposal(cam, ref, zeros(size(Seg.Img)), 0, 160, 3 );    
  end

    if 0
      % need to do that first to have a nice oob thingy as well:
      DispImg__x1 = getDisparitySGM_proposal(cam, ref, u,stereoT_2d_new, 0, 160, 3 );
      disp_write(DispImg__x1, sprintf('I:/cvogel/dispImagesTraining_videoNew/disparityEstimate%03d_%02d.png', par.imgNr, par.subImg ));
      stereoT_2d = -DispImg__x1;
    end

    % all i need in the end: - and cam.K and .. loading
    % last 1 : use both views to fit
%    [N_lin, Rt_lin] = initSeg_2dFlow(ref, cam, Seg, stereoT_2d, stereoT_2d, flowL_2d, flowR_2d, par, 1);
     [N_lin, Rt_lin] = initSeg_2dFlowTest(ref, cam, Seg, stereoT_2d, stereoT_2d, flowL_2d, flowR_2d, 1, par.fitSegs);	
%    [N_lin, Rt_lin] = initSeg_2dFlowTukey(ref, cam, Seg, stereoT_2d, stereoT_2d, flowL_2d, flowR_2d, 1);  

%     M = size(cam.I(1).I,1);N = size(cam.I(1).I,2);
%     u  = ones(M,N,3);u(:,:,1) = repmat( [1:N],  M, 1 );u(:,:,2) = repmat( [1:M]', 1, N );
%     plotAnalysis(ref, cam(1), N_lin, Rt_lin, Seg, u, 10, par, 1, sprintf('%03d_genPropSGM', par.imgNr));
%{
    fprintf('full2d\n');
    [oErr, noErr] = getKittiErrSF ( stereoT_2d(:,:,1), flowL_2d(:,:,1), flowL_2d(:,:,2) );
    kittiStr = sprintf('%s\nDispPix-occ 2/3/4/5 %.3f & %.3f & %.3f & %.3f\nFlowPix-occ 2/3/4/5 %.3f & %.3f & %.3f & %.3f\n', kittiStr, oErr.err2, oErr.err3, oErr.err4, oErr.err5, oErr.err2f, oErr.err3f, oErr.err4f, oErr.err5f);
    kittiStr = sprintf('%s\nDispPix-noc 2/3/4/5 %.3f & %.3f & %.3f & %.3f\nFlowPix-noc 2/3/4/5 %.3f & %.3f & %.3f & %.3f\n', kittiStr, noErr.err2, noErr.err3, noErr.err4, noErr.err5, noErr.err2f, noErr.err3f, noErr.err4f, noErr.err5f);

    fprintf('full2d - lpFit\n');
    [oErr, noErr] = getKittiErr3dSF ( Seg, ref, cam, N_lin, Rt_lin );
    kittiStr = sprintf('%s\nDispPix-occ 2/3/4/5 %.3f & %.3f & %.3f & %.3f\nFlowPix-occ 2/3/4/5 %.3f & %.3f & %.3f & %.3f\n', kittiStr, oErr.err2, oErr.err3, oErr.err4, oErr.err5, oErr.err2f, oErr.err3f, oErr.err4f, oErr.err5f);
    kittiStr = sprintf('%s\nDispPix-noc 2/3/4/5 %.3f & %.3f & %.3f & %.3f\nFlowPix-noc 2/3/4/5 %.3f & %.3f & %.3f & %.3f\n', kittiStr, noErr.err2, noErr.err3, noErr.err4, noErr.err5, noErr.err2f, noErr.err3f, noErr.err4f, noErr.err5f);
  %}  
    N_prop  = cat( 2, N_prop,  N_lin);
    RT_prop = cat( 3, RT_prop, Rt_lin);
    
   if ~flowstereo2d_cen % 2nd p set with less different mvps - works also with less
     [N_lin, Rt_lin] = initSeg_2dFlowTest(ref, cam, Seg, stereoT_2d, stereoT_2d, flowL_2d, flowR_2d, 1, par.fitSegs/2 );
     N_prop  = cat( 2, N_prop,  N_lin);
     RT_prop = cat( 3, RT_prop, Rt_lin);
   end    
    
  end %   if 2dflowstereo

%  fid = fopen(sprintf('%s/ProposalError%03d_%s.txt', par.sFolder, par.imgNr,  date), 'w', 'n');
%  if fid~=-1
%    fwrite(fid, kittiStr, 'char');
%    fclose(fid);
%  end

  RT_prop(1:3,4,find(abs( RT_prop(3,4,:) )>100)) = 0; % there is no hyper speed 
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
  oracle.flowR = flowR_2d;
  
end
