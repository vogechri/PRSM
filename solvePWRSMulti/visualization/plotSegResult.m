function plotSegResult (cam, ref, Seg, par, dmX, allD, dataS, ending, ref_inv, cam_inv, dataA )

  %  plotAnalysis(ref, cam, N_lin(:,newIds+1), Rt_lin(:,:,newIds+1), Seg, u, 20, par, sprintf('%03d_segRough', par.imgNr));

  if ~exist('ending','var')
    ending = '';
  end
  
  if numel(cam) >=2
    plotDepthMapsVC( cat(3, allD(:,:,1), allD(:,:,2) ), cat(3, ref.I(1).I, cam(1).I(1).I), 'SegLR0_ds03', par, 0, -0.5 );
    plotDepthMapsVC( cat(3, allD(:,:,4), allD(:,:,5) ), cat(3, ref.I(2).I, cam(1).I(2).I), 'SegLR1_ds03', par, 0, -0.5 );
    if size(allD,3) > 5 && numel(cam) >=2
      plotDepthMapsVC( cat(3, allD(:,:, 3), allD(:,:,6) ), cat(3, cam(2).I(1).I, cam(2).I(2).I), 'SegLRCam2_ds03', par, 0, -0.5 );
    end
	strE = sprintf( '%s', ending);
    plotDataE( ref.I(1).I, cam(1).I(1).I, dmX(:,:,1), dmX(:,:,2), 22, par, sprintf('Label_%s_l0r0',strE) );
    plotDataE( ref.I(1).I, cam(2).I(1).I, dmX(:,:,3), dmX(:,:,4), 22, par, sprintf('Label_%s_l0r0C2',strE) );    

    plotDataE( ref.I(1).I, ref.I(2).I, dmX(:,:,5), dmX(:,:,6), 23, par, sprintf('Label_%s_l0l1',strE) );
    plotDataE( cam(1).I(1).I, cam(1).I(2).I, dmX(:,:,7), dmX(:,:,8), 24, par, sprintf('Label_%s_r0r1',strE) );
    plotDataE( cam(2).I(1).I, cam(2).I(2).I, dmX(:,:,9), dmX(:,:,10), 24, par, sprintf('Label_%s_r0r1C2',strE) );

    plotDataE( ref.I(2).I, cam(1).I(2).I, dmX(:,:,11), dmX(:,:,12), 25, par, sprintf('Label_%s_l1r1',strE) );
    plotDataE( ref.I(2).I, cam(2).I(2).I, dmX(:,:,13), dmX(:,:,14), 25, par, sprintf('Label_%s_l1r1C2',strE) );

    strE = 'VCOnlyData';

    plotDataS( Seg, ref.I(1).I, cam(1).I(1).I, dataS(:,:,1), dataS(:,:,2), 22, par, sprintf('Data_%s_l0r0',strE) );
    plotDataS( Seg, ref.I(1).I, cam(2).I(1).I, dataS(:,:,3), dataS(:,:,4), 22, par, sprintf('Data_%s_l0r0C2',strE) );    

    plotDataS( Seg, ref.I(1).I,    ref.I(2).I,    dataS(:,:,5), dataS(:,:, 6), 23, par, sprintf('Data_%s_l0l1',strE) );
    plotDataS( Seg, cam(1).I(1).I, cam(1).I(2).I, dataS(:,:,7), dataS(:,:, 8), 24, par, sprintf('Data_%s_r0r1',strE) );
    plotDataS( Seg, cam(2).I(1).I, cam(2).I(2).I, dataS(:,:,9), dataS(:,:,10), 24, par, sprintf('Data_%s_r0r1C2',strE) );

    plotDataS( Seg, ref.I(2).I, cam(1).I(2).I, dataS(:,:,11), dataS(:,:,12), 25, par, sprintf('Data_%s_l1r1',strE) );
    plotDataS( Seg, ref.I(2).I, cam(2).I(2).I, dataS(:,:,13), dataS(:,:,14), 25, par, sprintf('Data_%s_l1r1C2',strE) );
    
  else
  strE = sprintf( '%s', ending);
  plotDataE( ref.I(1).I, cam(1).I(1).I, dmX(:,:,1), dmX(:,:,2), 22, par, sprintf('Label_%s_l0r0',strE) );
  plotDataE( ref.I(1).I, ref.I(2).I, dmX(:,:,3), dmX(:,:,4), 23, par, sprintf('Label_%s_l0l1',strE) );
  plotDataE( cam(1).I(1).I, cam(1).I(2).I, dmX(:,:,5), dmX(:,:,6), 24, par, sprintf('Label_%s_r0r1',strE) );
  plotDataE( ref.I(2).I, cam(1).I(2).I, dmX(:,:,7), dmX(:,:,8), 25, par, sprintf('Label_%s_l1r1',strE) );
  
  strE = sprintf( 'VCOnlyData%s', ending);
  plotDataS( Seg, ref.I(1).I, cam(1).I(1).I, dataS(:,:,1), dataS(:,:,2), 22, par, sprintf('Data_%s_l0r0',strE) );
  plotDataS( Seg, ref.I(1).I, ref.I(2).I, dataS(:,:,3), dataS(:,:,4), 23, par, sprintf('Data_%s_l0l1',strE) );
  plotDataS( Seg, cam(1).I(1).I, cam(1).I(2).I, dataS(:,:,5), dataS(:,:,6), 24, par, sprintf('Data_%s_r0r1',strE) );
  plotDataS( Seg, ref.I(2).I, cam(1).I(2).I, dataS(:,:,7), dataS(:,:,8), 25, par, sprintf('Data_%s_l1r1',strE) );  
  
  if exist('dataA','var')
    strE = sprintf( 'VCOnlyAuto%s', ending);
    plotDataS( Seg, ref.I(1).I, cam(1).I(1).I, dataA(:,:,1), dataA(:,:,2), 22, par, sprintf('Auto_%s_l0r0',strE) );
    plotDataS( Seg, ref.I(1).I, ref.I(2).I, dataA(:,:,3), dataA(:,:,4), 23, par, sprintf('Auto_%s_l0l1',strE) );
    plotDataS( Seg, cam(1).I(1).I, cam(1).I(2).I, dataA(:,:,5), dataA(:,:,6), 24, par, sprintf('Auto_%s_r0r1',strE) );
    plotDataS( Seg, ref.I(2).I, cam(1).I(2).I, dataA(:,:,7), dataA(:,:,8), 25, par, sprintf('Auto_%s_l1r1',strE) );
  end
  if size(dmX,3)>10 && exist('ref_inv','var')
    plotDataE( ref_inv, ref_inv, dmX(:,:,9), dmX(:,:,10), 26, par, sprintf('Label_%s_l0lli',strE) );
    plotDataE( cam_inv, cam_inv, dmX(:,:,11), dmX(:,:,12), 27, par, sprintf('Label_%s_r0ri',strE) );
    plotDataE( ref_inv, cam_inv, dmX(:,:,13), dmX(:,:,14), 28, par, sprintf('Label_%s_liri',strE) );
  end
  if size(dmX,3)>14 && exist('backPack','var')
    plotDataE( ref.I(2).I, backPack.Il,  dmX(:,:,15), dmX(:,:,16), 26, par, sprintf('Label_%s_l1l2',strE) );
    plotDataE( cam(1).I(2).I, backPack.Ir,  dmX(:,:,17), dmX(:,:,18), 27, par, sprintf('Label_%s_r1r2',strE) );
    plotDataE( backPack.Il, backPack.Ir, dmX(:,:,19), dmX(:,:,20), 28, par, sprintf('Label_%s_l2r2',strE) );
  end

  plotDepthMapsVC( cat(3, allD(:,:,1), allD(:,:,2) ), cat(3, ref.I(1).I, cam(1).I(1).I), sprintf( 'SegT%d%s', 0, ending), par, 0, -0.5 );
  plotDepthMapsVC( cat(3, allD(:,:,3), allD(:,:,4) ), cat(3, ref.I(2).I, cam(1).I(2).I), sprintf( 'SegT%d%s', 1, ending), par, 0, -0.5 );
%  plotDepthMapsVC( cat(3, allD(:,:,5), allD(:,:,6) ), cat(3, ref.I(2).I, cam.I(2).I), 'SegR0_ds03', par, 0, -0.5 ); 
%  plotDepthMapsVC( cat(3, allD(:,:,7), allD(:,:,8) ), cat(3, cam.I(1).I, cam.I(2).I), 'SegR1_ds03', par, 0, -0.5 );  
if size(allD,3) > 5 && exist('ref_inv','var')
  plotDepthMapsVC( cat(3, allD(:,:, 5), allD(:,:,6) ), cat(3, ref_inv, cam_inv), sprintf( 'SegT%d%s', 2, ending), par, 0, -0.5 );
%  plotDepthMapsVC( cat(3, allD(:,:,11), allD(:,:,12) ), cat(3, cam_inv.I(1).I, cam_inv.I(2).I), 'SegRR2_ds03', par, 0, -0.5 );
%  plotDepthMapsVC( cat(3, allD(:,:,13), allD(:,:,14) ), cat(3, cam.I(1).I, cam.I(2).I), 'SegRR_ds03', par, 0, -0.5 );
end
  end  
% if allD > 4  
%   plotDepthMapsFromViews (allD(:,:,5), allD(:,:,1), allD(:,:,3), ref_inv.I(2).I, ref.I(1).I, ref.I(2).I, 'SegL_ds03', par);
%   plotDepthMapsFromViews (allD(:,:,6), allD(:,:,2), allD(:,:,4), cam_inv.I(2).I, cam.I(1).I, cam.I(2).I, 'SegR_ds03', par);
% else
%   plotDepthMapsFromViews (zeros(size(allD(:,:,1))), allD(:,:,1), allD(:,:,3), ref.I(1).I, ref.I(1).I, ref.I(2).I, 'VCOnlyPropSegL_ds03', par, 0, -0.5);
%   plotDepthMapsFromViews (zeros(size(allD(:,:,1))), allD(:,:,2), allD(:,:,4), cam.I(1).I, cam.I(1).I, cam.I(2).I, 'VCOnlyPropSegR_ds03', par, 0, -0.5);  
% end

  %flow2d =  reconstruc2dFlowHom( ref, cam, N_linM(:,newIds+1), Rt_linM(:,:,newIds+1), Seg, 0 );
  %return
  
  %plotDepthMapsFromViews (dL1i, dR0, dL1, ref_inv.I(2).I, ref.I(1).I, ref.I(2).I, 'SegL_ds03', par);
  %plotDepthMapsFromViews (dR1i, dR0, dR1, cam_inv.I(2).I, cam.I(1).I, cam.I(2).I, 'SegR_ds03', par);
  
  %plotDepthMapsFromViews (dR0, dL1, dR1, cam.I(1).I, ref.I(2).I, cam.I(2).I, 'SegiRough_ds04', par);
  %plotDepthMapsFromViews (dR0, dL1i, dR1i, cam.I(1).I, ref_inv.I(2).I, cam_inv.I(2).I, 'SegiRoughBack_ds04', par, 0, -0.4);
  
  %plotDepthMapsFromViews (dR0, dL2, dR2, cam.I(1).I, backPack.Il, backPack.Ir, 'SegiRoughNext_ds04', par, 0, -0.4);

%   for i = 1:size(N_lin, 2)
%     Rt_lin_inv(:,:,i) =  inv(Rt_lin(:,:,i)) * Rt_cam;
%   end
%   plotAnalysis(ref_inv, cam_inv, N_lin(:,newIds+1), Rt_lin_inv(:,:,newIds+1), Seg, u, 20, par, sprintf('%03d_segRoughBack', par.imgNr));
end
