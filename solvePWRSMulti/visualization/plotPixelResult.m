function plotPixelResult( ref, cam, par, dmX, allD, N_lin, Rt_lin, u, ...
                          Seg1, Seg2, Seg3, Seg4, Seg5, Seg6, ref_invI, cam_invI, colors)

if ~exist('colors','var')
  colors = 0;
end                      
%plotSegOverlayConsistent_3F(ref.I(1).I, Seg1, cam(1).I(1).I, Seg2, ref.I(2).I, Seg3, cam(1).I(2).I, Seg4, ref_inv.I(2).I, Seg5, cam_inv.I(2).I, Seg6, 100, par, 0.5, sprintf('3FramesSegi%03d', par.imgNr));                        
                        
if exist('Seg5','var') && exist('ref_inv','var')
  plotSegOverlayConsistent_3F(ref.I(1).I, Seg1, cam(1).I(1).I, Seg2, ref.I(2).I, Seg3, cam(1).I(2).I, Seg4, ref_invI, Seg5, cam_invI, Seg6, 100, par, 0.5, sprintf('Pix_%03d', par.imgNr), colors);
else
  plotSegOverlayConsistent(ref.I(1).I, Seg1, cam(1).I(1).I, Seg2, ref.I(2).I, Seg3, cam(1).I(2).I, Seg4, 100, par, 0.5, sprintf('Pix_%03d', par.imgNr), colors);
end
%plotDepthMapsFromViews (allD(:,:,5), allD(:,:,1), allD(:,:,3), ref_inv.I(2).I, ref.I(1).I, ref.I(2).I, 'PixL_ds03', par);
%plotDepthMapsFromViews (allD(:,:,6), allD(:,:,2), allD(:,:,4), cam_inv.I(2).I, cam(1).I(1).I, cam(1).I(2).I, 'PixR_ds03', par);

strE = 'Pix_Ncut';
plotDataE( ref.I(1).I, cam(1).I(1).I, dmX(:,:,1), dmX(:,:,2), 22, par, sprintf('Label_%s_lr',strE) );
plotDataE( ref.I(1).I, ref.I(2).I, dmX(:,:,3), dmX(:,:,4), 23, par, sprintf('Label_%s_ll',strE) );
plotDataE( cam(1).I(1).I, cam(1).I(2).I, dmX(:,:,5), dmX(:,:,6), 24, par, sprintf('Label_%s_rr',strE) );
plotDataE( ref.I(2).I, cam(1).I(2).I, dmX(:,:,7), dmX(:,:,8), 25, par, sprintf('Label_%s_lrt',strE) );

  if size( dmX, 3) > 9 && exist('ref_inv','var')
    plotDataE( ref.I(1).I, ref_invI, dmX(:,:, 9), dmX(:,:,10), 23, par, sprintf('Label_%s_il',strE) );
    plotDataE( cam.I(1).I, cam_invI, dmX(:,:,11), dmX(:,:,12), 24, par, sprintf('Label_%s_ir',strE) );
    plotDataE( ref_invI,   cam_invI, dmX(:,:,13), dmX(:,:,14), 25, par, sprintf('Label_%s_ilir',strE) );
  end

  plotDepthMapsVC( cat(3, allD(:,:,1), allD(:,:,2) ), cat(3, ref.I(1).I, cam(1).I(1).I), 'PixT0', par, 0, -0.5 );
  plotDepthMapsVC( cat(3, allD(:,:,3), allD(:,:,4) ), cat(3, ref.I(2).I, cam(1).I(2).I), 'PixT1', par, 0, -0.5 );  
if size(allD,3)> 5 && exist('ref_invI','var')
  plotDepthMapsVC( cat(3, allD(:,:,5), allD(:,:,6) ), cat(3, ref_invI, cam_invI), 'PixT2', par, 0, -0.5 );
end
if size(allD,3)> 7
%  plotDepthMapsVC( cat(3, allD(:,:,7), allD(:,:,8) ), cat(3, ref_inv.I(2).I, cam_inv.I(2).I), 'PixLRT3_ds03', par, 0, -0.5 );
end

%plotAnalysis(ref, cam(1), N_lin, Rt_lin, Seg1, u, 20, par, sprintf('%03d_perPixel', par.imgNr));