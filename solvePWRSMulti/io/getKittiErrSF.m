function [occErr, noccErr, epes] = getKittiErrSF ( disp, uDisp, vDisp, p, plotF )

global flow2DGt;
global flow2DGt_noc;

fNr = 23;
plotMyFigs = 0;
if exist ('plotF','var') 
  plotMyFigs = plotF;
end

err2 = disp_error(flow2DGt(:,:,1), -disp, 2);
err3 = disp_error(flow2DGt(:,:,1), -disp, 3);
err4 = disp_error(flow2DGt(:,:,1), -disp, 4);
err5 = disp_error(flow2DGt(:,:,1), -disp, 5);

err2f = flow_error(flow2DGt(:,:,2:end), cat(3, uDisp, vDisp), 2);
err3f = flow_error(flow2DGt(:,:,2:end), cat(3, uDisp, vDisp), 3);
err4f = flow_error(flow2DGt(:,:,2:end), cat(3, uDisp, vDisp), 4);
err5f = flow_error(flow2DGt(:,:,2:end), cat(3, uDisp, vDisp), 5);

err2n = disp_error(flow2DGt_noc(:,:,1), -disp, 2);
err3n = disp_error(flow2DGt_noc(:,:,1), -disp, 3);
err4n = disp_error(flow2DGt_noc(:,:,1), -disp, 4);
err5n = disp_error(flow2DGt_noc(:,:,1), -disp, 5);

err2fn = flow_error(flow2DGt_noc(:,:,2:end), cat(3, uDisp, vDisp), 2);
err3fn = flow_error(flow2DGt_noc(:,:,2:end), cat(3, uDisp, vDisp), 3);
err4fn = flow_error(flow2DGt_noc(:,:,2:end), cat(3, uDisp, vDisp), 4);
err5fn = flow_error(flow2DGt_noc(:,:,2:end), cat(3, uDisp, vDisp), 5);

cprintf('blue', 'DispPix 2/3/4/5 %.3f/%.3f/%.3f/%.3f\n', err2, err3, err4, err5);
cprintf('blue', 'FlowPix 2/3/4/5 %.3f/%.3f/%.3f/%.3f\n', err2f, err3f, err4f, err5f);
cprintf('blue', 'DispPix-noc 2/3/4/5 %.3f/%.3f/%.3f/%.3f\n', err2n, err3n, err4n, err5n);
cprintf('blue', 'FlowPix-noc 2/3/4/5 %.3f/%.3f/%.3f/%.3f\n', err2fn, err3fn, err4fn, err5fn);

epeD     = getEndPointError_Kitti(-disp, flow2DGt(:,:,1));
epe_nocD = getEndPointError_Kitti(-disp, flow2DGt_noc(:,:,1));
epe      = getEndPointError_Kitti(cat(3, uDisp, vDisp), flow2DGt(:,:,2:end));
epe_noc  = getEndPointError_Kitti(cat(3, uDisp, vDisp), flow2DGt_noc(:,:,2:end));
epes.epe = epe;
epes.epe_noc = epe_noc;
epes.epeD = epeD;
epes.epe_nocD = epe_nocD;

cprintf('blue', 'FlowPix-Epe %.3f noc-%.3f\n', epe, epe_noc);
cprintf('blue', 'DispPix-Epe %.3f noc-%.3f\n', epeD, epe_nocD);

occErr.err2 = err2;
occErr.err3 = err3;
occErr.err4 = err4;
occErr.err5 = err5;

occErr.err2f = err2f;
occErr.err3f = err3f;
occErr.err4f = err4f;
occErr.err5f = err5f;

noccErr.err2 = err2n;
noccErr.err3 = err3n;
noccErr.err4 = err4n;
noccErr.err5 = err5n;

noccErr.err2f = err2fn;
noccErr.err3f = err3fn;
noccErr.err4f = err4fn;
noccErr.err5f = err5fn;

%E = (flow2DGt_noc(:,:,1) +disp).* (flow2DGt_noc(:,:,1)>0);
%figure(9); imagesc( E, [-5, 5]), colorbar;

if plotMyFigs
  E = (flow2DGt(:,:,1) +disp).* (flow2DGt(:,:,1)>0);
  f=figure(fNr);set(f, 'visible','off'); imshow( E, [-5, 5]), colorbar, colormap(jet);
  export_fig( sprintf('%s/%03d_disp_Error.png', p.sFolder, p.imgNr), '-m1');

  E = flow_error_img(flow2DGt(:,:,2:end), cat(3, uDisp, vDisp));
  f=figure(fNr+1); set(f, 'visible','off'); imshow(E, [0, 5] ), colorbar, colormap(jet);
  export_fig( sprintf('%s/%03d_flow_Error.png', p.sFolder, p.imgNr), '-m1');
end
