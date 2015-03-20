function calltestTGV( inr, cEps, lambda, warps, pyramid_factor, innerIts, ring, dataTerm, doStereo )

% flow csad, 12.333, a0 1.5 ?!
% stereo cen, 9.333 or better 10, a0 5 g on ?!
  addpath('./io/');
  addpath('../KittiIO/');
  path(path,'../devkit/matlab/');
  doTesting = 0;
  p.imageName = sprintf('KITTI_%03d', inr );

  global Linux_;Linux_=0;
  
  if Linux_
    if ~doStereo
      [I1, I2, flowGT, flowGT_noc, p.imageName] = ...
        loadKITTIImage(inr, '../data/data_stereo_flow/data_stereo_flow/training/',  1);
    else
      [I1, I2, flowGT, flowGT_noc, p.imageName] = ...
        loadKITTIImage_Stereo(inr, '../data/data_stereo_flow/data_stereo_flow/training/',  1);
    end
    if doTesting
      [I1, I2, p.imageName] = ...
        loadKITTIImageTesting(inr, '../data/data_stereo_flow/data_stereo_flow/data_stereo_flow/testing/',  1 , doStereo);
    end
  else
    if ~doStereo
      [I1, I2, flowGT, flowGT_noc, p.imageName] = ...
        loadKITTIImage(inr, '../data/data_stereo_flow/training/',  10, 1);
    else
      [I1, I2, flowGT, flowGT_noc, p.imageName] = ...
        loadKITTIImage_Stereo(inr, '../data/data_stereo_flow/training/', 10, 1);
    end
    if doTesting
      [I1, I2, p.imageName] = ...
        loadKITTIImageTesting(inr, '../../../data/data_stereo_flow/testing/',  doStereo);
    end
  end

  flow = TGV_flow(cEps, lambda, warps, pyramid_factor, I1, I2, innerIts, ring, dataTerm, doStereo);
  
  if ~doStereo
%     if doTesting
%       flow_write(cat(3, flow, zeros(size(flow,1), size(flow,2))), sprintf('%s%s.png', p.sFolder, p.imageName) );
%     else
%       flow_write(cat(3, flow, flowGT(:,:,3)), sprintf('%s%s.png', p.sFolder, p.imageName) );
%     end
%    flow_write(cat(3, flow, flowGT(:,:,3)), sprintf('%s%s.png', p.sFolder, p.imageName) );
    err2f = flow_error(flowGT(:,:,1:end), flow, 2);
    err3f = flow_error(flowGT(:,:,1:end), flow, 3);
    err4f = flow_error(flowGT(:,:,1:end), flow, 4);
    err5f = flow_error(flowGT(:,:,1:end), flow, 5);
    
    err2fn = flow_error(flowGT_noc(:,:,1:end), flow, 2);
    err3fn = flow_error(flowGT_noc(:,:,1:end), flow, 3);
    err4fn = flow_error(flowGT_noc(:,:,1:end), flow, 4);
    err5fn = flow_error(flowGT_noc(:,:,1:end), flow, 5);

    epeErr  = getEndPointError(cat(3, flow, ones(size(I1))), flowGT);
    epeErrN = getEndPointError(cat(3, flow, ones(size(I1))), flowGT_noc);    
  else
%    disp_write(-flow(:,:,1), sprintf('%s%s.png', p.sFolder, p.imageName) );

    err2f = disp_error(flowGT(:,:,1), -flow(:,:,1), 2);
    err3f = disp_error(flowGT(:,:,1), -flow(:,:,1), 3);
    err4f = disp_error(flowGT(:,:,1), -flow(:,:,1), 4);
    err5f = disp_error(flowGT(:,:,1), -flow(:,:,1), 5);

    err2fn = disp_error(flowGT_noc(:,:,1), -flow(:,:,1), 2);
    err3fn = disp_error(flowGT_noc(:,:,1), -flow(:,:,1), 3);
    err4fn = disp_error(flowGT_noc(:,:,1), -flow(:,:,1), 4);
    err5fn = disp_error(flowGT_noc(:,:,1), -flow(:,:,1), 5);
    
    epeErr  = getEndPointError(cat(3, -flow, ones(size(I1))), flowGT);
    epeErrN = getEndPointError(cat(3, -flow, ones(size(I1))), flowGT_noc);
  end

    
  kittiStr = sprintf('DispPix-occ 2/3/4/5 %.3f & %.3f & %.3f & %.3f\nFlowPix-occ 2/3/4/5 %.3f & %.3f & %.3f & %.3f\n', err2f, err3f, err4f, err5f, err2f, err3f, err4f, err5f);
  kittiStr = sprintf('%s\nDispPix-noc 2/3/4/5 %.3f & %.3f & %.3f & %.3f\nFlowPix-noc 2/3/4/5 %.3f & %.3f & %.3f & %.3f\n', kittiStr, err2fn, err3fn, err4fn, err5fn, err2fn, err3fn, err4fn, err5fn);
  kittiStr = sprintf('%s\n EPE %.3f & EPE(noc) %.3f\n', kittiStr, epeErr, epeErrN);
  
  if ~doTesting
%     fid = fopen(sprintf('%s/RESULTS_%03d_%s.txt', p.sFolder, p.imgNr, date), 'w', 'n');
%     if (fid ~= -1)
%       fwrite(fid, kittiStr, 'char');
%       fclose(fid);
%     end
    fprintf(kittiStr);
  end

%  fprintf(1, 'Elapsed total time: %f\n', totalTime);
