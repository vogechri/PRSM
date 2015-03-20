function processFrameKitti( subFrames, folder, testImages, io_folder )

% with error masks etc?
kittipics =1;

%processFrameKitti(10, './JournalPlots', 74, 'C:\Users\vogechri\Desktop\work\init\Journal\ECCV\figs\figs\AllewOn_auto0_halfEvalEnergy_refx0_3f_lsaux_s04_S1J1E1_sw01_as07_tp025GH_3fps_pwds')

% journal plot
%processFrameKitti(10, './JournalPlots', 74, 'C:\Users\vogechri\Desktop\work\init\Journal\ECCV\figs\figs\AllewOn_auto0_halfEvalEnergy_refx1_3f_lsaux_s04_S1J1E1_sw01_as07_tp025GH_3fps_pwds')


% example
%processFrameKitti(10, './iccv_occ', 171, 'C:\Users\vogechri\Desktop\work\ICCV\results\PWRSFCode\ICCV_Centred_ego2nd_smo31_occ_noTrunc_vis')

%        for i=7:14
%         processVideoChal( i, sprintf('vid_ChalProd_v3_J20_full_next', i));
%        end

%      for i=5:14
%        processVideoChal( i, sprintf('vid_Final2', i));
%      end


%    for i=5:14
%       processVideoChal( i, sprintf('vid_test_nextProp_redo', i));
%     end
%   processVideoChal( 6, sprintf('vid_ChalProd_v3_J5_full_smo5_next', 6));

if ~isdeployed
  path(path,'./io/');
  path(path,'./stuff/');
  path(path,'./mex/');
  path(path,'../../export_fig');%plotting
  path(path, './stuff/sc/');%plotting
end

doVideo = 0;
%max_flow = 300;
%max_disp = 200;

max_flow = 65;
max_disp = 65;

p.sFolder = folder;
p.testing = 0;
p.frames=2;
% my storage of challenges:

iptsetpref('ImshowBorder','tight');
if doVideo
  writerObj  = VideoWriter('./eccvHCI_2.avi','Motion JPEG AVI');
  writerObj.FrameRate = 10;
  open(writerObj);
end

if ~exist(p.sFolder,'dir')
  mkdir(p.sFolder);
end


for si = 1:numel(subFrames)
  
  subFrame = subFrames(si);
  p.subImg = subFrame;  
  %{
  p.subImg = subFrame;
  if p.subImg == 5
    io_folder = '../../data/ChallengingSequences/RainBlur/sequence';
    testImages = [15:27];
  end
  if p.subImg == 6
    io_folder = '../../data/ChallengingSequences/BlinkingArrow/sequence';
    testImages = [3:28];
  end
  if p.subImg == 7
    io_folder = '../../data/ChallengingSequences/RainFlares/sequence';
    testImages = [3:28];
  end
  if p.subImg == 8
    io_folder = '../../data/ChallengingSequences/CarTruck/sequence';
    testImages = [9:28];
  end
  if p.subImg == 9
    io_folder = '../../data/ChallengingSequences/NightAndSnow/sequence';   % 0:64
    testImages = [3:20];
  end
  if p.subImg == 10
    io_folder = '../../data/ChallengingSequences/SunFlare/sequence';       % 0:30
    testImages = [];
  end
  if p.subImg == 11
    io_folder = '../../data/ChallengingSequences/ShadowOnTruck/sequence';  % 0:30
    testImages = [3:28];
  end
  if p.subImg == 12
    io_folder = '../../data/ChallengingSequences/CrossingCars/sequence';    % 0:30
    testImages = [5:28];
  end
  if p.subImg == 13
    io_folder = '../../data/ChallengingSequences/FlyingSnow/sequence';   % 0:30
    testImages = [3:28];
  end
  if p.subImg == 14
    io_folder = '../../data/ChallengingSequences/WetAutobahn/sequence';   % 0:30
    testImages = [12:28];
  end
  %}
  for testImg_ = 1:numel(testImages)
    
    p.imgNr = testImages(testImg_);
    
    ending = sprintf('%02d_%06d', p.subImg, p.imgNr);
    
    if kittipics
    [cam, ref, imageName, flow2DGt, flow2DGt_noc] = loadKittiFlow('../../../data/' , p.imgNr, p);
    else
if p.subImg == 5
  pic_folder = '../../data/ChallengingSequences/RainBlur/sequence';
end
if p.subImg == 6
  pic_folder = '../../data/ChallengingSequences/BlinkingArrow/sequence';
end
if p.subImg == 7
  pic_folder = '../../data/ChallengingSequences/RainFlares/sequence';
end
if p.subImg == 8
  pic_folder = '../../data/ChallengingSequences/CarTruck/sequence';
end

% my storage of challenges:
if p.subImg == 9
  pic_folder = '../../data/ChallengingSequences/NightAndSnow/sequence';   % 0:64
end
if p.subImg == 10
  pic_folder = '../../data/ChallengingSequences/SunFlare/sequence';       % 0:30
end
if p.subImg == 11
 pic_folder = '../../data/ChallengingSequences/ShadowOnTruck/sequence';  % 0:30
end
if p.subImg == 12
 pic_folder = '../../data/ChallengingSequences/CrossingCars/sequence';    % 0:30
end
if p.subImg == 13
 pic_folder = '../../data/ChallengingSequences/FlyingSnow/sequence';   % 0:30
end
if p.subImg == 14
 pic_folder = '../../data/ChallengingSequences/WetAutobahn/sequence';   % 0:30
end
    [cam, ref, imageName, flow2DGt, flow2DGt_noc] = loadChallangingFlow(pic_folder, p.imgNr, p);
    end
    
    
    I_key1 = ref.I(1).I;
    I_key1 = uint8(I_key1*255);
    I_key3 = repmat(I_key1, [1,1,3]);
    
p.subImg = 0;    
    disp   = -disp_read(sprintf('%s/disp/%06d_%02d.png', io_folder, p.imgNr, p.subImg));
    flow2d =  flow_read(sprintf('%s/flow/%06d_%02d.png', io_folder, p.imgNr, p.subImg));
p.subImg =10;
    
    if ~doVideo
      
      uDisp = flow2d(:,:,1);
      vDisp = flow2d(:,:,2);
      fNr = 1;
      disp = max(-max_disp, disp);
      %-% disparity image
      tmp = max( 0, -disp);
      f = figure(fNr+testImg_);set(f, 'visible','off');
      sc(tmp, 'hicontrast');
      tmp = getframe(f);tmp=tmp.cdata;
      tmp = 0.5*tmp+0.5*I_key3;
      close(f);f = figure(fNr+testImg_);set(f, 'visible','off');
      imshow(tmp), axis image, axis off;
      export_fig( sprintf('%s/%s_disp_HQ.png', p.sFolder, ending), '-m1');
      
      %  max_flow = max([abs(uDisp(:)); abs(vDisp(:))]);

      F_mag = sqrt(uDisp.*uDisp+vDisp.*vDisp);
      F_dir = atan2(vDisp,uDisp);
      I_flow2d = flow_map(F_mag(:),F_dir(:),ones(numel(uDisp),1),max_flow,8);% 8 is normal convention
      I_flow2d = reshape(I_flow2d, [size(uDisp, 1), size(uDisp, 2), 3]);
      %    I_flow2d = flowToColor(cat(3, uDisp, vDisp))./255;
      
      %    I_flow2d = flow_map(F_mag(:),F_dir(:),ones(1,numel(uDisp)),max_flow,10*max_flow);I_flow2d = reshape(I_flow2d, [size(uDisp, 1), size(uDisp, 2), 3]);figure(9), imshow(I_flow2d)
      
      
      if exist('occIds2', 'var' )
        occludedArea = false(size(Seg.Img));
        if numel(occIds2) == numel(Seg.Img)
          occludedArea( occIds2 ~=0 ) = true;
        else
          for i=1:numel(occIds2)
            occludedArea (Seg.Img == occIds2(i)) = true;
          end
        end
      end
      
      f = figure(fNr+18);set(f, 'visible','off');
      imshow(I_flow2d, 'Border','tight'),axis image,axis off;
      tmp = getframe(f);tmp=tmp.cdata;
      
      if exist('occIds2', 'var' )
        tmp( repmat( occludedArea, [1,1,3]) ) = 0;% white, black also ok
      end
      
      tmp = 0.5*tmp+0.5*I_key3;
      close(f);f = figure(fNr+18);set(f, 'visible','off');
      imshow(tmp), axis image, axis off;
      export_fig( sprintf('%s/%s_flow.png', p.sFolder, ending), '-m1');
      close(f);close force all
      %  plotAnalysis(ref, cam, N_proj(:,ids_prevFrame), Rt_proj(:,:,ids_prevFrame), Seg, u, 20, par, sprintf('%03d_projectedPrevSolution', par.imgNr));
      
      tmp=avi_Frame(ref.I(1).I,cam.I(1).I, disp, flow2d, ref.I(2).I,cam.I(2).I, 1);
      % abused:
      %{
%      tmp = getframe(f);
      tmp=tmp.cdata;
      imshow(tmp), axis image, axis off;
      export_fig( sprintf('%s/%s_DispFlow.png', p.sFolder, ending), '-m1');
%}
      
      %%%%%%%%%%%% try to add some info on error:

      addpath(' C:\Users\vogechri\Desktop\work\ofdata\GCPR\devkit_stereo_flow\devkit\matlab\specialSession');
      %load GT:
%      flow2DGt, flow2DGt_noc
    
      F_err = flow_error_image(flow2DGt(:,:,2:4),flow2d);
%      figure,imshow([flow_to_color([flow2DGt(:,:,2:4);flow2d]);F_err]);
%      tmpall = cat(2,tmp,  F_err*255 );
%{
      F_err_v = flow2d(:,1:end-1,1:2)-flow2d(:,2:end,1:2);
      F_err_h = flow2d(1:end-1,:,1:2)-flow2d(2:end,:,1:2);
      F_err_v = round(sqrt(F_err_v(:,:,1).*F_err_v(:,:,1)+F_err_v(:,:,2).*F_err_v(:,:,2)));
      F_err_h = round(sqrt(F_err_h(:,:,1).*F_err_h(:,:,1)+F_err_h(:,:,2).*F_err_h(:,:,2)));

      figure(1),hist(F_err_v(:), [0:max(F_err_v(:))]);axis tight;axis off
      h = findobj(gca,'Type','patch');
      set(h,'FaceColor','g','EdgeColor','w');
      export_fig( sprintf('%s/%s_ver.png', p.sFolder, ending), '-m1');      
      figure(2),hist(F_err_h(:), [0:max(F_err_h(:))]);axis tight;axis off
      h = findobj(gca,'Type','patch');
      set(h,'FaceColor','r','EdgeColor','w')
      export_fig( sprintf('%s/%s_hor.png', p.sFolder, ending), '-m1');
  %}    
      F_errD = disp_error_image(flow2DGt(:,:,1), -disp);
      tmpall = cat(1,tmp,  cat(2,F_errD*255, F_err*255 ));
      imwrite(tmpall, sprintf('%s/%s_DispFlowErr.png', p.sFolder, ending));

      imwrite(tmp, sprintf('%s/%s_DispFlow.png', p.sFolder, ending));
      
      %  h = figure(2);imshow(videoTmp, 'Border','tight');
      %  frame = getframe;
      %addframe
      %  close(h);
    else
      frame=avi_Frame(ref.I(1).I,cam.I(1).I, disp, flow2d);
      writeVideo(writerObj,frame);
    end
  end
end
if doVideo
  close(writerObj);
end
% img=45;
% subset = [45,116,132,138,147,172 , 2,3,5,6,7,8, 13,15];
%
% addpath('C:\Users\vogechri\Desktop\work\OpticalFlowModified');
% %addpath('C:\Sceneflow\OpticalFlowModified');
%
% for i = 1:numel(subset)
%   img = subset(i);
%   for subI = 0:19
%     frame = processStereo( 0.4, 1, 1, 4, 3, 0.9, sprintf('C:/Users/vogechri/Desktop/work/data/postICCV/OpenLab/final/', img), 'Test.inf', img, 0, 0, 99, 500, 1, 0, 0.5, 0, 0, -1, 0, subI);
% %    frame = getframe;
%     writeVideo(writerObj,frame);
%   end
% end
%
% rmpath('C:\Users\vogechri\Desktop\work\OpticalFlowModified');
% %rmpath('C:\Sceneflow\OpticalFlowModified');
% close(writerObj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function I = flow_map(F_mag,F_dir,F_val,max_flow,n)
    
    I(:,1) = mod(F_dir/(2*pi),1);
    I(:,2) = F_mag * n / max_flow;
    I(:,3) = n - I(:,2);
    %I(:,3) = F_mag * n / max_flow;
    I(:,[2 3]) = min(max(I(:,[2 3]),0),1);
    I = hsv2rgb(reshape(I,[],1,3));
    
    %I(:,:,1) = I(:,:,1).* F_val;
    %I(:,:,2) = I(:,:,2).* F_val;
    %I(:,:,3) = I(:,:,3).* F_val;
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function frame = avi_Frame(I1,I2, disp, flowL_2d, I1t,I2t, useForPlotting)
    
    maxDisp = 65;
%    maxDisp = 120;
    maxVecSize = 65;
    frameSize = size(I1); % kitti : cams rectification leads to diff in size
    frameSize (2) = frameSize (2) *2; % two pics on top
    if exist('I1t', 'var')
%    frameSize (1) = frameSize (1) *2; % two pics on top
%      I_topPic = imresize( cat(1, cat(2, I_key1, uint8(I_cam*255) ), cat(2, uint8(I1t*255), uint8(I2t*255) )), [size(tmp,1), size(tmp,2)] );
    end
    
    if ~exist ('useForPlotting' , 'var')
      useForPlotting = false;
    end
    
    I_cam = I2;
    I_key = I1;
    if size(I_key,3) > 1
      I_cam  = mean(I_cam,3);
      I_key1 = mean(I_key,3);
    else
      I_key1 = I_key;
    end
    I_key1 = uint8(I_key1*255);
    I_key3 = repmat(I_key1, [1,1,3]);
    I_cam3 = repmat(uint8(I_cam*255), [1,1,3]);
    

    fprintf(sprintf('processing:  I:/videoResult/flow/%06d_%02d.png\n', p.imgNr, p.subImg));
    
    %  flowL_2d   = flow_read_Kitti(sprintf('I:/videoResult/flow/%06d_%02d.png', p.imgNr, p.subImg));
    %  disp       = disp_read_Kitti(sprintf('I:/videoResult/disp/%06d_%02d.png', p.imgNr, p.subImg));
    
    %  f=figure(9),imshow(-disp, 'Border','tight'), colormap(gray);
    
    %  tmp = sc(-disp, 'hicontrast');%contrast
    tmp = sc(sqrt(-disp), [0, sqrt(maxDisp)], 'hicontrast');%contrast
    
    %  tmp = getframe(f);tmp=tmp.cdata;
    tmp = uint8(255*tmp);
    tmp = 0.5*tmp+0.5*I_key3;
    %  close(f);
    f = figure(1);set(f, 'visible','on');
    %  imshow(tmp), axis image, axis off;
    
    % with flow 2d
    %  flowL_2d(:,:,1:2) = maxVecSize*bsxfun( @rdivide, flowL_2d(:,:,1:2), max (maxVecSize, sqrt(flowL_2d(:,:,1).^2 + flowL_2d(:,:,2).^2)) );
    
    length = sqrt(flowL_2d(:,:,1).^2+flowL_2d(:,:,2).^2);
    flowL_2d(:,:,1:2) = bsxfun( @times, sqrt(sqrt(min(maxVecSize, length))), bsxfun(@rdivide, flowL_2d(:,:,1:2), length));
    
    
    imshow(tmp, 'Border','tight'),colormap('default'), hold on, plotflow2(flowL_2d(:,:,1:2), 'vector',0);
    %  imshow(repmat(partImg1, [1,1,3])),colormap('default'), hold on, plotflow2(flowPart, 'vector',0, scale);
    
    tmp = getframe(f);tmp=tmp.cdata;
    I_topPic = imresize( cat(2, I_key1, uint8(I_cam*255) ), [size(tmp,1), size(tmp,2)] );
%    I_topPic = imresize( cat(2, I_key1, uint8(I_cam*255) ), [size(tmp,1), size(tmp,2)] );

    if exist('I1t', 'var')
      I_topPic = imresize( cat(1, cat(2, I_key1, uint8(I_cam*255) ), cat(2, uint8(I1t*255), uint8(I2t*255) )), [size(tmp,1), size(tmp,2)] );
    end

    tmp = cat(2, repmat(I_topPic, [1,1,3]), tmp);
    close(f);
    
    videoTmp(:,:,1) = imresize( tmp(:,:,1), frameSize );
    videoTmp(:,:,2) = imresize( tmp(:,:,2), frameSize );
    videoTmp(:,:,3) = imresize( tmp(:,:,3), frameSize );
    
    if useForPlotting
      frame = videoTmp;
    else
      h = figure(2);imshow(videoTmp, 'Border','tight');
      frame = getframe;
      %addframe
      close(h);
    end    

  end
end