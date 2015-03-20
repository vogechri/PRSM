function processVideoChal( subFrames, folder )

%        for i=7:14
%         processVideoChal( i, sprintf('vid_ChalProd_v3_J20_full_next', i));
%        end

%      for i=5:14
%        processVideoChal( i, sprintf('vid_Final3', i));
%      end

%processVideoChal( [5,6,7,8,9,11,12,13,14], sprintf('vid_Final3'));

%    for i=5:14
%       processVideoChal( i, sprintf('vid_test_nextProp_redo', i));
%     end
%   processVideoChal( 6, sprintf('vid_ChalProd_v3_J5_full_smo5_next', 6));

if ~isdeployed
  path(path,'./io/');
  path(path,'./io/other/');

  path(path,'./stuff/');
  path(path,'./mex/');
  path(path,'../../../export_fig');%plotting
  path(path, './stuff/sc/');%plotting
end

doVideo = 0;

p.sFolder = folder;

p.frames=2;
% my storage of challenges:

iptsetpref('ImshowBorder','tight');
if doVideo
%  writerObj  = VideoWriter('./eccvHCI_2.avi','Motion JPEG AVI');
  writerObj  = VideoWriter('./eccvHCI_2.mp4','MPEG-4');  
%  set(gca,'nextplot','replacechildren');
%  set(gcf,'Renderer','zbuffer');% ?
  writerObj.FrameRate = 10;
  open(writerObj);
end

if ~exist(p.sFolder,'dir')
  mkdir(p.sFolder);
end


for si = 1:numel(subFrames)
  
  subFrame = subFrames(si);
  
  p.subImg = subFrame;
  if p.subImg == 5
    io_folder = '../../data/ChallengingSequences/RainBlur/sequence';
    testImages = [3:27];% just 13 sucks -> call it a failure , use full sequence
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
  
%  io_folder = 'C:\Users\vogechri\Desktop\work\Journal\ECCV\solvePWRSMulti\SeriesDiss';
  io_folder = 'C:\Users\vogechri\Desktop\relatedThesis\gitThesis\christoph-thesis\flagPics\SeriesDiss_ps10GT_RNT';
  io_folder = 'C:\Users\vogechri\Desktop\work\Journal\ECCV\solvePWRSMulti\PureGTx';
  io_folder = 'C:\Users\vogechri\Desktop\relatedThesis\gitThesis\christoph-thesis\flagPics\JointFull512_check';
  io_folder = 'C:\Users\vogechri\Desktop\work\Source\SceneFlowICCV11\results\NewD\th125_l20_s16_o200_d1_i2_CSAD_1stD2F3J_pcgSq111_cap115_newData_med_full_3datas0X';

  testImages = [110,123,74];%[2,3,4,6,8,9];

  p.subImg = 10;
  
  for testImg_ = 1:numel(testImages)
    
    p.imgNr = testImages(testImg_);
    
    ending = sprintf('%02d_%06d', p.subImg, p.imgNr);
    
    p.testing =0;
%    [ImgL, flowGT, R_l, data_cam, imageName, ~, ~, ~, ~, ~, ~, flow2DGt, flow2DGt_noc] = loadKittiFlow('../../data' , imgNr, p);    

    p.subImg = 10;
    [cam, ref, imageName, flow2DGt, flow2DGt_noc] = loadKittiFlow('../../../data', p.imgNr, p);
    p.subImg = 0;
%    [cam, ref, imageName, flow2DGt, flow2DGt_noc] = loadTestFlow( p.imgNr, p);
%    [cam, ref, imageName, flow2DGt, flow2DGt_noc] = loadChallangingFlow(io_folder, p.imgNr, p);
    
    I_key1 = ref.I(1).I;
    I_key1 = uint8(I_key1*255);
    I_key3 = repmat(I_key1, [1,1,3]);
    
    % pic 8: bad flow
    
    %000004_05
    %  disp   = -disp_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/ChalPros_v3_J20_full_mm120_fromIt1/disp/%06d_%02d.png', p.imgNr, p.subImg));
    %  flow2d =  flow_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/ChalPros_v3_J20_full_mm120_fromIt1/flow/%06d_%02d.png', p.imgNr, p.subImg));
    %  disp   = -disp_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/chalResj20/disp/%06d_%02d.png', p.imgNr, p.subImg));
    %  flow2d =  flow_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/chalResj20/flow/%06d_%02d.png', p.imgNr, p.subImg));
    
    % flow BROKEN:
    %  disp   = -disp_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/FebChal_Rough_CrossingCars55/disp/%06d_%02d.png', p.imgNr, p.subImg));
    %  flow2d =  flow_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/FebChal_Rough_CrossingCars55/flow/%06d_%02d.png', p.imgNr, p.subImg));
    %  disp   = -disp_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/FebChal_Rough_WetAutobahnR/disp/%06d_%02d.png', p.imgNr, p.subImg));
    %  flow2d =  flow_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/FebChal_Rough_WetAutobahnR/flow/%06d_%02d.png', p.imgNr, p.subImg));
    
    %  disp   = -disp_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/FebChal_v4_J5/disp/%06d_%02d.png', p.imgNr, p.subImg));
    %  flow2d =  flow_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/FebChal_v4_J5/flow/%06d_%02d.png', p.imgNr, p.subImg));
    
    %  disp   = -disp_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/ChalProd_v3_J20_full_mm120_fromX14/disp/%06d_%02d.png', p.imgNr, p.subImg));
    %  flow2d =  flow_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/ChalProd_v3_J20_full_mm120_fromX14/flow/%06d_%02d.png', p.imgNr, p.subImg));
    
    % which is better ?
    %  disp   = -disp_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/ChalProd_v3_J5_full_smo5_next/disp/%06d_%02d.png', p.imgNr, p.subImg));
    %  flow2d =  flow_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/ChalProd_v3_J5_full_smo5_next/flow/%06d_%02d.png', p.imgNr, p.subImg));
    
    %  disp   = -disp_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/ChalProd_v3_J20_full_next/disp/%06d_%02d.png', p.imgNr, p.subImg));
    %  flow2d =  flow_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/ChalProd_v3_J20_full_next/flow/%06d_%02d.png', p.imgNr, p.subImg));
    
    %  disp   = -disp_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/ChalProd_v3_J20_smo35_next/disp/%06d_%02d.png', p.imgNr, p.subImg));
    %  flow2d =  flow_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/ChalProd_v3_J20_smo35_next/flow/%06d_%02d.png', p.imgNr, p.subImg));
    
    
    
    
    %  disp   = -disp_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/ChalProd_J20_nextProp/disp/%06d_%02d.png', p.imgNr, p.subImg));
    %  flow2d =  flow_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/ChalProd_J20_nextProp/flow/%06d_%02d.png', p.imgNr, p.subImg));
    
    % actually not using the 2nd prop - which worked for crossing above - but
    % due to ego stuff sucked at cartruck
    %  disp   = -disp_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/ChalProd_J20_nextProp_redo/disp/%06d_%02d.png', p.imgNr, p.subImg));
    %  flow2d =  flow_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/ChalProd_J20_nextProp_redo/flow/%06d_%02d.png', p.imgNr, p.subImg));
 
    
%    disp   = -disp_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/final/disp/%06d_%02d.png', p.imgNr, p.subImg));
%    flow2d =  flow_read(sprintf('C:/Users/vogechri/Desktop/work/init/Challenges/final/flow/%06d_%02d.png', p.imgNr, p.subImg));


    disp   = -disp_read(sprintf('%s/disp/%06d_%02d.png', io_folder , p.imgNr, p.subImg));
    flow2d =  flow_read(sprintf('%s/flow/%06d_%02d.png', io_folder , p.imgNr, p.subImg));    
    
    if ~doVideo
      
      uDisp = flow2d(:,:,1);
      vDisp = flow2d(:,:,2);
      fNr = 1;
      disp = max(-65, disp);
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
      max_flow = 60;
      max_flow = 85; % kitti
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
      
      avi_Frame(ref.I(1).I,cam.I(1).I, disp, flow2d, 1);
      export_fig( sprintf('%s/%s_joint.png', p.sFolder, ending), '-m1');
      close force all
      
      
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
  function frame = avi_Frame(I1,I2, disp, flowL_2d, video)

    maxDisp = 85;
    maxVecSize = 65;
    frameSize = size(I1); % kitti : cams rectification leads to diff in size
    frameSize (2) = frameSize (2) *2; % two pics on top
        
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
    I_cam3 =  repmat(uint8(I_cam*255), [1,1,3]);
    
    fprintf(sprintf('processing:  I:/videoResult/flow/%06d_%02d.png\n', p.imgNr, p.subImg));
    
    %  flowL_2d   = flow_read_Kitti(sprintf('I:/videoResult/flow/%06d_%02d.png', p.imgNr, p.subImg));
    %  disp       = disp_read_Kitti(sprintf('I:/videoResult/disp/%06d_%02d.png', p.imgNr, p.subImg));
    
    %  f=figure(9),imshow(-disp, 'Border','tight'), colormap(gray);


    tmp = sc((-disp), [0, (maxDisp)], 'hicontrast');%contrast    
    
%   tmp = sc(-disp, 'hicontrast');%contrast
    tmp = sc(sqrt(-disp), [0, sqrt(maxDisp)], 'hicontrast');%contrast
%     tmp = sc(-disp, 'hicontrast');imshow(tmp, 'Border','tight'), colorbar
%     tmp = sc(-disp, [12, 20], 'hicontrast');imshow(tmp, 'Border','tight'), colorbar
     
    %  tmp = getframe(f);tmp=tmp.cdata;
    tmp = uint8(255*tmp);
    tmp = 0.5*tmp+0.5*I_key3;
    %  close(f);
    f = figure(1);set(f, 'visible','on');
    %  imshow(tmp), axis image, axis off;
    
    % with flow 2d
    %  flowL_2d(:,:,1:2) = maxVecSize*bsxfun( @rdivide, flowL_2d(:,:,1:2), max (maxVecSize, sqrt(flowL_2d(:,:,1).^2 + flowL_2d(:,:,2).^2)) );
    
    length = sqrt(flowL_2d(:,:,1).^2+flowL_2d(:,:,2).^2);
    % K
    flowL_2d(:,:,1:2) = bsxfun( @times, sqrt(sqrt(min(maxVecSize, length))), bsxfun(@rdivide, flowL_2d(:,:,1:2), length));    
    % flag
%    flowL_2d(:,:,1:2) = bsxfun( @times, (sqrt(min(maxVecSize, length))), bsxfun(@rdivide, flowL_2d(:,:,1:2), length));
    
    imshow(tmp, 'Border','tight'),colormap('default'), hold on, plotflow2(flowL_2d(:,:,1:2), 'vector',0);
    %  imshow(repmat(partImg1, [1,1,3])),colormap('default'), hold on, plotflow2(flowPart, 'vector',0, scale);

    if exist('video','var')    
      return;
    end
    tmp = getframe(f);tmp=tmp.cdata;
    I_topPic = imresize( cat(2, I_key1, uint8(I_cam*255) ), [size(tmp,1), size(tmp,2)] );
%    I_topPic = imresize( cat(2, I_key1, uint8(I_cam*255) ), [size(tmp,1), size(tmp,2)] );
    tmp = cat(2, repmat(I_topPic, [1,1,3]), tmp);
    close(f);
    
    videoTmp(:,:,1) = imresize( tmp(:,:,1), frameSize );
    videoTmp(:,:,2) = imresize( tmp(:,:,2), frameSize );
    videoTmp(:,:,3) = imresize( tmp(:,:,3), frameSize );
    
    h = figure(2);imshow(videoTmp, 'Border','tight');
    
    if exist('video','var')
      frame = getframe;
      %addframe
      close(h);
    else      
    end

  end
end