function [cam, ref, imageName, flow2DGt, flow2DGt_noc] = loadTestFlow( testImg, p)
%[ImgL, flowGT, R_l, data_supp, imageName, strVz, p, stereoD, validMap, IltReal, Ilt1Real] = ...
%


% load images
  % images should be in [0, 1]
  %p.frames = 0; % not more possible here
  validMap = 0;
  stereoD = 0;
%   data_supp = struct( 'It0', {}, 'It1', {}, 'R', {}, 'F', {}, 'Ft', {}, ...
%     'epi', {}, 'pop', {}, 'popRef', {}, 'Rot', {}, 'Tra' , {} );

%  data_supp = struct( 'I', {}, 'R', {}, 'F', {}, 'Ft', {}, 'epi', {}, ...
%    'pop', {}, 'popRef', {}, 'Rot', {}, 'Tra' , {});
  data_supp = struct( 'I', {}, 'R', {}, 'F', {}, 'Ft', {}, 'epi', {}, ...
    'pop', {}, 'popRef', {},'Kl', {}, 'Kr', {}, 'Rl', {}, 'Rr', {}, ...
    'Tl', {}, 'Tr', {});

  
  ImgL      = {};
  ImgR      = {};
  
  imageName = 'Test';
  strVz = '../../data/SceneFlow-data/box14T3';
  
  folder = '../../../data/SceneFlow-data/';
  %{
  if testImg == 1
    strVz = sprintf( '%s/box16T0', folder);
    imageName = 'Teddy';
  end

  if testImg == 2
    strVz = sprintf( '%s/box16T1', folder);
    imageName = 'Whale';
  end

  if testImg == 3
      strVz = sprintf( '%s/box16T2', folder);
    imageName = 'Leaf';
  end

  if testImg == 4
    strVz = sprintf( '%s/box16T3', folder);
    imageName = 'Cars';
  end
  
  if testImg == 5
    strVz = sprintf( '%s/box10', folder);
    imageName = 'No_Moving_box10';
  end
  
  if testImg == 6
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box13';
    imageName = 'Zoom_box13';
  end

  if testImg == 7
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box11';
    imageName = 'Overlap_box11';
  end

  if testImg == 8
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box12New';
    imageName = 'Overlap_box12New';
  end

  if testImg == 9
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box17_unsteadyFlow';
    imageName = 'UnsteadyFlow_box17';
  end
  %}
  if testImg == 10
    strVz = '../../data/SceneFlow-data/box19_small_unsteady';
    imageName = 'box19_small_unsteady';
  end
  
  if testImg == 11
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box20_FrontoParallel';
    imageName = 'box20_FrontoParallel';
  end

  if testImg == 12
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box21_NearFronto';
    imageName = 'box21_NearFronto';
  end
 
  if testImg == 13
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box22_NearFronto2';
    imageName = 'box22_NearFronto2';
  end
  
  if testImg == 14
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box23_smallFronto';
    imageName = 'box23_smallFronto';
  end
  
  if testImg == 15
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box24_smallFronto';
    imageName = 'box24_smallFronto';
  end
  
  if testImg == 16
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box25_smallFronto';
    imageName = 'box25_smallFronto';
  end
  
  if testImg == 17
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box26_tiny';
    imageName = 'box26_tiny';
  end
  
  if testImg == 18
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box27_tiny_near';
    imageName = 'box27tiny_near';
  end
  
  if testImg == 19
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box28_tiny_near_dy';
    imageName = 'box28tiny_near_dy';
  end
  
  if testImg == 20
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box29_tiny_near_rot';
    imageName = 'box29tiny_near_rot';
  end

  if testImg == 21
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box30_near';
    imageName = 'box30_near';
  end
  
  if testImg == 22
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box31_near';
    imageName = 'box31_near';
  end

  % 32 and 33 are 38 background 6.5 foreground
  if testImg == 23
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box_32_tiny_noMove';
    imageName = 'box32-tiny-noMove';
  end

  if testImg == 24
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box_33_tiny_noMove2';
    imageName = 'box33-tiny-noMove2';
  end

  if testImg == 25
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box_34_tiny_moveXY';
    imageName = 'box_34_tiny_moveXY';
  end

  if testImg == 26
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box_35_tiny_moveXY2';
    imageName = 'box_35_tiny_moveXY2';
  end

  if testImg == 27 % movement is small compared to box_37
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box_36_tiny_noMove_largeDepth';
    imageName = 'box_36_tiny_noMove_largeDepth';
  end
  
  % box_37 is 72-90 background 29 foreground
  if testImg == 28
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box_37_tiny_noMove_largeDepth2';
    imageName = 'box37-tiny-noMove-largeDepth2';
  end

    % 38 is 75-80 background 30 foreground
  if testImg == 29
    strVz = '/home/christop/Matlab/data/SceneFlow-data/box_38_tiny_noMove_largeDepth3';
    imageName = 'box_38_tiny_noMove_largeDepth3';
  end

  % sphere1 is 72-90 background 33 foreground
  if testImg == 30
    strVz = '/home/christop/Matlab/data/SceneFlow-data/sphere_1_Test';
    imageName = 'sphere1';
  end
  
  % sphere2 is 155 background 65 foreground
  if testImg == 31
%    strVz = '/home/christop/Matlab/data/SceneFlow-data/sphere_2_test';
%    imageName = 'sphere2';
    strVz = sprintf( '%s/highResTest', folder);
    imageName = 'highResTest';
  end
  % depth 162:bg 80:fg, needs higher init d_ = 50
  if testImg == 32
    strVz = sprintf( '%s/box_40_noMove_', folder);
    imageName = 'box40noMove';
  end
  % right x - 2
  if testImg == 33
    strVz = sprintf( '%s/box_40_movex_', folder);
    imageName = 'box40MoveX';
  end
  % y-left - 2.4
  if testImg == 34
    strVz =sprintf( '%s/box_40_movey_', folder);
    imageName = 'box40MoveY';
  end
  % z-left - -9.6
  if testImg == 35
    strVz = sprintf( '%s/box_40_movez_', folder);
    imageName = 'box40MoveZ';
  end
  
  % depth bg 80 near 40
  if testImg == 36
    strVz = sprintf( '%s/box_41_noMove', folder);
    imageName = 'box41noMove';
  end

  if testImg == 37
    strVz = sprintf( '%s/box_41_movex', folder);
    imageName = 'box41MoveX';
  end
    
  if testImg == 38
    strVz = sprintf( '%s/box_41_movey', folder);
    imageName = 'box41MoveY';
  end
        
  if testImg == 39
    strVz = sprintf( '%s/box_41_movez', folder);
    imageName = 'box41MoveZ';
  end
  % depth bg 40 near 17
  if testImg == 40
    strVz = sprintf( '%s/box_42_near_noMove', folder);
    imageName = 'box42noMove';
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  if testImg == 41
    strVz = sprintf( '%s/Combi1', folder);
    imageName = 'Combi1';
%    strVz = sprintf( '%s/box_42_near_MoveX', folder);
%    imageName = 'box42MoveX';
  end
    
  if testImg == 42
    strVz = sprintf( '%s/Combi2', folder);
    imageName = 'Combi2';
%    strVz = sprintf( '%s/box_42_near_MoveY', folder);
%    imageName = 'box42MoveY';
  end
        
  if testImg == 43
    strVz = sprintf( '%s/Rot5', folder);
    imageName = 'Rot5';
%    strVz = sprintf( '%s/box_42_near_MoveZ', folder);
%    imageName = 'box42MoveZ';
  end
  
  % depth bg 80 near 40
  if testImg == 44
    strVz = sprintf( '%s/Rot6', folder);
    imageName = 'Rot6';
%    strVz = sprintf( '%s/box_43_near_noMove', folder);
%    imageName = 'box40noMove';
  end
  % right x - 2
  if testImg == 45
    strVz = sprintf( '%s/Zoom3', folder);
    imageName = 'Zoom3';
%    strVz = sprintf( '%s/box_43_near_movex', folder);
%    imageName = 'box40MoveX';
  end
  % y-left - 2.4
  if testImg == 46
    strVz = sprintf( '%s/Zoom9', folder);
    imageName = 'Zoom9';
%    strVz = sprintf( '%s/box_43_near_movey', folder);
%    imageName = 'box40MoveY';
  end
  % z-left - -9.6
  if testImg == 47
%    strVz = sprintf( '%s/box_43_near_movez', folder);
%    imageName = 'box40MoveZ';
    strVz = sprintf( '%s/SphereRot2', folder);
    imageName = '3DSphere2';
  end
  if testImg == 48
    %strVz = sprintf( '%s/test3DSphere', folder);
    %imageName = 'test3DSphere';
    strVz = sprintf( '%s/SphereRot', folder);
    imageName = '3DSphere';
  end
  if testImg == 49
    strVz = sprintf( '%s/highResTest', folder);
    imageName = 'highResTest';
  end 
  
  if testImg == 50
    strVz = sprintf( '%s/box_highRes_dx', folder);
    imageName = 'highRes_dx';
  end
  
  if testImg == 51
    strVz = sprintf( '%s/box_highRes_dy', folder);
    imageName = 'highRes_dy';
  end

  if testImg == 52
    strVz = sprintf( '%s/box_highRes_dr', folder);
    imageName = 'highRes_dr';
  end

  if testImg == 53
    strVz = sprintf( '%s/box_highRes_dr2', folder);
    imageName = 'highRes_dr2';
  end

  if testImg == 54
    strVz = sprintf( '%s/box_highRes_dr3', folder);
    imageName = 'highRes_dr3';
  end
  
  %%% 1-9 flags:
  if testImg == 1
    strVz = sprintf( '%s/flag1', folder);
    imageName = 'flag1';
  end
  if testImg == 2
    strVz = sprintf( '%s/flag2', folder);
    imageName = 'flag2';
  end
  if testImg == 3
    strVz = sprintf( '%s/flag3', folder);
    imageName = 'flag3';
    strVz = sprintf( '%s/flag1New/againA/', folder);
    imageName = 'flag1New/againA';
  end  
  if testImg == 4
    strVz = sprintf( '%s/flag4', folder);
    imageName = 'flag4';
  end

  if testImg == 5
    strVz = sprintf( '%s/flag5', folder);
    imageName = 'flag5';
  end  
  if testImg == 6
    strVz = sprintf( '%s/flag6', folder);
    imageName = 'flag6';
  end
  
  % ok
  if testImg == 7
    strVz = sprintf( '%s/flag7', folder);
    imageName = 'flag7';
  end 
  % now fine: large improve
  if testImg == 8
    strVz = sprintf( '%s/flag8', folder);
    imageName = 'flag8';
  end  
  % good: 
  if testImg == 9
    strVz = sprintf( '%s/flax1', folder);
    imageName = 'flag9';
  end
  %%%%%%%% Session 1 %%%%%%%%%%%
  if testImg == 100
    strVz = sprintf( '%s/RealWorld_1', folder);
    imageName = 'RealWorld_1';
  end

  if testImg == 101
    strVz = sprintf( '%s/RealWorld_2', folder);
    imageName = 'RealWorld_2';
  end
  
  if testImg == 102
    strVz = sprintf( '%s/RealWorld_3', folder);
    imageName = 'RealWorld_3';
  end
  
  if testImg == 103
    strVz = sprintf( '%s/RealWorld_4', folder);
    imageName = 'RealWorld_4';
  end
  %%%%%%%% Session 2 %%%%%%%%%%%
  if testImg == 104
    strVz = sprintf( '%s/Lego1', folder);
    imageName = 'Lego1';
  end

  if testImg == 105
    strVz = sprintf( '%s/Lego2', folder);
    imageName = 'Lego2';
  end
  
  if testImg == 106
    strVz = sprintf( '%s/Lego3', folder);
    imageName = 'Lego3';
  end
  
  if testImg == 107
    strVz = sprintf( '%s/Lego4', folder);
    imageName = 'Lego4';
  end
  
  if testImg == 108
    strVz = sprintf( '%s/Lego5', folder);
    imageName = 'Lego5';
  end
  
  if testImg == 109
    strVz = sprintf( '%s/Lego6', folder);
    imageName = 'Lego6';
  end

  if testImg == 110
    strVz = sprintf( '%s/Dose1', folder);
    imageName = 'Dose1';
  end

  if testImg == 111
    strVz = sprintf( '%s/Dose2', folder);
    imageName = 'Dose2';
  end

  if testImg == 112
    strVz = sprintf( '%s/Cat', folder);
    imageName = 'Cat';
  end
  
  if testImg == 113
    strVz = sprintf( '%s/CatMug', folder);
    imageName = 'CatMug';
  end

  if testImg == 114
    strVz = sprintf( '%s/Small2', folder);
    imageName = 'Small2';
  end

  if testImg == 115
    strVz = sprintf( '%s/Book1', folder);
    imageName = 'Book1';
  end
    
  if testImg == 116
    strVz = sprintf( '%s/Book2', folder);
    imageName = 'Book2';
  end

  if testImg == 117
    strVz = sprintf( '%s/MultiObject1', folder);
    imageName = 'MultiObject1';
  end

  if testImg == 118
    strVz = sprintf( '%s/MultiObject2', folder);
    imageName = 'MultiObject2';
  end
  %%%%%%%% Session 3 %%%%%%%%%%%
  if testImg == 119
    strVz = sprintf( '%s/LegoMove1', folder);
    imageName = 'LegoMove1';
  end
  
  if testImg == 120
    strVz = sprintf( '%s/LegoMove2', folder);
    imageName = 'LegoMove2';
  end
  
  if testImg == 121
    strVz = sprintf( '%s/LegoStraight', folder);
    imageName = 'LegoStraight';
  end
  
  if testImg == 122
    strVz = sprintf( '%s/LegoStraightTurn', folder);
    imageName = 'LegoStraightTurn';
  end
  
  if testImg == 123
    strVz = sprintf( '%s/LegoTurn1', folder);
    imageName = 'LegoTurn1';
  end
  
  if testImg == 124
    strVz = sprintf( '%s/LegoTurn2', folder);
    imageName = 'LegoTurn2';
  end
  
  if testImg == 125
    strVz = sprintf( '%s/LegoTurn3', folder);
    imageName = 'LegoTurn3';
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Basha images : 
  if testImg == 55
    strVz = sprintf( '%s/CarpetBall_Full', folder);
    imageName = 'CarpetBall';
  end
  if testImg == 56
    strVz = sprintf( '%s/cars1Basha', folder);
    imageName = 'cars1Basha';
  end
  if testImg == 57
    strVz = sprintf( '%s/catBasha', folder);
    imageName = 'CatBasha';
  end
  if testImg == 58
    strVz = sprintf( '%s/mariaBasha', folder);
    imageName = 'Maria';
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  if testImg > 200 && testImg <= 300
    strVz = sprintf( '%s/Erdbeer', folder);
    imageName = sprintf('Erdbeer_%03d', testImg-200);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  if testImg == 301
    strVz = sprintf( '%s/HoguetBall/Ball', folder);
    imageName = 'HoguetBall';
  end
  if testImg == 302
    strVz = sprintf( '%s/HoguetBall/Venus', folder);
    imageName = 'Venus';
  end
  if testImg == 303
    strVz = sprintf( '%s/HoguetBall/Teddy', folder);
    imageName = 'Teddy';
  end
  if testImg == 304
    strVz = sprintf( '%s/HoguetBall/Cones', folder);
    imageName = 'Cones';
  end
  if testImg == 401
    strVz = sprintf( '%s/Box0_5Frames/', folder);
    imageName = 'Box0_5Frame';
  end
  if testImg == 402
    strVz = sprintf( '%s/Box1_5Frames/', folder);
    imageName = 'Box1_5Frame';
  end
  if testImg == 403
    strVz = sprintf( '%s/Box4_5Frames/', folder);
    imageName = 'Box4_5Frame';
  end
  if testImg == 404
    strVz = sprintf( '%s/Box5_5Frames/', folder);
    imageName = 'Box5_5Frame';
  end

  % WOW GL encodes the image coordiantes by -255/256 .. 255/256
  % when having 256 pixels in each dimension or -(N-1/N) .. N-1/N
  % and -(M-1/M) .. M-1/M with a stepsize of 2/N of course (N pixel).
  % This leads to weird recalculations everywhere!
  
  if testImg >= 55 && testImg <= 58
    % load images
    ImgL{1} = double(imread(sprintf('%s/%s_t0_0.png', strVz, imageName)))/255;
    ImgL{2} = double(imread(sprintf('%s/%s_t1_0.png', strVz, imageName)))/255;
    if testImg == 55
      flowGT = readCarpetBall( size(ImgL{1},1), size(ImgL{1},2), sprintf('%s%s', strVz, '/GroundTruth/' ));
    else
      [IN IM] = size(ImgL{1});
      flowGT = zeros(IN, IM, 4);
    end
    ImgL{1} = ImgL{1} (:,:,1);
    ImgL{2} = ImgL{2} (:,:,1);
    
    for j=1:p.nCams
      filename_t0 = sprintf('%s/%s_t0_%d.png', strVz, imageName, j);
      filename_t1 = sprintf('%s/%s_t1_%d.png', strVz, imageName, j);
      fid = fopen(filename_t0, 'r');
      if (fid < 0)
        p.nCams = j;
        break;
      end
      fclose(fid);
      ImgR{1} = double(imread(filename_t0))/255;
      ImgR{2} = double(imread(filename_t1))/255;
      ImgR{1} = ImgR{1} (:,:,1);
      ImgR{2} = ImgR{2} (:,:,1);
      data_supp(j).I = ImgR;
%      data_supp(j).I{1} = ImgR{1};
%      data_supp(j).I{2} = ImgR{2};
    end
    ImgR = data_supp(1).I; % later copied into
%    ImgR{1}  = ImgR{1};
%    ImgR{2}  = ImgR{2};
  elseif testImg < 200
    if testImg < 100
      flowGT = readFlowFile( sprintf('%s/data.sfl', strVz) );
      ImgL{1}  = double(imread(sprintf('%s/leftCamT0.png', strVz)))/255;
      ImgR{1}  = double(imread(sprintf('%s/rightCamT0.png', strVz)))/255;
      ImgL{2} = double(imread(sprintf('%s/leftCamT1.png', strVz)))/255;
      ImgR{2} = double(imread(sprintf('%s/rightCamT1.png', strVz)))/255;
    else % the names are flipped
      ImgR{1} = double(imread(sprintf('%s/leftCamT0.bmp', strVz)))/255;
      ImgL{1} = double(imread(sprintf('%s/rightCamT0.bmp', strVz)))/255;
      ImgR{2} = double(imread(sprintf('%s/leftCamT1.bmp', strVz)))/255;
      ImgL{2} = double(imread(sprintf('%s/rightCamT1.bmp', strVz)))/255;

      [IN IM] = size(ImgL{1});
      flowGT = zeros(IN, IM, 4);
    end

    ImgL{1} = ImgL{1}(:,:,1);
    ImgL{2} = ImgL{2}(:,:,1);
    ImgR{1} = ImgR{1}(:,:,1);
    ImgR{2} = ImgR{2}(:,:,1);
    
  elseif testImg > 400

    if p.frames == 2
      dataIds = [  3,1,4,0];
      imgIds  = [2,3,1,4,0];
    elseif p.frames == 1
      dataIds = [  3,1];
      imgIds  = [2,3,1];
    else
      dataIds = [  3];
      imgIds  = [2,3];
    end
    flowGT = readFlowFile( sprintf('%s/data%d.sfl', strVz, dataIds(1)) );
    for i = 2:numel(dataIds)
      sPos     = 1+min(i-1,1)+3*(i-1);
      ePos     = 1+3*(i);
      flowTemp = readFlowFile( sprintf('%s/data%d.sfl', strVz, dataIds(i)) );
      flowGT(:,:,sPos:ePos) = flowTemp(:,:,2:4);
    end
    
    for i=numel(dataIds):-1:3
      flowGT(:,:,i*3+1) = flowGT(:,:,i*3+1) - flowGT(:,:,(i-2)*3+1);
      flowGT(:,:,i*3  ) = flowGT(:,:,i*3  ) - flowGT(:,:,(i-2)*3  );
      flowGT(:,:,i*3-1) = flowGT(:,:,i*3-1) - flowGT(:,:,(i-2)*3-1);
    end
    
    for i=1:numel(imgIds)
      ImgL{i} = double(imread(sprintf('%s/leftCamT%d.png', strVz, imgIds(i))))/255;
      ImgR{i} = double(imread(sprintf('%s/rightCamT%d.png', strVz, imgIds(i))))/255;
      ImgL{i} = ImgL{i}(:,:,1);
      ImgR{i} = ImgR{i}(:,:,1);
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if testImg == 301
    
    [ImgL{1}, ImgL{2}, flowGT, R_l, data_supp, imageName, strVz] = ...
      loadHoguetFlow(strVz , p);
    [M N]  = size(ImgL{1});

%    [stereoD validMap] = buildDepthFromStereo(N, M, ceil (N/512), ImgL{1}, data_supp(1).It0, R_l, data_supp(1).R, data_supp(1).epi, ...
%      data_supp(1).Rot, data_supp(1).Tra, p, [0 0 0 0], [0 0 0 0], 40, 10000, 128, 1, 0, 0, testImg);

%    [stereoD validMap] = buildDepthFromStereo(N, M, 1, ImgL{1}, data_supp(1).I{1}, ...
%      R_l, data_supp(1).R, data_supp(1).epi, data_supp(1).Rot, data_supp(1).Tra,...
%      [0 0 0 0], [0,0,0,0], 40, 10000, 128, 1, -1, 0, 0);

    stereoD = flowGT(:,:,1);
    validMap = ones(M,N);
    
    ImgR{1} = data_supp(1).It0;
    ImgR{2} = data_supp(1).It1;
    R_r = data_supp(1).R;
    F = data_supp(1).F;
    data_supp(1).Ft = zeros(3);
    epi = data_supp(1).epi;
    
    popl = data_supp(1).popRef;
    popr = data_supp(1).pop;
    Rot  = data_supp(1).Rot;
    Tra  = data_supp(1).Tra;

%    return;
  end

    if testImg >= 302 && testImg <= 304

    [ImgL{1}, ImgL{2}, flowGT, R_l, data_supp, imageName, strVz] = ...
      loadStereoFlow(strVz , p, imageName);
    [M N]  = size(ImgL{1});

%    [stereoD validMap] = buildDepthFromStereo(N, M, ceil (N/512), ImgL{1}, data_supp(1).I{1}, R_l, data_supp(1).R, epi, ...
%      Rot, Tra, [0 0 0 0], [0 0 0 0], 40, 10000, 80, 9, 32, 1, testImg);
if testImg == 304
    [stereoD validMap] = buildDepthFromStereo(N, M, 1, ImgL{1}, data_supp(1).I{1}, ...
      R_l, data_supp(1).R, data_supp(1).epi, data_supp(1).Rot, data_supp(1).Tra,...
      [0 0 0 0], [0,0,0,0], 10, 100, 64, 5, -1, 0, testImg);
elseif testImg == 303
    [stereoD validMap] = buildDepthFromStereo(N, M, 1, ImgL{1}, data_supp(1).I{1}, ...
      R_l, data_supp(1).R, data_supp(1).epi, data_supp(1).Rot, data_supp(1).Tra,...
      [0 0 0 0], [0,0,0,0], 20, 100, 48, 5, -1, 0, testImg); 
else
    [stereoD validMap] = buildDepthFromStereo(N, M, 1, ImgL{1}, data_supp(1).I{1}, ...
      R_l, data_supp(1).R, data_supp(1).epi, data_supp(1).Rot, data_supp(1).Tra,...
      [0 0 0 0], [0,0,0,0], 20, 85, 32, 7, -1, 0, testImg);
end
%    stereoD = flowGT(:,:,1);
    validMap = ones(M,N);
    
    ImgR{1} = data_supp(1).I{1};
    ImgR{2} = data_supp(1).I{2};
    R_r = data_supp(1).R;
    F = data_supp(1).F;
    data_supp(1).Ft = zeros(3);
    epi = data_supp(1).epi;
    
    popl = data_supp(1).popRef;
    popr =data_supp(1).pop;
    Rot = data_supp(1).Rot;
    Tra = data_supp(1).Tra;

%    return;
  end

  
  if testImg < 55 || testImg > 400 % MY ARTIFICIAL
    [Pl, Pr, F, epi, R_l, R_r, Rot, Tra] = readMatrixFile(sprintf('%s/matrix.sfl',strVz));
    
    Rr   = Pr(:,1:3);
    popl = Pl(:,end);
    popr = Pr(:,end);

    data_supp(1).R = Rr;
    data_supp(1).F = F;
    data_supp(1).Ft = zeros(3);
    data_supp(1).epi = epi;
    data_supp(1).popRef = popl;
    data_supp(1).pop = popr;
    data_supp(1).Rot = Rot;
    data_supp(1).Tra = Tra;
    
    [data_supp(1).Kl, data_supp(1).Rl, data_supp(1).Tl, pp] = cameraParameters ( Pl );
    [data_supp(1).Kr, data_supp(1).Rr, data_supp(1).Tr, pp] = cameraParameters ( Pr );
    data_supp(1).Tl = - inv(data_supp(1).Rl) * data_supp(1).Tl;
    data_supp(1).Tr = - inv(data_supp(1).Rr) * data_supp(1).Tr;
    data_supp(1).Rl = inv(data_supp(1).Rl);
    data_supp(1).Rr = inv(data_supp(1).Rr);
    data_supp(1).I = ImgR;
    [N M] = size(ImgR{1});
    Viewport = [ (N)/2 0 (N+1)/2 ; 0 (M)/2 (M+1)/2 ; 0  0  1 ];
    T = [ 2/(N) 0 -(N+1)/(N) ; 0 2/(M) -(M+1)/(M) ; 0  0  1];
    data_supp(1).F = T' * data_supp(1).F * T;
    data_supp(1).R = Viewport * data_supp(1).R;
    data_supp(1).epi = Viewport * data_supp(1).epi;
    data_supp(1).pop = Viewport * data_supp(1).pop;
    data_supp(1).Kl = Viewport * data_supp(1).Kl;
    data_supp(1).Kr = Viewport * data_supp(1).Kr;

    
    [cam, ref] = generateStructures ( data_supp, ImgL, data_supp.R );
    % nothing better here :
    cam.Inew = {ref.I(2).I, cam.I(2).I};
    cam.Iold = {ref.I(1).I, cam.I(1).I};
    % well : how ?    
    [flow2DGt(:,:,1), flow2DGt(:,:,2), flow2DGt(:,:,3), flow2DGt(:,:,4) ] ...
      = convert3Dto2D(ref, cam, ... 
      imresize(-flowGT(:,:,1), size(ImgR{1}), 'bilinear'), ...
      imresize(-flowGT(:,:,2), size(ImgR{1}), 'bilinear'), ...
      imresize(-flowGT(:,:,3), size(ImgR{1}), 'bilinear'), ...
      imresize(-flowGT(:,:,4), size(ImgR{1}), 'bilinear') );
    % well either or depth as - depth ?
%    flow2DGt(:,:,1) = flow2DGt(:,:,1);
    flow2DGt_noc = flow2DGt;

    u  = ones(M,N,3);u(:,:,1) = repmat( [1:N],  M, 1 );u(:,:,2) = repmat( [1:M]', 1, N );
    dispOff = -flow2DGt(:,:,1) + u(:,:,1) < 1 | -flow2DGt(:,:,1) + u(:,:,1) > N;
    flowOff = flow2DGt(:,:,2) + u(:,:,1) < 1 | flow2DGt(:,:,2) + u(:,:,1) > N | ...
      flow2DGt(:,:,3) + u(:,:,2) < 1 | flow2DGt(:,:,3) + u(:,:,2) > M;
    flow2DGt_noc(dispOff) = -1;
    flow2DGt_off = ones(M,N);
    flow2DGt_off(flowOff) = 0;
    flow2DGt_noc(:,:,4) = flow2DGt_off;
    return;
    
    Rot = data_supp(1).Rr;
    Tra = data_supp(1).Tr;
    popr = Pr(:,end);
    popl = Pl(:,end);
    stereoD = flowGT(:,:,1);
    validMap = ones (size(stereoD));
  elseif testImg >= 200 && testImg < 300 % PTAM

    % read images
    nr = testImg - 200;
    imgName = sprintf('im%04d', nr );
    
    p.eqOn.lrt0 = 1; % left right at timestep 0 (depth only)
    p.eqOn.lrt1 = 0;p.eqOn.ll = 0;p.eqOn.rr = 0;
    p.eqOn.lt0rt1 = 0; p.eqOn.lt1rt0 = 0;

    if 1
      [R_l, K, data_supp, ImgL, Rot, Tra] = ...
        readMatrixFileTrackingUnDistort(strVz, nr, data_supp, imgName, p);
    else
      
      ImgL{1}  = double(imread(sprintf('%s/%s_%02d.pgm', strVz, imgName, 0)))/255;
      %Ilt1 = double(imread(sprintf('%s/%s_%02d.pgm', strVz, imgName ,0)))/255;
      ImgL{2} = ImgL{1};
      % read matrix
      [M N] = size (ImgL{1});
      [R_l, K, data_supp] = readMatrixFileTracking(strVz, nr, N, M, data_supp, p.nCams);
      w = 0.911538;
      [ImgL{1} minX minY maxX maxY] = computeUndistortedImagesPTAM ( ImgL{1},  K, w );
      
      for j = 1:min(p.nCams, 4)
        ImgR{1} = double(imread(sprintf('%s/%s_%02d.pgm', strVz, imgName, j)))/255;
        data_supp(j).I{1} = computeUndistortedImagesPTAM ( ImgR{1}, K, w );
        data_supp(j).I{2} = data_supp(j).I{1};
      end

    %{
    % cut out out of buonds region: cut out AND new K
    Viewport = [ (N)/2 0 (N+1)/2 ; 0 (M)/2 (M+1)/2 ; 0  0  1 ];
    for j = 1:min(p.nCams, 4)
      data_supp(j).I{1} = data_supp(j).I{1}( minY+1:maxY-1, minX+1:maxX-1 );
      data_supp(j).I{2} = data_supp(j).I{2}( minY+1:maxY-1, minX+1:maxX-1 );
      
    end
    ImgL{1} = ImgL{1} ( minY+1:maxY-1, minX+1:maxX-1 );
    ImgL{2} = ImgL{2}( minY+1:maxY-1, minX+1:maxX-1 );
%}
    end

    % no distortion any more
    D1 = [0 0 0 0]; D2 = [0 0 0 0];
    [M N]         = size (ImgL{1});

    flowGT        = zeros(M,N,4);
    flowGT(:,:,1) = 1 * ones(M,N); % in meter as initialization

%    [estD] = buildDepthFromStereo(N, M, 1, ImgL{1}, data_supp(4).I{1}, R_l, ...
%      data_supp(4).R, data_supp(4).Rot, data_supp(4).Tra, ...
%      [0 0 0 0], [0 0 0 0], 0.55, 1.4, 80, 5);

%    flowGT(:,:,1) = estD;

    flowGT = readFlowFile( sprintf('%s/Erdbeer_011_2DFlow.sfl', strVz) );

    % for debugging (see below)
    popr = data_supp(1).pop;
    popl = data_supp(1).popRef;
    R_r  = data_supp(1).R;
    epi  = data_supp(1).epi;
    F    = data_supp(1).F;
    ImgR{1} = data_supp(1).I{1};
    ImgR{2} = data_supp(1).I{1};
    
    %% Seems to be working:
    %{
  %same as above
  for j=1:size(data_supp,2)
    T = [ (N)/2 0 (N+1)/2 ; 0 (M)/2 (M+1)/2 ; 0  0  1];
    %T = [ 2/(N) 0 -(N+1)/(N) ; 0 2/(M) -(M+1)/(M) ; 0  0  1];
    data_supp(j).F = inv(T)' * data_supp(j).F * inv(T);
    data_supp(j).F = data_supp(j).F / norm(data_supp(j).F,'fro');
  end

  % invertZ -1
%  drawEpiLines( ImgL{2}, data_supp(1).I{2}, data_supp(1).F', 7, 0, 0);
%  drawEpiLines( data_supp(1).I{2}, ImgL{2}, data_supp(1).F, 6, 0, 0);
  
  % invertZ : 1 NO
  imgId = 4;
  drawEpiLines( ImgL{1}, data_supp(imgId).I{1}, data_supp(imgId).F, 7, 0, 0);
  drawEpiLines( data_supp(imgId).I{1}, ImgL{1}, data_supp(imgId).F', 6, 0, 0);
    %}
    
    % probe inpainting:
    %{
  [M, N, D] = size(flowGT);
  Viewport = [ (N)/2 0 (N+1)/2 ; 0 (M)/2 (M+1)/2 ; 0  0  1];

  % keep out when transforming normalized coordinates
  R_l_ = Viewport * R_l;
  [u, lambda_uv] = Compute_lambdaUV(R_l_, N, M);
%  flowGT(:,:,1) = 1 * ones(M,N); % in meter ?
  for j=1:p.nCams
    R_r_ = Viewport * data_supp(j).R;
    epi_ = Viewport * data_supp(j).epi;

    ur = Compute_ur_all_d(R_l_, R_r_, epi_, u, lambda_uv, flowGT(:,:,1) );
    Irt_warped  = interp2(data_supp(j).I{1} , ur (:,:,1), ur (:,:,2), 'linear');
    figure(1); imshow(ImgL{1}), figure(1+j), imshow(Irt_warped);
  end
    %}
    
  elseif testImg >= 100 && testImg < 300
    
    [M N] = size (ImgL{1});
    flowGT = zeros(M,N,4);
    
    %{
      sxl = 2614;
      syl = 2616;
      x0l = 649.60;
      y0l = 536.20;
      
      sxr = 2615;
      syr = 2615;
      x0r = 629.69;
      y0r = 518.65;
      
      R = rodrigues( [-0.00529,   -0.18582,  -0.00983]);
      T = [167.34330   -1.05333  19.37890];
      Rt = [R , T'];
      %      Rt = [inv(R), -inv(R)*T'];
      
      Viewport = [ (N)/2 0 (N+1)/2 ; 0 (M)/2 (M+1)/2 ; 0  0  1 ];
      %FlipY = [ 1 0 0 ; 0 -1 M+1; 0 0 1];
      FlipY = [ 1 0 0 ; 0 1 0; 0 0 1];
      
      K_l   = [sxl 0 x0l; 0 syl y0l; 0 0 1];
      K_r   = [sxr 0 x0r; 0 syr y0r; 0 0 1];
      R_l = [1 0 0 0; 0 1 0 0; 0 0 1 0]; % original line
      R_r = Rt;

      D1 = [0.02501 -0.59719 -0.00305 -0.00316];
      D2 = [-0.08304   0.40370   -0.00063   0.00073];
      
      [ImgL{1} ImgR{1} K1 K2] = Undistort_MexWindows( M, N, uint8(ImgL{1}*255), uint8(ImgR{1}*255), K_l, K_r, D1, D2);
      ImgL{1} = ImgL{1}/255;
      ImgR{1} = ImgR{1}/255;
      [ImgL{2} ImgR{2} K1 K2] = Undistort_MexWindows( M, N, uint8(ImgL{2}*255), uint8(ImgR{2}*255), K_l, K_r, D1, D2);
      ImgL{2} = (ImgL{2})/255;
      ImgR{2} = (ImgR{2})/255;
      % NOW THERE IS NO DISTORTION ANY MORE
      D1 = zeros (size(D1));
      D2 = zeros (size(D1));
    %}
    
    [ R_l R_r popl popr epi F Dl Dr ImgL{1} ImgL{2} ImgR{1} ImgR{2} Rot Tra data_supp] = ...
      computeMatrixFileRealWorldSession1 ( M, N, [649.60 536.20], [629.69 518.65], ...
      [2614 2616], [2615 2615], [-0.00529, -0.18582, -0.00983], ...
      [167.34330   -1.05333  19.37890], ImgL{1}, ImgL{2}, ImgR{1}, ImgR{2}, data_supp, ...
      [0.02501 -0.59719 -0.00305 -0.00316], [-0.08304 0.40370 -0.00063 0.00073]);
    
    if ~exist( sprintf('%s/leftCamT0_unwarped.png', strVz), 'file')
      imwrite(ImgL{1}, sprintf('%s/leftCamT0_unwarped.png', strVz));
      imwrite(ImgR{1}, sprintf('%s/rightCamT0_unwarped.png', strVz));
      imwrite(ImgL{2}, sprintf('%s/leftCamT1_unwarped.png', strVz));
      imwrite(ImgR{2}, sprintf('%s/rightCamT1_unwarped.png', strVz));
    end
    
    if 1 || ~exist(sprintf('%s/RealWorld101.sfl', strVz), 'file')
%      [stereoD validMap] = buildDepthFromStereo(N, M, ceil (N/512), ImgL{1}, ImgR{1}, R_l, R_r, epi, ...
%        Rot, Tra, [0 0 0 0], [0 0 0 0], 750, 1800, 128, 9, -16);%-16);

      [stereoD validMap] = buildDepthFromStereo(N, M, ceil (N/512), ImgL{1}, ImgR{1}, R_l, R_r, epi, ...
        Rot, Tra, [0 0 0 0], [0 0 0 0], 750, 1800, 80, 9, 32, 1, testImg);

      
      if size(stereoD,1) ~= M
        [ flowGT(:,:,1) ] = interp2OccFull(stereoD, M, N, validMap);
      else
        flowGT(:,:,1) = stereoD;
      end
      Viewport = [ (N)/2 0 (N+1)/(2) ; 0 (M)/(2) (M+1)/(2) ; 0  0  1 ];

      % writing images to check the quality of the stereo crap
      %{
      Rr__  = Viewport * R_r;
      Rl__  = Viewport * R_l;
      epi__ = Viewport * epi;
      [u, lambda_uv] = Compute_lambdaUV(Rl__, N, M);
      [ur, ~, ulx, ~, ~, ~, ~, urx] = computeImageProperties( Rl__, Rr__, ...
        u, epi__, lambda_uv, flowGT(:,:,1), flowGT(:,:,2), flowGT(:,:,3), flowGT(:,:,4), ...
        0, 0, 0, 1);

      % already done, just to check that it sucks
%      offPixel = ur(:,:,1) < 0 | ur(:,:,1) > N |ur(:,:,2) < 0 | ur(:,:,2) > M;
%      d = TvL1Matrix_Mex(1000, M, N, int32(offPixel),  flowGT(:,:,1),  flowGT(:,:,1));
      [Irt1_warped, ~, ~, Irt_warped, ~, ~, Ilt1_warped ] = ...
        computeWarpedImages( ImgL{2}, ImgR{1}, ImgR{2}, 0, 0, 0, 0, 0, 0, ulx, ur, urx, 0);
      figure(1), imshow(ImgL{1}), figure(3), imshow(Irt_warped)%, figure(2), imshow(ImgR{1}),
%}
      writeFlow(Viewport*R_l, flowGT, sprintf('%s/RealWorld101.sfl', strVz) , 0);
    else % computed before
      flowGT = readFlowFile( sprintf('%s/RealWorld101.sfl', strVz) );
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  elseif testImg >= 55 && testImg < 59% imgNr in [55..58]
    [M N] = size (ImgL{1});
    
    %    [Pl, Pr, F, epi, R_l, R_r] = readMatrixFileCarpet( strVz, M, N);
    [R_l, data_supp, Rot, Tra] = readMatrixFileCarpetFull( data_supp, sprintf( '%s%s', strVz, '/calibration/'), N, M, size(data_supp,2));
    
    popr = data_supp(1).pop;
    popl = data_supp(1).popRef;
    
    % for debugging (see below)
    R_r = data_supp(1).R;
    epi = data_supp(1).epi;
    F = data_supp(1).F;

    %{
    %    print GT imagesc carpet Ball
    for j = 1:size(data_supp,2)
      [M N] = size (ImgL{1});
      Viewport = [ (N)/2 0 (N+1)/2 ; 0 (M)/2 (M+1)/2 ; 0  0  1 ];
      Rr__  = Viewport * data_supp(j).R;
      Rl__  = Viewport * R_l;
      epi__ = Viewport * data_supp(j).epi;
      [u, lambda_uv] = Compute_lambdaUV(Rl__, N, M);
      write4WarpedImages( Rl__, Rr__, epi__, u, lambda_uv, ...
        flowGT(:,:,1), flowGT(:,:,2), flowGT(:,:,3), flowGT(:,:,4), flowGT(:,:,1), flowGT(:,:,2), flowGT(:,:,3), flowGT(:,:,4),...
        ImgL{1}, ImgL{2}, data_supp(j).I{1}, data_supp(j).I{2}, './carpetWarped', sprintf('camera_%d', j), 1 );
    end
%}    
% debug
%  
%{
    for j=1:size(data_supp,2)
      add = size(data_supp,2);
      T = [ (N)/2 0 (N+1)/2 ; 0 (M)/2 (M+1)/2 ; 0  0  1];
      %T = [ 2/(N) 0 -(N+1)/(N) ; 0 2/(M) -(M+1)/(M) ; 0  0  1];
      data_supp(j).F = inv(T)' * data_supp(j).F * inv(T);
      data_supp(j).F = data_supp(j).F / norm(data_supp(j).F,'fro');
      drawEpiLines( ImgL{1}, data_supp(j).I{1}, data_supp(j).F, 7+add+j, 0, 0);
      drawEpiLines( data_supp(j).I{1}, ImgL{1}, data_supp(j).F',6+j, 0, 0);
    end
%}

    if testImg ~=55
      if ~exist(sprintf('%s/StereoMaria.sfl', strVz), 'file')
%        [estD] = buildDepthFromStereo(N, M, ceil (N/512), ImgL{1}, data_supp(1).I{1}, ...
%          R_l, data_supp(1).R, data_supp(1).epi, data_supp(1).Rot, data_supp(1).Tra,...
%          [0 0 0 0], [0,0,0,0], 950, 2000, 64, 7, 0);
        
        [stereoD validMap] = buildDepthFromStereo(N, M, 1, ImgL{1}, data_supp(1).I{1}, ...
          R_l, data_supp(1).R, data_supp(1).epi, data_supp(1).Rot, data_supp(1).Tra,...
          [0 0 0 0], [0,0,0,0], 980, 2200, 96, 11, 10, 1, testImg);    

        if testImg == 58 % maria
          [stereoD2 validMap2] = buildDepthFromStereo(N, M, 1, ImgL{1}, data_supp(2).I{1}, ...
            R_l, data_supp(2).R, data_supp(2).epi, data_supp(2).Rot, data_supp(2).Tra,...
            [0 0 0 0], [0,0,0,0], 980, 2200, 80, 11, -150, 1, -testImg);

          stereoD3 = stereoD .* validMap + (1-validMap) .* stereoD2;
          validMap3 = max( validMap , validMap2);
          stereoD = TvL1Matrix_Mex(1000, M, N, int32(1-validMap3), stereoD3, stereoD3);
        
        end        
        if size(stereoD,1) ~= M
          [ flowGT(:,:,1) ] = interp2OccFull(stereoD, M, N, validMap);
        else
          flowGT(:,:,1) = stereoD;
        end
        writeFlow(R_l, flowGT, sprintf('%s/StereoMaria.sfl', strVz) , 0);
%        flowGT(:,:,1) = imresize(estD, [M N], 'bilinear');
      else % computed before
        % no vlidity map so far - maybe store as well
        flowGT = readFlowFile( sprintf('%s/StereoMaria.sfl', strVz) );
      end
    end
%    flowGT = readFlowFile( sprintf('%s/2DFlowGT.sfl', strVz) );
  end

  data_supp(1).I = ImgR;
  data_supp(1).R = R_r;
  data_supp(1).F = F;
  data_supp(1).Ft = zeros(3);
  data_supp(1).epi = epi;

  data_supp(1).popRef = popl;
  data_supp(1).pop    = popr;
  data_supp(1).Rot    = Rot;
  data_supp(1).Tra    = Tra;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % rescale the flowGT to stop at a different resolution than the given one
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [M N D]  = size (flowGT);
  [IM IN]  = size(ImgL{1});
  IltReal  = ImgL{1};
  Ilt1Real = ImgL{2};
  
  if p.stopAtResolution > 0 && p.stopAtResolution ~= 1
    if p.reStartLevel >= 0 && p.reStartLevel ~= 99
      %    N * par.pyramid_factor^k < stopAtResolution; % determine k
      %    log(N) + log (par.pyramid_factor) * k < log (stopAtresolution);
      % -> k < (log(stopAtResolution) - log(N)) /  log(par.pyramid_factor);
      k1 = floor( (log(p.stopAtResolution) - log(N)) / log(p.pyramid_factor) );
      k2 = floor( (log(p.stopAtResolution) - log(M)) / log(p.pyramid_factor) );
      k = max(k1,k2);%k = max(0, max(k1,k2));
      N = ceil(N * p.pyramid_factor^k);
      M = ceil(M * p.pyramid_factor^k);
    else
      M = ceil(M/N * p.stopAtResolution);
      N = p.stopAtResolution;
    end
    %    N = IN; M = IM;
    % hack for higher resolution images from graound truth for nicer display properties
    if size(validMap,1) > 1 && testImg ~= 301
      [ flowGT_(:,:,1) ] = interp2OccFull(stereoD, M, N, validMap);
    else
      flowGT_(:,:,1) = imresize(flowGT(:,:,1), [M N], 'bilinear');
    end
    for i=2:size(flowGT,3)
      flowGT_(:,:,i) = imresize(flowGT(:,:,i), [M N], 'bilinear');
    end
    flowGT = flowGT_;
    clear flowGT_;
  end
  if (IN ~= N || IM ~= M)
    for i=1:numel(ImgL)
      ImgL{i}  = imresize(ImgL{i}, [M N], 'bilinear');
    end
%    ImgL{2} = imresize(ImgL{2}, [M N], 'bilinear');
    for j = 1:size(data_supp,2)
      for i=1:numel(data_supp(j).I)
        data_supp(j).I{i} = imresize(data_supp(j).I{i}, [M N], 'bilinear');
      end
%      data_supp(j).I{2} = imresize(data_supp(j).I{2}, [M N], 'bilinear');
    end
  end
  
  % Swap images according to desired reference view
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if p.refView >= 0
    
    p.refView = mod (p.refView, 3*size(data_supp,2) );
    cam = floor (p.refView / 3*size(data_supp,2)) + 1;
    pos = mod ( p.refView, 3 );
    
    % modify GT:
    [M, N, D] = size(flowGT);
    Viewport = [ (N)/2 0 (N+1)/2 ; 0 (M)/2 (M+1)/2 ; 0  0  1];
    
    % keep out when transforming normalized coordinates
    R_l_ = Viewport * R_l;
    [u, lambda_uv] = Compute_lambdaUV(R_l_, N, M);
    R_r_ = Viewport * data_supp(cam).R;
    epi_ = Viewport * data_supp(cam).epi;
    
    global occlusionDepth;
    occlusionDepth = 0.01 * occlusionDepth;
    
    if pos == 0
      [uGT, ~, ~, oMapGT, uDepthGT] = Compute_ur_all_d_weight(R_l_, R_r_, epi_, ...
        u, lambda_uv, flowGT(:,:,1), 1, 0, 1 );
    elseif pos == 1
      [uGT, ~, ~, ~, ~, ~, ~, ~, ~, oMapGT, uDepthGT] = ...
        Compute_ulx_all_d_weight(R_l_, u, lambda_uv, flowGT(:,:,1), ...
        flowGT(:,:,2), flowGT(:,:,3), flowGT(:,:,4), 1, 0, 1 );
    else
      [uGT, ~, ~, ~, ~, ~, ~, ~, ~, oMapGT, uDepthGT] = ...
        Compute_urx_all_d_weight(R_l_, R_r_, u, epi_, lambda_uv, ...
        flowGT(:,:,1), flowGT(:,:,2), flowGT(:,:,3), flowGT(:,:,4), 1, 0, 1 );
    end
    occlusionDepth = 100 * occlusionDepth;
    
    addBorder = 0.5;
    m = (uGT(:,:,1) > N+addBorder) | (uGT(:,:,1) < 1-addBorder) | ...
      (uGT(:,:,2) > M+addBorder) | (uGT(:,:,2) < 1-addBorder);
    oMapGT (m) = 0;
    
    % ignore cam for the first time
    if pos == 0 || pos == 2 % switch right, left
      [Kl Rl] = rq(R_l);
      [Kr Rr] = rq(data_supp(1).R);
      
      T = ones (3,1);
      T(:) = data_supp(1).Tra(:);
      
      Rt    = [inv(Rr) inv(Rr)*(-T)] * [Rr T; 0 0 0 1];
      RtNew = [inv(Rr) inv(Rr)*(-T)] * [Rl [0 0 0]'; 0 0 0 1];
      
      RotNew = RtNew(:, 1:3);
      TraNew = RtNew(:, end);
      
      % the flow must be rotated, since we made RotNew to be the identity
      % But it might be that the Rotation in the start is not (possible some -1 on diagonal)
      RFlow = RotNew * inv(Rl);
      
      Pref = Kr * Rt;
      Pnew = Kl * RtNew;
      
      [R_l, data_supp(1).popl, data_supp(1).epi, data_supp(1).F, ...
        data_supp(1).R, data_supp(1).popr] = getEpiF(Pref, Pnew);
      
      data_supp(1).Rot = RotNew;
      data_supp(1).Tra = TraNew;
      
      % switch images
      tempIlt0 = ImgL{1};
      tempIlt1 = ImgL{2};
      ImgL{1}  = data_supp(1).I{1};
      ImgL{2} = data_supp(1).I{2};
      data_supp(1).I{1} = tempIlt0;
      data_supp(1).I{2} = tempIlt1;
      
    else
      RFlow  = eye(3);
    end
    
    % swap images
    if pos >= 1
      
      tempIlt0 = ImgL{1};
      tempIrt0 = data_supp(1).I{1};
      
      ImgL{1}  = ImgL{2};
      data_supp(1).I{1} = data_supp(1).I{2};
      
      ImgL{2} = tempIlt0;
      data_supp(1).I{2} = tempIrt0;
      
    end
    
    flowGTNew = zeros( size( flowGT ) );
    missedMap = zeros( M,N );
    for i=1:N
      for j=1:M
        if oMapGT(j,i)
          jj = round(uGT(j,i,1));
          ii = round(uGT(j,i,2));
          
          flow_ij = flowGT(j,i,2:4);flow_ij = squeeze(flow_ij);
          flow_ij = RFlow * flow_ij;
          
          flowGTNew(ii, jj, 1)   = uDepthGT(j,i);
          flowGTNew(ii, jj, 2:4) = flow_ij;
          
          missedMap(ii,jj) = 1;
          
        end
      end
    end
    
    % fill holes
    for i = 1:4
      dispV = ones(N*M, 1);
      temp = squeeze(flowGTNew(:, :, i));
      dispV(:) = temp(:);
      [dispV, fX_, its_] = minimize(dispV, @tv_l1_minimizeFunctionD_noL1, ...
        [1000 150000000*N], temp, 1/1000000000, missedMap, ones(M,N), -1);
      temp(:) = dispV(:);
      flowGTNew(:, :, i) = temp;
    end
    
    flowGT =  flowGTNew;
    if pos > 0 % timestep reversed
      flowGT(:,:,2:4) = -flowGTNew(:,:,2:4);
    end
  end
