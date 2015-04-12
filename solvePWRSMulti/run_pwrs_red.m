%%% Demonstration function, reads kitti images and delivers error measures
%%% and if storeOutput == 1 stores the results in the folder storeFolder
%%% To run on other input data please change/modifiy the function
%%% loadKittiFlow. As a sceneflow method the code needs the calibration of
%%% the cameras besides the image data, thus such a function loading the 
%%% data was prefered.
%%% 
%%% Please read through the comments directly below to receive an
%%% explanation of the paramters. The image data is encoded by scene-number: imgNr
%%% and by frame number: subImg.
%%% Please note that consecutive images in time are assumed to have
%%% consecutive subImg numbers. The three frame method assumes computes
%%% a solution for the middle frame and assumes the past and next frames to
%%% exist. Further to reuse proposals from the previous solution the
%%% function must have been run on the frame number before. 
%%% Three frames can be turned off by setting p.use3Frames = false; below.
%%%
%%% returns the 2d scene flow: disparity and flow and disparity difference
%%% also saves the disparities and the flow if storeOutput==1, 
%%% these can be read with flow_read( 'filename' )
%%%
%%% dt: data-weight, ds: smooth-weight (\lambda), ts: relative motion smoothness
%%% dj,tj: cutoff smoothness in pixel (eta); ps: per pixel patchsize: N_s;
%%% lr: local replacement: on/off; re: refinement on/off
%%% segWeight (mu); pseg,pjit,pego iccv13 precomputation on/off
%%% as: theta_oob/occ; tp: theta_mvp; ots: tolerated distance (epsilon)
function flow2d = run_pwrs_red(imgNr, storeFolder, subImg, ...
 dt, ds, ts, dj, tj, ps, lr, re, testing, segWeight, as, tp, ots )

%example: pic 138
%run_pwrs_red( [138], './2Frames/', 10, 0.4, 0.045, 1.0, 20, 20, 25 );
% example: set: p.saveProposals =true; p.use3Frames = true;p.usePrevProps = true;
%run_pwrs_red( [151,27], './3Frames/', [8,9,10], 0.4, 0.045, 1.0, 20, 20, 25 );

global doKittiErrors;doKittiErrors =0; % turn on/off permanent error evaluation
storeOutput = 1; % store the output as flow/stereo files in the folder: storeFolder 
% saved files can  be read with flow_read( 'filename' )

% folder with image data to read from -- ausmes certain structure
%dataFolder =  '../../../Desktop/work/data/';
%dataFolder = '../../data/data_stereo_flow/';
dataFolder =  '../kittidata/'; % local folder -- just a few images added for testing

p.subImg  = 10; % the frame number to proces, eg 10 -> pics 9 (if 3-frame version)
% 10 and 11 are loaded. To suppress the camera motion, eg. from a stereo
% rig mounted on top of the car however a solution of the previous frame
% must be available: .usePrevProps = true and p.use3Frames   = true
% 
% save proposals for future frames in p.tempFolder
p.saveProposals = true;
% 'temporal' folder for multi frame version -- save&read proposals of previous frames
%p.tempFolder   = '/cluster/scratch_xp/public/vogechri/pastproposalsTrainNewTime/';
p.tempFolder   = './scratch/pastproposals';

% stores results in this folder
p.storeFolder  = './scratch/test';
p.sFolder      = p.storeFolder;

% store proposals from previous frames here to reuse in the next frame
%pFolder = '/cluster/scratch_xp/public/vogechri/propsTrainNewTime/';
pFolder = './props/'; % store generated proposals here and reuse 
% formerly used to save paramters and misc. info about run

p.use3Frames   = true; % well use 3 not 2 frames -- assume pics loaded in ref/cam structure
p.usePrevProps = true; % use proposals from the previous time step -> loaded from tempFolder
% the procedure falls back to the standard 2 frame procedure in case there
% is no video data available
%
p.computeRflow = true;          % use flow from right camera 
p.generateMoreProposals = true; % demonstrates how to append additional proposals
%
p.fitSegs  = 1000;  % reduce # of fitted proposals to this value, default 1000 
p.ps  = 25; % patchsize of per pixel optimization: 25: 50x50 pixel, default 25
p.tj  = 20; % truncation value of motion smoothness, default 20
p.dj  = 20; % truncation of disparity smoothness, default 20
p.ts  = 1;  % smoothness of the motion field (p.ds*p.ts), default 0.1
p.ds  = 0.045;% smoothness of disparity field default 0.045 or 0.05
p.dt  = 0.4; % data penalty, default 0.4
p.segWeight    = 0.1; % default 0.1: segmentation weight (mu in journal)

%%%% view-consistency parameters %%%%
p.autoScale  = 0.75; % occlusion/out-of-bounds penalty (+0.1) default 0.7 or 0.75 -- so 0.7+0.1 == 0.8 == maxdata*0.5
p.vcPottsSeg = 0.15; % theta_mvp on segment level, default 0.15 or 0.1
p.vcPottsPix = 0.25; % theta_mvp on pixel level, default 0.25 
p.vcEpsSeg   = 0.15; % epsilon in vc data term on segment level, default 0.15 or 0.1
p.vcEpsPix   = 0.015;% epsilon in vc data term on pixel level, default 0.15 or 0.1

% for speedup disable, these 2 extension to the basic approach however deliver better results 
p.locRep=1;% optional : naive local replacement for view-consistent implementation;
p.refine=1;% optional : hierarchical refinement

%%%% size of expansion area in segments - should/could be wrt size of images
p.gx=8; % 8 kitti - can be adjusted to image size / relative size
p.gy=5; % 5 kitti
p.gridSize= 16; % kitti default 16 but depends on image size trades accuracy with speed
p.testing = 0; % whether to load kitti images from the training or test set

if exist('ots','var')
  p.vcEpsSeg = ots; % per segment tolerance epsilon
end
if exist('tp','var')
  p.vcPottsSeg = tp;
end
if exist('as','var')
  p.autoScale = as;
end
if exist('segWeight','var')
  p.segWeight = segWeight;
end
if exist('testing','var')
  p.testing = testing;
end
if exist('subImg','var')
  p.subImg = subImg;
end
if exist('storeFolder', 'var')
  p.storeFolder = storeFolder;
  p.sFolder = p.storeFolder;
end

if exist('ps', 'var')
  p.ps= ps;
end
if exist('tj', 'var')
  p.tj=tj;
end
if exist('dj', 'var')
  p.dj = dj;
end
if exist('ts', 'var')
  p.ts =ts;
end
if exist('ds', 'var')
  p.ds = ds;
end
if exist('dt', 'var')
  p.dt = dt;
end
if exist('lr', 'var')
p.locRep = lr;
end
if exist('re', 'var')
p.refine = re;
end
p.frames = 0;

global flow2DGt;
global flow2DGt_noc;
global Linux_;
Linux_ = 1;% different folders to load data from ?!

if ~isdeployed
  path(path,'./io/');  
  path(path,'./io/other/');
  path(path,'./mex/');
  path(path, './egomotion/');
  path(path,'./visualization/');
  path(path, './Segmentation');
  path(path, './weighting');
  path(path, './proposals');
  path(path,'./stuff/');
  path(path,'./pwrsf/');
  path(path,'../export_fig');%plotting - read the Readme
  path(path, '../sc/');%plotting - read the Readme
  path(path,'./stereo/');
  path(path,'../createFlow/');
  path(path,'../KittiIO/');
  path(path,'./ViewMappings');
end

testImages = imgNr;
subImages  = subImg;

permStr = '';
permStr = sprintf('%sdt:%1.2f ds:%.3f ts:%.3f dj:%.0f tj:%.0f ps:%d \n\n', ...
  permStr, p.dt, p.ds, p.ts, p.dj, p.tj, p.ps );
fprintf(2,permStr);

if ~exist(p.sFolder,'dir')
  mkdir(p.sFolder);
end

for testImg_ = 1:numel(testImages)  
  p.imgNr = testImages(testImg_);
  
  for subImg_ = 1:numel(subImages)
    p.subImg = subImages(subImg_ );
 
% can sometimes be unwanted, prevents from overwriting results:
%    if exist( sprintf('%s/RESULTS_K%03d_%02d_%s.txt', p.sFolder, p.imgNr, p.subImg, date), 'file' )
%      continue;
%    end
    
    doKitti = 1;
    if doKitti
      [cam, ref, imageName, flow2DGt, flow2DGt_noc] = loadKittiFlow(dataFolder , p.imgNr, p);
    else
      %%%% room to load your data: flow2DGt and flow2DGt_noc can be set to
      %%%% zeros if no GT is available; or set doKittiErrors to 0;
    end
    p.imageName = imageName;
    
    % construct the initial segmentation: cubes almost as good as super-pixel
    Seg = SegmentImageCube( ref.I(1).I, p.gridSize, p );%, 13); % 2nd parameter defines desired patchsize: the higher the larger the initial patches
    Seg = correctSegCenters(Seg);
    Seg = setWeights_patchSmooth( Seg, cam(1).Kl );

    tic
    N_prop=0; RT_prop=0;
    % function provides proposals by variational flow/stereo, SGM, other methods can be used (also additionally)

    if ~exist( sprintf( '%s/PropSolution%03d_%02d.mat', pFolder, p.imgNr, p.subImg), 'file');
      [N_prop, RT_prop]  = generateProposals(p, cam, ref, Seg );
      if ~exist(pFolder,'dir')
        mkdir(pFolder);
      end
      save( sprintf( '%s/PropSolution%03d_%02d.mat', pFolder, p.imgNr, p.subImg), 'Seg', 'N_prop', 'RT_prop');
    else
      load( sprintf( '%s/PropSolution%03d_%02d.mat', pFolder, p.imgNr, p.subImg));
    end      

    if size(N_prop,1) < 4
      N_prop = cat(1, N_prop, ones(1,size(N_prop, 2)));
    end    
    toc
    
    % new simplified version - does it work ?
    [flow2d, Energy] = pwrsf_v4 ( ref, cam, p, Seg, N_prop, RT_prop );
    
    doKittiErrors =1;
    [occErr, noccErr, epes] = getKittiErrSF ( flow2d(:,:,1), flow2d(:,:,2), flow2d(:,:,3) ); %, p, 1 );
    
    kittiStr = sprintf('DispPix-occ 2/3/4/5 %.3f & %.3f & %.3f & %.3f\nFlowPix-occ 2/3/4/5 %.3f & %.3f & %.3f & %.3f\n', occErr.err2, occErr.err3, occErr.err4, occErr.err5, occErr.err2f, occErr.err3f, occErr.err4f, occErr.err5f);
    kittiStr = sprintf('%s\nDispPix-noc 2/3/4/5 %.3f & %.3f & %.3f & %.3f\nFlowPix-noc 2/3/4/5 %.3f & %.3f & %.3f & %.3f\n', kittiStr, noccErr.err2, noccErr.err3, noccErr.err4, noccErr.err5, noccErr.err2f, noccErr.err3f, noccErr.err4f, noccErr.err5f);
    kittiStr = sprintf('%s\nDispEPE %.3f & %.3f\nFlowEPE %.3f & %.3f\n', kittiStr, epes.epe_nocD, epes.epeD, epes.epe_noc, epes.epe);
    
    fid = fopen(sprintf('%s/RESULTS_K%03d_%02d_%s.txt', p.sFolder, p.imgNr, p.subImg, date), 'w', 'n');
    if fid~=-1
      fwrite(fid, kittiStr, 'char');
      fclose(fid);
    end
    
    if ~exist(sprintf('%s/disp/', p.sFolder), 'dir')
      mkdir(sprintf('%s/disp/', p.sFolder));
    end
    if ~exist(sprintf('%s/flow/', p.sFolder), 'dir')
      mkdir(sprintf('%s/flow/', p.sFolder));
    end

    % store results:
    if storeOutput
      flow_write( cat(3, squeeze(flow2d(:,:,2)), squeeze(flow2d(:,:,3)), ones(size(squeeze(flow2d(:,:,3))))), sprintf('%s/flow/%06d_%02d.png', p.sFolder, p.imgNr, p.subImg ));
      disp_write( squeeze(-flow2d(:,:,1)), sprintf('%s/disp/%06d_%02d.png', p.sFolder, p.imgNr, p.subImg ));
    end

  end
end
end
