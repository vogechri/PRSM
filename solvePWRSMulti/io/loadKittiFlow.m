function [cam, ref, imageName, flow2DGt, flow2DGt_noc] = ...
  loadKittiFlow(strFolder, nr, p)

% loading the gt or some other stuff - no need for that only for painting

testing = p.testing;
ImgL      = cell(1,2*p.frames+2);
ImgR      = cell(1,2*p.frames+2);

strNumber0 = sprintf( '%06d_%02d', nr, p.subImg );
strNumber1 = sprintf( '%06d_%02d', nr, p.subImg+1 );
imageName  = sprintf('Kitti_%s', strNumber0);

localFolder = strFolder;
if ~testing

sprintf('%s/%s/%s.png', localFolder, '/data_stereo_flow/training/image_0', strNumber0 )

  ImgL{1} = double(imread(sprintf('%s/%s/%s.png', localFolder, '/data_stereo_flow/training/image_0', strNumber0 )))/255;
  ImgR{1} = double(imread(sprintf('%s/%s/%s.png', localFolder, '/data_stereo_flow/training/image_1', strNumber0 )))/255;
  ImgL{2} = double(imread(sprintf('%s/%s/%s.png', localFolder, '/data_stereo_flow/training/image_0', strNumber1 )))/255;
  ImgR{2} = double(imread(sprintf('%s/%s/%s.png', localFolder, '/data_stereo_flow/training/image_1', strNumber1 )))/255;
  calibFile = sprintf('%s/%s/%06d.txt', localFolder, '/data_stereo_flow/training/calib', nr );
else
  ImgL{1} = double(imread(sprintf('%s/%s/%s.png', localFolder, '/data_stereo_flow/testing/image_0', strNumber0 )))/255;
  ImgR{1} = double(imread(sprintf('%s/%s/%s.png', localFolder, '/data_stereo_flow/testing/image_1', strNumber0 )))/255;
  ImgL{2} = double(imread(sprintf('%s/%s/%s.png', localFolder, '/data_stereo_flow/testing/image_0', strNumber1 )))/255;
  ImgR{2} = double(imread(sprintf('%s/%s/%s.png', localFolder, '/data_stereo_flow/testing/image_1', strNumber1 )))/255;
  calibFile = sprintf('%s/%s/%06d.txt', localFolder, '/data_stereo_flow/testing/calib', nr );
end

data_supp = loadCalibrationKITTI_new( calibFile );
data_supp.I = ImgR;

gt_Name = sprintf('%06d_10', nr);

if testing
  DispImg = ones(size(ImgL{1}));
  flowImg = repmat(DispImg, [1,1,3]);

else
  stereoFolder = sprintf('%s/data_stereo_flow/training/disp_noc', strFolder );
  flowFolder   = sprintf('%s/data_stereo_flow/training/flow_noc', strFolder );
  DispImg = disp_read (sprintf('%s/%s.png',stereoFolder,gt_Name));
  flowImg = flow_read (sprintf('%s/%s.png',flowFolder,gt_Name));
end
flow2DGt_noc = cat (3, DispImg, flowImg);

if testing

else
  stereoFolder = sprintf('%s/data_stereo_flow/training/disp_occ', strFolder );
  flowFolder   = sprintf('%s/data_stereo_flow/training/flow_occ', strFolder );
  DispImg = disp_read (sprintf('%s/%s.png',stereoFolder,gt_Name));
  flowImg = flow_read (sprintf('%s/%s.png',flowFolder,gt_Name));
end

flow2DGt     = cat (3, DispImg, flowImg);
[cam, ref] = generateStructures ( data_supp, ImgL, data_supp.R );

if p.subImg>0
  strNumberprev = sprintf( '%06d_%02d', nr, p.subImg-1 );
  try
    ImgOld = cell(1,2*p.frames+2);
    ImgOld{1} = double(imread(sprintf('%s/%s/%s.png', localFolder, '/data_stereo_flow/training/image_0', strNumberprev )))/255;
    ImgOld{2} = double(imread(sprintf('%s/%s/%s.png', localFolder, '/data_stereo_flow/training/image_1', strNumberprev )))/255;
    cam.Iold = ImgOld;
  catch
    fprintf('Failed to load images: %s', sprintf('%s/%s/%s.png', localFolder, '/data_stereo_flow/training/image_0', strNumberprev ) );
    fprintf('Failed to load images: %s', sprintf('%s/%s/%s.png', localFolder, '/data_stereo_flow/training/image_1', strNumberprev ) );
  end
end
if p.subImg<19
  strNumberprev = sprintf( '%06d_%02d', nr, p.subImg+2 );
  try
  Imgnew = cell(1,2*p.frames+2);
  Imgnew{1} = double(imread(sprintf('%s/%s/%s.png', localFolder, '/data_stereo_flow/training/image_0', strNumberprev )))/255;
  Imgnew{2} = double(imread(sprintf('%s/%s/%s.png', localFolder, '/data_stereo_flow/training/image_1', strNumberprev )))/255;
  cam.Inew = Imgnew;
  catch
    fprintf('Failed to load images: %s', sprintf('%s/%s/%s.png', localFolder, '/data_stereo_flow/training/image_0', strNumberprev ) );
    fprintf('Failed to load images: %s', sprintf('%s/%s/%s.png', localFolder, '/data_stereo_flow/training/image_1', strNumberprev ) );
  end  
end