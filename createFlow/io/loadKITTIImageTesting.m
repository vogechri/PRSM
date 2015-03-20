function [I1, I2, imageName] = ...
  loadKITTIImageTesting( imgNr, folder, reduzFactor, doStereo)

if doStereo
  imgName1 = sprintf( '%06i_10.png', imgNr);
  imgName2 = sprintf( '%06i_10.png', imgNr);
  
%  I1 = imread( sprintf('%s/image_0/%s', folder, imgName1 ) );
%  I2 = imread( sprintf('%s/image_0/%s', folder, imgName2 ) );
 
% images at rt+1 
  I1 = imread( sprintf('%s/image_0/%s', folder, imgName1 ) );
  I2 = imread( sprintf('%s/image_1/%s', folder, imgName2 ) );
else
  
  imgName1 = sprintf( '%06i_10.png', imgNr);
  imgName2 = sprintf( '%06i_11.png', imgNr);
  
%  I1 = imread( sprintf('%s/image_0/%s', folder, imgName1 ) );
%  I2 = imread( sprintf('%s/image_0/%s', folder, imgName2 ) );
 
% images at rt+1 
  I1 = imread( sprintf('%s/image_0/%s', folder, imgName1 ) );
  I2 = imread( sprintf('%s/image_0/%s', folder, imgName2 ) );
end

%  flowGT_noc = flow_read(sprintf('%s/flow_noc/%s', folder, imgName1 ));
%  flowGT_occ = flow_read(sprintf('%s/flow_occ/%s', folder, imgName1 ));
  
  if exist ('reduzFactor', 'var') && reduzFactor~=1
    I1 = imresize( I1, 1/reduzFactor );
    I2 = imresize( I2, 1/reduzFactor );    
    flowGT_noc = cat(3, fGt1, fGt2, fGt3);
  end
   
  imageName = sprintf('KITTI_%d', imgNr);
  
