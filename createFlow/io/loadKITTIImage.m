function [I1, I2, flowGT_occ, flowGT_noc, imageName] = ...
  loadKITTIImage( imgNr, folder, frame, left )
 
  reduzFactor = 1;

  imgName1 = sprintf( '%06i_10.png', imgNr);
  imgName2 = sprintf( '%06i_11.png', imgNr);
  
  I1 = imread( sprintf('%s/image_0/%s', folder, imgName1 ) );
  I2 = imread( sprintf('%s/image_0/%s', folder, imgName2 ) );
  
  flowGT_noc = flow_read(sprintf('%s/flow_noc/%s', folder, imgName1 ));
  flowGT_occ = flow_read(sprintf('%s/flow_occ/%s', folder, imgName1 ));
  
  if exist ('reduzFactor', 'var') && reduzFactor~=1
    I1 = imresize( I1, 1/reduzFactor );
    I2 = imresize( I2, 1/reduzFactor );
   
    fGt1 = imresize( flowGT_noc(:,:,1), 1/reduzFactor, 'nearest' );
    fGt2 = imresize( flowGT_noc(:,:,2), 1/reduzFactor, 'nearest' );
    fGt3 = round(imresize( flowGT_noc(:,:,3), 1/reduzFactor, 'nearest' ));
    
    flowGT_noc = cat(3, fGt1./reduzFactor, fGt2./reduzFactor, fGt3);
    
    fGt1 = imresize( flowGT_occ(:,:,1), 1/reduzFactor, 'nearest' );
    fGt2 = imresize( flowGT_occ(:,:,2), 1/reduzFactor, 'nearest' );
    fGt3 = round(imresize( flowGT_occ(:,:,3), 1/reduzFactor, 'nearest' ));
    
    flowGT_occ = cat(3, fGt1./reduzFactor, fGt2./reduzFactor, fGt3);
  end
  
%  mask = zeros( size(I1) );
  
imageName = sprintf('KITTI_%03d_%02d', imgNr, frame);
  
