function [I1, I2, flowGT_occ, flowGT_noc, imageName] = ...
  loadKITTIImage_Stereo( imgNr, folder, frame, reduzFactor )
 
  reduzFactor =1;
  
  imgNameGT = sprintf( '%06i_10.png', imgNr);
  imgName1  = sprintf( '%06i_10.png', imgNr);
  %  imgName2 = sprintf( '%06i_10.png', imgNr);
  
  I1 = double(imread( sprintf('%s/image_0/%s', folder, imgName1 ) ));
  I2 = double(imread( sprintf('%s/image_1/%s', folder, imgName1 ) ));
  
  flowGT_noc = disp_read(sprintf('%s/disp_noc/%s', folder, imgNameGT ));
  flowGT_occ = disp_read(sprintf('%s/disp_occ/%s', folder, imgNameGT));
  
  if exist ('reduzFactor', 'var') %&& reduzFactor~=1
    I1 = imresize( I1, 1/reduzFactor );
    I2 = imresize( I2, 1/reduzFactor );
   
    fGt1 = imresize( flowGT_noc(:,:,1), 1/reduzFactor, 'nearest' );
    fGt2 = zeros(size(I1));%imresize( flowGT_noc(:,:,2), 1/reduzFactor );
    fGt3 = round(imresize( double(flowGT_noc>0), 1/reduzFactor, 'nearest' ));
    
    flowGT_noc = cat(3, fGt1./reduzFactor, fGt2, fGt3);
    
    fGt1 = imresize( flowGT_occ(:,:,1), 1/reduzFactor, 'nearest' );
    fGt2 = zeros(size(I1));%fGt2 = imresize( flowGT_occ(:,:,2), 1/reduzFactor );
%    fGt3 = round(imresize( flowGT_occ(:,:,3), 1/reduzFactor ));
    fGt3 = round(imresize( double(flowGT_occ>0), 1/reduzFactor, 'nearest' ));
%    fGt3 = double(flowGT_occ>0);%
    flowGT_occ = cat(3, fGt1./reduzFactor, fGt2, fGt3);
  end
  
imageName = sprintf('KITTI_%03d_%02d', imgNr, frame);
  
