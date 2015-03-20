function colors = plotSegOverlayConsistent_3F(Img, S, Img2, S2, Img3, S3, Img4, S4, Img5, S5, Img6, S6, fnr, par, alpha, fname, colors )

  Img  = Img ./max(max(Img));
  Img2 = Img2./max(max(Img2));
  Img3 = Img3./max(max(Img3));
  Img4 = Img4./max(max(Img4));  
  Img5 = Img5./max(max(Img5));
  Img6 = Img6./max(max(Img6));
  
  [N,M,D] = size(Img);
  
  if ~exist('fnr', 'var')
    fnr = 514;
  end
  
  if ~exist('fname', 'var')
    fname = '';
  end
  
  if ~exist('alpha', 'var')
      alpha = 0.5;
  end

  if ~exist('w', 'var')
      w = 0;
  end

  img  = double(S.Img);
  img2 = double(S2.Img);
  img3 = double(S3.Img);
  img4 = double(S4.Img);
  img5 = double(S5.Img);
  img6 = double(S6.Img); 
  
  if min(img(:)) < 1 || min(img2(:)) < 1|| min(img3(:)) < 1|| min(img4(:)) < 1 || min(img5(:)) < 1|| min(img6(:)) < 1
    img = img +1;
    img2 = img2 +1;
    img3 = img3 +1;
    img4 = img4 +1;
    img5 = img5 +1;
    img6 = img6 +1;    
  end

  if exist('colors','var') && ndims(colors) == 3 
    colR = colors(:,:,1);
    colG = colors(:,:,2);
    colB = colors(:,:,3);
  else  
    maxCol = max( [max(img(:)), max(img2(:)), max(img3(:)), max(img4(:)), max(img5(:)), max(img6(:))] );
    colR = rand(maxCol,1);
    colG = rand(maxCol,1);
    colB = rand(maxCol,1);
    colors = cat(3, colR, colG, colB );
  end
  
%  col = rand(maxCol,1);
  overLay = img;
  overLay(:,:,1) = colR(img);%col = rand(maxCol,1);
  overLay(:,:,2) = colG(img);%col = rand(maxCol,1);
  overLay(:,:,3) = colB(img);
  
  overLay2 = img2;
  overLay2(:,:,1) = colR(img2);%col = rand(maxCol,1);
  overLay2(:,:,2) = colG(img2);%col = rand(maxCol,1);
  overLay2(:,:,3) = colB(img2);  
  overLay3 = img3;
  overLay3(:,:,1) = colR(img3);%col = rand(maxCol,1);
  overLay3(:,:,2) = colG(img3);%col = rand(maxCol,1);
  overLay3(:,:,3) = colB(img3);  
  overLay4 = img4;
  overLay4(:,:,1) = colR(img4);%col = rand(maxCol,1);
  overLay4(:,:,2) = colG(img4);%col = rand(maxCol,1);
  overLay4(:,:,3) = colB(img4);
  overLay5 = img5;
  overLay5(:,:,1) = colR(img5);%col = rand(maxCol,1);
  overLay5(:,:,2) = colG(img5);%col = rand(maxCol,1);
  overLay5(:,:,3) = colB(img5);
  overLay6 = img6;
  overLay6(:,:,1) = colR(img6);%col = rand(maxCol,1);
  overLay6(:,:,2) = colG(img6);%col = rand(maxCol,1);
  overLay6(:,:,3) = colB(img6);  
  
  % less UIJ_AreThereWindowShowsPending - timeout waiting for window to show up
% disable drawings:
  theOverlay = overLay*alpha+(1-alpha)*repmat(mean(Img./max(Img(:)),3), [1,1,3]);
  f   = figure(fnr);  set(f, 'visible','off');imshow(theOverlay);axis off;
  if exist('par', 'var')
    export_fig( sprintf('%s/Segments_%d_overlay_%s_1.png', par.sFolder, par.imgNr, fname), '-m1');
    close(f);
  end
  theOverlay2 = overLay2*alpha+(1-alpha)*repmat(mean(Img2./max(Img2(:)),3), [1,1,3]);
  f   = figure(fnr+1);  set(f, 'visible','off');imshow(theOverlay2);axis off;
  if exist('par', 'var')
    export_fig( sprintf('%s/Segments_%d_overlay_%s_2.png', par.sFolder, par.imgNr, fname), '-m1');
    close(f);
  end 
  theOverlay3 = overLay3*alpha+(1-alpha)*repmat(mean(Img3./max(Img3(:)),3), [1,1,3]);
  f   = figure(fnr+1);  set(f, 'visible','off');imshow(theOverlay3);axis off;
  if exist('par', 'var')
    export_fig( sprintf('%s/Segments_%d_overlay_%s_3.png', par.sFolder, par.imgNr, fname), '-m1');
    close(f);
  end 
  theOverlay4 = overLay4*alpha+(1-alpha)*repmat(mean(Img4./max(Img4(:)),3), [1,1,3]);
  f   = figure(fnr+1);  set(f, 'visible','off');imshow(theOverlay4);axis off;
  if exist('par', 'var')
    export_fig( sprintf('%s/Segments_%d_overlay_%s_4.png', par.sFolder, par.imgNr, fname), '-m1');
    close(f);
  end
  theOverlay5 = overLay5*alpha+(1-alpha)*repmat(mean(Img5./max(Img5(:)),3), [1,1,3]);
  f   = figure(fnr+1);  set(f, 'visible','off');imshow(theOverlay5);axis off;
  if exist('par', 'var')
    export_fig( sprintf('%s/Segments_%d_overlay_%s_5.png', par.sFolder, par.imgNr, fname), '-m1');
    close(f);
  end
  theOverlay6 = overLay6*alpha+(1-alpha)*repmat(mean(Img6./max(Img6(:)),3), [1,1,3]);
  f   = figure(fnr+1);  set(f, 'visible','off');imshow(theOverlay6);axis off;
  if exist('par', 'var')
    export_fig( sprintf('%s/Segments_%d_overlay_%s_6.png', par.sFolder, par.imgNr, fname), '-m1');
    close(f);
  end
  
%  ff  = figure(fnr+2);imshow(overLay);axis on;
%  fff = figure(fnr+3);imshow(overLay2);axis on;  
%  ffff = figure(fnr+4);imshow(overLay3);axis on;  
%  fff = figure(fnr+2); imshow((1-alpha)*repmat(0.025*imgGy+0.025*imgGx+mean(Img,3), [1,1,3])+alpha*overLay), colormap(gray)


%%%%%%%%%%%%%%
%   if w
%     F = getframe(f);
%     imwrite(F.cdata, './SegmentedImage_overlay1.png', 'png');
%     
%     F = getframe(ff);
%     imwrite(F.cdata, './SegmentedImage.png', 'png');
%     
%     F = getframe(fff);
%     imwrite(F.cdata, './SegmentedImageAmplified.png', 'png');
% 
%     if plotEdges 
%       F = getframe(ffff);
%       imwrite(F.cdata, './SegmentedImageedgeMap.png', 'png');
%     end
%   end  

  
  pause(1);% less UIJ_AreThereWindowShowsPending errors ? since then windows really are closed
end
