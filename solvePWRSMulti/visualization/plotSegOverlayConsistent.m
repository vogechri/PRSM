function colors = plotSegOverlayConsistent(Img, S, Img2, S2, Img3, S3, Img4, S4, fnr, par, alpha, fname, colors )


  Img  = imadjust(Img);
  Img2 = imadjust(Img2);
  Img3 = imadjust(Img3);
  Img4 = imadjust(Img4);
  brightness = 2; % if one wants to see something
  doEdges = 1;
  Img  = Img ./max(max(Img));
  Img2 = Img2./max(max(Img2));
  Img3 = Img3./max(max(Img3));
  Img4 = Img4./max(max(Img4));  
  
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

  add =1;
  if min(img(:)) < 1 || min(img2(:)) < 1|| min(img3(:)) < 1|| min(img4(:)) < 1
    img = img +1;
    img2 = img2 +1;
    img3 = img3 +1;
    img4 = img4 +1;
    add =1;
  end
  
  if exist('colors','var') && ndims(colors) == 3 
    colR = colors(:,:,1);
    colG = colors(:,:,2);
    colB = colors(:,:,3);
  else
    % just a lot to make sure all props get a color !!!! not sgements 
    maxCol = 10*max( [max(img(:)), max(img2(:)), max(img3(:)), max(img4(:))] );
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

  
  % less UIJ_AreThereWindowShowsPending - timeout waiting for window to show up
% disable drawings:
  theOverlay = overLay*alpha+(1-alpha)*repmat(mean(Img./max(Img(:)),3), [1,1,3]) * brightness;
%%%%% kirchen glas effekt
if isfield(S, 'NeighIIds') && isfield(S, 'NeighIds') && doEdges
%  overLay = repmat(mean(Img./max(Img(:)),3), [1,1,3]);
  overLay = theOverlay;
  nIds = S.NeighIIds;
  mIds = S.NeighIds;
  
  ids = img (N * round(S.pixCenters(:,2)) + round(S.pixCenters(:,1)));
  theOverlay( N * round(S.pixCenters(:,2)) + round(S.pixCenters(:,1)) ) = colR( ids );
  theOverlay( N * round(S.pixCenters(:,2)) + round(S.pixCenters(:,1))+N*M ) = colG( ids );
  theOverlay( N * round(S.pixCenters(:,2)) + round(S.pixCenters(:,1))+N*M*2 ) = colB( ids );  
%{  
  for j=1:numel( nIds )
    nIIds = nIds{j};
    mIIds = mIds{j}+add;

    for i=1:numel( nIIds )
      index = (nIIds{i})+add;
      theOverlay( index )       = colR(mIIds(i));
      theOverlay( index +N*M)   = colG(mIIds(i));
      theOverlay( index +N*M*2) = colB(mIIds(i));
    end
  end
%}
end
%%%%%  
  
  f   = figure(fnr);  set(f, 'visible','off'); imshow(theOverlay);axis off;
  if exist('par', 'var')
    export_fig( sprintf('%s/Segments_%d_overlay_%s.png', par.sFolder, par.imgNr, fname), '-m1');
    close(f);
  end
  theOverlay2 = overLay2*alpha+(1-alpha)*repmat(mean(Img2./max(Img2(:)),3), [1,1,3]) * brightness;
%%%%% kirchen glas effekt
if isfield(S2, 'NeighIIds') && isfield(S2, 'NeighIds') && doEdges
%  overLay = repmat(mean(Img./max(Img(:)),3), [1,1,3]);
  overLay = theOverlay2;
  nIds = S2.NeighIIds;
  mIds = S2.NeighIds;
  ids = img2 (N * round(S2.pixCenters(:,2)) + round(S2.pixCenters(:,1)));
  theOverlay2( N * round(S2.pixCenters(:,2)) + round(S2.pixCenters(:,1)) ) = colR( ids );
  theOverlay2( N * round(S2.pixCenters(:,2)) + round(S2.pixCenters(:,1))+N*M ) = colG( ids );
  theOverlay2( N * round(S2.pixCenters(:,2)) + round(S2.pixCenters(:,1))+N*M*2 ) = colB( ids );  
  %{
  for j=1:numel( nIds )
    nIIds = nIds{j};
    mIIds = mIds{j}+add;

    for i=1:numel( nIIds )
      index = (nIIds{i})+add;
      theOverlay2( index )       = colR(mIIds(i));
      theOverlay2( index +N*M)   = colG(mIIds(i));
      theOverlay2( index +N*M*2) = colB(mIIds(i));
    end
  end
  %}
end
%%%%%  
  
  f   = figure(fnr+1); set(f, 'visible','off'); imshow(theOverlay2);axis off;
  if exist('par', 'var')
    export_fig( sprintf('%s/Segments_%d_overlay2_%s.png', par.sFolder, par.imgNr, fname), '-m1');
    close(f);
  end 
  
  theOverlay3 = overLay3*alpha+(1-alpha)*repmat(mean(Img3./max(Img3(:)),3), [1,1,3]) * brightness;
  
%%%%% kirchen glas effekt
if isfield(S3, 'NeighIIds') && isfield(S3, 'NeighIds') && doEdges
%  overLay = repmat(mean(Img./max(Img(:)),3), [1,1,3]);
  overLay = theOverlay3;
  nIds = S3.NeighIIds;
  mIds = S3.NeighIds;
  ids = img3 (N * round(S3.pixCenters(:,2)) + round(S3.pixCenters(:,1)));  
  theOverlay3( N * round(S3.pixCenters(:,2)) + round(S3.pixCenters(:,1)) ) = colR( ids );    
  theOverlay3( N * round(S3.pixCenters(:,2)) + round(S3.pixCenters(:,1))+N*M ) = colG( ids );    
  theOverlay3( N * round(S3.pixCenters(:,2)) + round(S3.pixCenters(:,1))+N*M*2 ) = colB( ids );      
%{
  for j=1:numel( nIds )
    nIIds = nIds{j};
    mIIds = mIds{j}+add;

    for i=1:numel( nIIds )
      index = (nIIds{i})+add;
      theOverlay3( index )       = colR(mIIds(i));
      theOverlay3( index +N*M)   = colG(mIIds(i));
      theOverlay3( index +N*M*2) = colB(mIIds(i));
    end
  end
  %}
end
%%%%%  
  f   = figure(fnr+1); set(f, 'visible','off'); imshow(theOverlay3);axis off;
  if exist('par', 'var')
    export_fig( sprintf('%s/Segments_%d_overlay3_%s.png', par.sFolder, par.imgNr, fname), '-m1');
    close(f);
  end 
  
  theOverlay4 = overLay4*alpha+(1-alpha)*repmat(mean(Img4./max(Img4(:)),3), [1,1,3]) * brightness;
%%%%% kirchen glas effekt
if isfield(S4, 'NeighIIds') && isfield(S4, 'NeighIds') && doEdges
%  overLay = repmat(mean(Img./max(Img(:)),3), [1,1,3]);
  overLay = theOverlay4;
  nIds = S4.NeighIIds;
  mIds = S4.NeighIds;
  ids = img4 (N * round(S4.pixCenters(:,2)) + round(S4.pixCenters(:,1)));  
  theOverlay4( N * round(S4.pixCenters(:,2)) + round(S4.pixCenters(:,1)) ) = colR( ids );  
  theOverlay4( N * round(S4.pixCenters(:,2)) + round(S4.pixCenters(:,1))+N*M ) = colG( ids );  
  theOverlay4( N * round(S4.pixCenters(:,2)) + round(S4.pixCenters(:,1))+N*M*2 ) = colB( ids );    
  %{
  for j=1:numel( nIds )
    nIIds = nIds{j};
    mIIds = mIds{j}+add;

    for i=1:numel( nIIds )
      index = (nIIds{i})+add;
      theOverlay4( index )       = colR(mIIds(i));
      theOverlay4( index +N*M)   = colG(mIIds(i));
      theOverlay4( index +N*M*2) = colB(mIIds(i));
    end
  end
  %}
end
%%%%%  
  f   = figure(fnr+1); set(f, 'visible','off'); imshow(theOverlay4);axis off;
  if exist('par', 'var')
    export_fig( sprintf('%s/Segments_%d_overlay4_%s.png', par.sFolder, par.imgNr, fname), '-m1');
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
