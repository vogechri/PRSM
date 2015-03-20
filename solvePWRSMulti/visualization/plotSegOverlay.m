function plotSegOverlay(Img, S, fnr, par, alpha, fname, w )

  plotEdges = 0;
  plotSeeds = 0;
  plotWeights = 0;
  Img = Img./max(max(Img));
  
  [N,M,D] = size(Img);

  if ~isfield(S, 'pixEdges') && ~isfield(S, 'Edges')
    plotEdges = 0;
  end
  
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

  img = double(S.Img);
  
  add = 0;
  if min(img(:)) < 1
    add = 1;
    img = img +1;
  end    

  maxCol = max(img(:));
  colR = rand(maxCol,1);
  colG = rand(maxCol,1);
  colB = rand(maxCol,1);
  
%  col = rand(maxCol,1);
  overLay = img;
  overLay(:,:,1) = colR(img);%col = rand(maxCol,1);
  overLay(:,:,2) = colG(img);%col = rand(maxCol,1);
  overLay(:,:,3) = colB(img);
% less UIJ_AreThereWindowShowsPending - timeout waiting for window to show up
% disable drawings:
  theOverlay = overLay*alpha+(1-alpha)*repmat(mean(Img./max(Img(:)),3), [1,1,3]);
  f   = figure(fnr);  imshow(overLay*alpha+(1-alpha)*repmat(mean(Img./max(Img(:)),3), [1,1,3]));axis off;
  if exist('par', 'var')
    path(path,'C:\Users\test1\Documents\MATLAB\export_fig');
    export_fig( sprintf('%s/Segments_overlay_%s_%d.png', par.sFolder, fname, par.imgNr), '-m1');
    close(f);
  end

  
  ff  = figure(fnr+1);imshow(overLay);axis on;


%  fff = figure(fnr+2); imshow((1-alpha)*repmat(0.025*imgGy+0.025*imgGx+mean(Img,3), [1,1,3])+alpha*overLay), colormap(gray)

%%%%%
if isfield(S, 'NeighIIds') && isfield(S, 'NeighIds') && 0

  overLay = repmat(mean(Img./max(Img(:)),3), [1,1,3]);
  overLay = theOverlay;
  nIds = S.NeighIIds;
  mIds = S.NeighIds;
  for j=1:numel( nIds )
    nIIds = nIds{j};
    mIIds = mIds{j}+add;

    for i=1:numel( nIIds )
      index = (nIIds{i})+add;
      overLay( index )       = colR(mIIds(i));
      overLay( index +N*M)   = colG(mIIds(i));
      overLay( index +N*M*2) = colB(mIIds(i));
    end
  end
  fff  = figure(fnr+2); imshow(overLay, 'Border','tight'), axis image,axis off;
  if exist('par', 'var')
    path(path,'C:\Users\test1\Documents\MATLAB\export_fig');
    export_fig( sprintf('%s/SegmentsO_%s_%d.png', par.sFolder, fname, par.imgNr), '-m2');
    close(fff);
  end
% tmp = getframe(f);tmp=tmp.cdata;
% close(fff);f = figure(fNr+17);set(f, 'visible','off');
% imshow(tmp), axis off;
% export_fig( sprintf('%s/wZ_HQ%s.png', par.sFolder, ending), '-m1');
% close(f);
end
%%%%%

if isfield(S, 'NeighIIds') && isfield(S, 'DisSim') && isfield(S, 'Weights') && plotWeights

  doDisSim = 0;
  
  fullTes = [];
  if doDisSim
    tes  = S.DisSim;% 1: high similarity, 0: low similarity
    tesW = S.Weights;
    for i=1:numel(tes)
      tes{i}(:,2) = tes{i}(:,2) .* tesW{i};
      fullTes = cat(1, fullTes, tes{i}(:,2));
    end
  else
    tes = S.Weights;
    for i=1:numel(tes)
%      tes{i} = 1./tes{i};
      fullTes = cat(1, fullTes, tes{i}(:,1));
      tes{i} = cat(2, tes{i}(:,1), tes{i}(:,1));
    end
  end
%  fullTes = 1./fullTes;
  maxDis = prctile(fullTes, 95);
  minDis = prctile(fullTes,  5);
  for i=1:numel(tes)
%    tes{i}(:,2) = 1-min(1, max(0, (tes{i}(:,2)-minDis)./maxDis));
    tes{i}(:,2) = min(1, max(0, (tes{i}(:,2)-minDis)./maxDis));
  end

  overLay = repmat(mean(Img./max(Img(:)),3), [1,1,3]);
  nIds = S.NeighIIds;
  mIds = S.NeighIds;
  dis  = tes;
  for j=1:numel( nIds )
    nIIds = nIds{j};
    mIIds = mIds{j}+add;
    ds = dis{j};
    
    for i=1:numel( nIIds )
      index = (nIIds{i})+add;
      d     = ds(i,2);
      rgb = getRGB_interpolation(d);
      overLay( index )       = rgb(1);
      overLay( index +N*M)   = rgb(2);
      overLay( index +N*M*2) = rgb(3);
    end
  end
  fff  = figure(fnr+2); imshow(overLay, 'Border','tight'), axis image,axis off;
  if exist('par', 'var')
    path(path,'C:\Users\test1\Documents\MATLAB\export_fig');
    export_fig( sprintf('%s/Segments_Sim_%s_%d.png', par.sFolder, fname, par.imgNr), '-m1');
    close(fff);
  end
end
%%%%%%%%%%%%%%

if plotEdges 
%  overLay = img;
% if isfield(S, 'pixEdges')
%   edges = S.pixEdges;
% else
  if isfield(S,'pixEdges')
    edges = S.pixEdges;
  else
    edges = S.Edges;
  end
  ffff  = figure(fnr+3); imshow(Img);axis off;hold on;
  for i=1:numel(edges)
    locEdge = edges{i}+add;
    for j=1:size( locEdge, 1 )
      plot ([locEdge(j,3), locEdge(j,5)], [locEdge(j,2), locEdge(j,4)], 'r-');
    end
  end
end

  if w
    F = getframe(f);
    imwrite(F.cdata, './SegmentedImage_overlay.png', 'png');
    
    F = getframe(ff);
    imwrite(F.cdata, './SegmentedImage.png', 'png');
    
    F = getframe(fff);
    imwrite(F.cdata, './SegmentedImageAmplified.png', 'png');

    if plotEdges 
      F = getframe(ffff);
      imwrite(F.cdata, './SegmentedImageedgeMap.png', 'png');
    end
  end  
  if plotSeeds
    % check Seeds : ok
    overLay = repmat(mean(Img,3), [1,1,3]);
    marks = S.Img+1;
    mark = zeros(N,M);
    mark( N * round(S.pixCenters(:,2)) + round(S.pixCenters(:,1)) ) = 1;
%    mark(S.Seeds) = 1;
    marks(logical(mark)) = numel(S.Seeds)+1;
    % marks(S.Seeds) = [1:numel(S.Seeds)];
    %  marks(S.Seeds) = numel(S.Seeds)+1;
    maxCol = max(marks(:));
    col = rand(maxCol,1);
    overLay(:,:,1) = col(marks);col = rand(maxCol,1);
    overLay(:,:,2) = col(marks);col = rand(maxCol,1);
    overLay(:,:,3) = col(marks);
    figure(fnr+4), imagesc(overLay);
  end
  
  pause(1);% less UIJ_AreThereWindowShowsPending errors ? since then windows really are closed
end
