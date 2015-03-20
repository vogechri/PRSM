%% moves centers of segments into segment if necessary
%% distinguishes between matlab (starting at 1,1) and cpp (starting at 0,0)
%% entries
%% after the segmentation the centers and edges are given in cpp style
%% later however (see e.g. bicubic/bilinear interpolation)
%% always the matlab rhing is used -> CHAOS
function Seg = correctSegCenters(Seg, matlabPix)

if ~exist('matlabPix', 'var')
  matlabPix = 0;% input is cpp centers starting at 0,0
end

if isfield(Seg, 'pixCenters')
  pc = Seg.pixCenters;
else
  pc = Seg.Centers;
end

if ~matlabPix
  pc = pc+1;
end

%% step 1 find centers out of segment:
fail= [];
for i = 1:numel(Seg.Ids)
  a=round(pc(i,:));
  if Seg.Img(a(1),a(2))~= i-1
    fprintf('Seg %d\n',i-1);
    fail(end+1)=i-1;
  end;
end
% for i = 1:numel(fail)
%   failSet = find( Seg.Img == fail(i));
%   [I,J]=ind2sub( size(Seg.Img), failSet );
%
%   cx=mean(I)-1;
%   cy=mean(J)-1;
%   cc= Seg.pixCenters(fail(i)+1,:);
%
%   fprintf('Centers: %f %f vs %f %f \n', cx,cy, cc(1), cc(2));
% end

%% step 2 for those not in the segment find better places
% pick an interior point closest to the center

for i = 1:numel(fail)
  
  ids = Seg.Ids{fail(i)+1}+1; % cpp ids start at 0
  
  nids = Seg.NeighIds{fail(i)+1}+1; % cpp ids start at 0
  % find border ids:
  borderIds = [];
  for j=1:numel(nids)
    nnids = Seg.NeighIds{nids(j)}+1; % cpp ids start at 0
    id = find(nnids == fail(i)+1); % cpp ids start at 0
    
    pixIds = Seg.NeighIIds{nids(j)}{id}+1; % cpp ids start at 0
    borderIds=cat(1, borderIds, pixIds);
  end
  
  innerPix = setdiff(ids, borderIds);
  if isempty(innerPix)
    innerPix = borderIds;
  end
  
  [cy,cx] = ind2sub(size(Seg.Img), innerPix);
  
  dist = double(cat(2, cy,cx)) - repmat( pc(fail(i)+1,:), numel(cx), 1);
  [~, mid]=min(sum( dist.^2, 2));
  
  newCenter = cat(2, cy(mid), cx(mid));
  
  cc = pc(fail(i)+1,:); % cpp ids start at 0
  fprintf('Matlab Centers: %f %f vs %f %f \n', newCenter(1), newCenter(2), cc(1), cc(2));
  if matlabPix
    pc(fail(i)+1,:) = newCenter; % back to matlab centers
  else
    pc(fail(i)+1,:) = newCenter-1; % back to cpp centers
  end
end

if isfield(Seg, 'pixCenters')
  for i = 1:numel(fail)
    Seg.pixCenters(fail(i)+1,:) = pc(fail(i)+1,:);
  end
else
  for i = 1:numel(fail)
    Seg.Centers(fail(i)+1,:) = pc(fail(i)+1,:);
  end
end
