%%% optimized for gray scale images - can not fully utilize 
%%% the infomation of color images. 
%%% a different superpixelation algorithm might yield better results here
%
%%% also beneficial is to use initial depth information here :
function S = SegmentImageCube( Img, ps, par )

if ~exist('ps', 'var') 
  ps     = round(sqrt(8*8*( 1242*375/(624*461)))); % the size of the patches in the segmentation algorithm
  ps = ps +1; % test an enlarged area, worse step1 BUT equal step2 NICE
end

%%%%%%%%%%%%%%%%%%%%%%%%
ISup = Img*255;

S = struct('Img', {}, 'Ids', {}, 'NeighIds', {}, 'NeighIIds', {}, 'Seeds', {}, 'PatchSize', {});
S(1).PatchSize = ps;

%dumbest idea ever: just a grid:
a = 1:ps:size(ISup,1);
b = 1:ps:size(ISup,2);
if a(end) < size(ISup,1)
  a(end+1) = size(ISup,1);
end
if b(end) < size(ISup,2)
  b(end+1) = size(ISup,2);
end
for i=1:numel(a)-1
  for j=1:numel(b)-1
    solution(a(i):a(i+1), b(j):b(j+1)) = (numel(a)-1)*(j-1)+i;
  end
end
solution = uint32(solution);
nSegs = numel(unique(solution));

[S(1).Img, S(1).NeighIds, S(1).Ids, S(1).NeighIIds, S(1).Seeds, S(1).Areas, S(1).Centers, S(1).Edges, S(1).DisSim] = ...
  SuperPixelMisc(uint32(ISup), ps, 100, 5, int32(solution), nSegs); % ok

%plotSegOverlay(Img, S, 400, par, 0.5, sprintf('cube%03d', 300));