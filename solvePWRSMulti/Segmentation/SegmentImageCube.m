%%% optimized for gray scale images - can not fully utilize 
%%% the infomation of color images. 
%%% a different superpixelation algorithm might yield better results here.
%%% Super pixels are computed by SuperPixelMisc, with no source or
%%% interface available. That function takes an initial segmentation, which
%%% defaults to fixed size squares. Optionally, this initial segmentation
%%% can be further split by passing regions in the left and right images
%%% that should be considered static (such as bodies rigidly attached to
%%% the camera) or invalid (such as rectification artifacts) in
%%% par.leftStatic, par.rightStatic, par.leftValid, par.rightValid 
%%% All images can be logical masks of image size. par.leftStatic can also
%%% be a double image, where each non-zero pixel would contain a prior 
%%% estimate of disparity. For more information on using a double 
%%% par.leftValid, see generateProposals. In this function, it is cast 
%%% down to a logical.
%
%%% also beneficial is to use initial depth information here :
%
% see also generateProposals

function S = SegmentImageCube( Img, ps, par)

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
% Split segments that include two classes of pixels (valid/invalid or
% static/dynamic view) by moving all masked pixels to a new range and then
% re-numbering
if isfield(par,'leftValid') 
    solution = splitLabelImageWMask(solution, ~par.leftValid);
end
if isfield(par,'rightValid')
    solution = splitLabelImageWMask(solution, ~par.rightValid);
end
if isfield(par,'leftStatic')
    solution = splitLabelImageWMask(solution, logical(par.leftStatic));
end
if isfield(par,'rightStatic')
    solution = splitLabelImageWMask(solution, logical(par.rightStatic));
end
nSegs = numel(unique(solution));

[S(1).Img, S(1).NeighIds, S(1).Ids, S(1).NeighIIds, S(1).Seeds, S(1).Areas, S(1).Centers, S(1).Edges, S(1).DisSim] = ...
  SuperPixelMisc(uint32(ISup), ps, 100, 5, int32(solution), nSegs); % ok
end
function labelImage = splitLabelImageWMask(labelImage, mask)
% move all pixels in mask to new range, then re-index 
    labelImage(mask) = (max(labelImage(:))+1)+labelImage(mask);
    [~,~,labelImage(:)] = unique(labelImage);
end