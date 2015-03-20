%%% optimized for gray scale images - can not fully utilize 
%%% the infomation of color images. 
%%% a different superpixelation algorithm might yield better results here
%
%%% also beneficial is to use initial depth information here :
function S = SegmentImage( Img, ps, depth )

[M,N,C] = size(Img);
if ~exist('ps', 'var') 
  ps     = round(sqrt(8*8*( 1242*375/(624*461)))); % the size of the patches in the segmentation algorithm
  ps = ps +1; % test an enlarged area, worse step1 BUT equal step2 NICE
end

%%%%%%%%%%%%%%%%%%%%%%%%
ISup = Img*255;

u = ones(M,N,2);
u(:,:,1) = repmat( [1:N],  M, 1 );
u(:,:,2) = repmat( [1:M]', 1, N );

if exist('depth', 'var')
%  disp  = 1./depth;
  disp  = depth;
  disp  = (255/(max(disp(:))-min(disp(:))))   * (disp  - min(disp(:)));
  dispImg = cat(3, ISup, disp );
else
  dispImg = ISup;
end

S = struct('Img', {}, 'Ids', {}, 'NeighIds', {}, 'NeighIIds', {}, 'Seeds', {}, 'PatchSize', {});
S(1).PatchSize = ps;

% more than 1 channel: weight per channel and combine by (addtion/multiplication)?
if exist('depth', 'var')
  [weight_x, weight_y, weight_xy, weight_ixy] = getEdgeWeightsForSegmentation (dispImg(:,:,2)./255, 35, 0.1);
else
  [weight_x, weight_y, weight_xy, weight_ixy] = getEdgeWeightsForSegmentation (dispImg./255, 35, 0.1);
end

minW = min(weight_xy(:));
maxW = max(weight_xy(:));
%ps, iterations, unary, max edge weight (fix 01 and border), nothing, edgeweighs
rw=1000;solution = SuperPixelStandard(uint32(dispImg), ps, 3, 0, 10*rw*maxW, 0, rw*(weight_x-minW), rw*(weight_y-minW), rw*(weight_xy-minW), rw*(weight_ixy-minW));

nSegs = numel(unique(solution));

[S(1).Img, S(1).NeighIds, S(1).Ids, S(1).NeighIIds, S(1).Seeds, S(1).Areas, S(1).Centers, S(1).Edges, S(1).DisSim] = ...
  SuperPixelMisc(uint32(dispImg), ps, 100, 5, int32(solution), nSegs); % ok
