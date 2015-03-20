%%% computed edgeweights between segments and pixels based on the input
%%% images and edges therein
function ew = getEdgeWeights(ref, cam, sW, useThreeFrames, pastFr)

[weight_x , weight_y , weight_xy , weight_ixy ] = getEdgeWeightsForSegmentation (ref.I(1).I, 35, 0.1);
[weight_x2, weight_y2, weight_xy2, weight_ixy2] = getEdgeWeightsForSegmentation (cam(1).I(1).I, 35, 0.1);
[weight_x3, weight_y3, weight_xy3, weight_ixy3] = getEdgeWeightsForSegmentation (ref.I(2).I, 35, 0.1);
[weight_x4, weight_y4, weight_xy4, weight_ixy4] = getEdgeWeightsForSegmentation (cam(1).I(2).I, 35, 0.1);

if numel(cam) ==1
  xWeights  = sW*cat(3, weight_x,   weight_x2,   weight_x3,   weight_x4);
  yWeights  = sW*cat(3, weight_y,   weight_y2,   weight_y3,   weight_y4);
  xyWeights = sW*cat(3, weight_xy,  weight_xy2,  weight_xy3,  weight_xy4);
  yxWeights = sW*cat(3, weight_ixy, weight_ixy2, weight_ixy3, weight_ixy4);

else
  [weight_x5, weight_y5, weight_xy5, weight_ixy5] = getEdgeWeightsForSegmentation (cam(2).I(2).I, 35, 0.1);
  [weight_x6, weight_y6, weight_xy6, weight_ixy6] = getEdgeWeightsForSegmentation (cam(2).I(2).I, 35, 0.1);
  
  xWeights  = sW*cat(3, weight_x,   weight_x2,   weight_x5,   weight_x3,   weight_x4,   weight_x6);
  yWeights  = sW*cat(3, weight_y,   weight_y2,   weight_y5,   weight_y3,   weight_y4,   weight_y6);
  xyWeights = sW*cat(3, weight_xy,  weight_xy2,  weight_xy5,  weight_xy3,  weight_xy4,  weight_xy6);
  yxWeights = sW*cat(3, weight_ixy, weight_ixy2, weight_ixy5, weight_ixy3, weight_ixy4, weight_ixy6);
end

%
if useThreeFrames && exist ('pastFr', 'var')

  [weight_x5, weight_y5, weight_xy5, weight_ixy5] = getEdgeWeightsForSegmentation (pastFr{1}, 35, 0.1);
  [weight_x6, weight_y6, weight_xy6, weight_ixy6] = getEdgeWeightsForSegmentation (pastFr{2}, 35, 0.1);

  xWeights  = cat(3, xWeights,  sW*cat(3, weight_x5,   weight_x6));
  yWeights  = cat(3, yWeights,  sW*cat(3, weight_y5,   weight_y6));
  xyWeights = cat(3, xyWeights, sW*cat(3, weight_xy5,  weight_xy6));
  yxWeights = cat(3, yxWeights, sW*cat(3, weight_ixy5, weight_ixy6));
end

ew.x  = xWeights;
ew.y  = yWeights;
ew.xy = xyWeights;
ew.yx = yxWeights;