%% computes the segment-segment weights for the smoothness cost
%% also converts the edges and centers into matlab and image plane
%% coordinates [Matlab coordinates are needed since the projection 
%%              matrices assume these coordinates]
function Seg = setWeights_patchSmooth( Seg, Kl )

% should write to Seg.Weights = |e_ij| / [ dist(c_i/e_ij) + dist(c_j/e_ij)]
% later add sismilarity weight, e.g. e(- sim/sigma(patch) )

% loop over 

Kl = inv(Kl);
edges   = Seg.Edges;
Seg.pixEdges = edges;
weights = cell( numel(edges), 1 );
Seg.pixCenters = Seg.Centers;
newCenters = Seg.Centers;

av_edge = 0;
nEdges = 0;

for i=1:numel(edges)

  locEdges = edges{i}+1;%% statr from 1 not from 0
  center = squeeze( Seg.Centers(i,:) )+1;
  locWeights = zeros( size(locEdges ,1), 1);
  
  temp = Kl * [center(2); center(1); 1];% flip
  newCenters(i,:) = [temp(1), temp(2)];
  
  for j=1:size (locEdges ,1)
%    nId = locEdges(j,1);
    
%    pEdge = [locEdges(j,2), locEdges(j,3)];

    % ortho to edge
    o_edge  = [-locEdges(j,3)+locEdges(j,5), locEdges(j,2)-locEdges(j,4)];

    temp = Kl * [locEdges(j,3); locEdges(j,2); 1];% flip
    locEdges(j,2) = temp(1); locEdges(j,3) = temp(2);
    temp = Kl * [locEdges(j,5); locEdges(j,4); 1];% flip
    locEdges(j,4) = temp(1); locEdges(j,5) = temp(2);

    eLength = norm(o_edge);
    % center2 = squeeze( Seg.Centers(nId, :) )+1;
    
    % perpendicular distance:

    % d1 = abs(o_edge * (center  - pEdge)');
    % d2 = abs(o_edge * (center2 - pEdge)');

    % this is mean: truly it is norm(o_edge) / (d1+d2)
    % but d1 = eLength*(o_edge/norm(o_edge)) * perpDist_center_to_edge
    % so we get eLength^2/ 
    % eLength*((perpDist_center1_to_edge) + (perpDist_center2_to_edge)
    %    locWeights(j) = eLength*eLength / (d1+d2);
  
    
    % NEW TEST
    locWeights(j) = eLength / 9;
 
    av_edge = av_edge + eLength / 9;
    nEdges = nEdges +1;
  end
  edges{i}   = locEdges';  % used defines topology
  weights{i} = locWeights; % not used anyway - used is the edge length
end

% must scale smoothness by this here
% aver_edge = av_edge / nEdges;

Seg.Weights = weights;
Seg.Edges   = edges;
Seg.Centers = newCenters;

% scale the weights:
% todo : LATER
%{
scale  = exp(-1);
scaleS = 1-scale;
doDisSim = 1;
if doDisSim
  tes  = Seg.DisSim;% 1: high similarity, 0: low similarity
  for i=1:numel(tes)
    Seg.Weights{i} = Seg.Weights{i} .* ((tes{i}(:,2)) * scaleS + scale) * 2;
  end
end
%}