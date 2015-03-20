% returns the gradients of the scalar field u
% these are defined on the dual elements of the grid
% which is here (linear basis elements) triangles
% note that there are roughly twice as many triangles than vertices 
function [dx, dy, d_xy, dixy] = getEdgeWeightsForSegmentation (u, varWindow, maxDiff)

minVal    = 0.05;
% simple version:

oldWay =1;
if oldWay

g = edge_Tensor(u);

dx = g(:,:,1);dx = dx(:,1:end-1);
dy = g(:,:,4);dy = dy(1:end-1,:);

d_xy = (g(:,:,1) + g(:,:,4) + g(:,:,2) + g(:,:,3))./2;d_xy = d_xy(1:end-1,1:end-1);
dixy = (g(:,:,1) + g(:,:,4) - g(:,:,2) - g(:,:,3))./2;dixy = dixy(1:end-1,1:end-1);

dx = max( minVal, dx );
dy = max( minVal, dy );

dixy = max( minVal, dixy );
d_xy = max( minVal, d_xy );

else

dx   = edge_Tensor_Dir(u, 0.5, 0,   [1,0]);
dy   = edge_Tensor_Dir(u, 0.0, 0.5, [0,1]);
d_xy = edge_Tensor_Dir(u, 0.5, 0.5, [1./sqrt(2),1./sqrt(2)]);
dixy = edge_Tensor_Dir(u, 0.5, 0.5, [1./sqrt(2),-1./sqrt(2)]);

dx = max( minVal, dx );
dy = max( minVal, dy );

dixy = max( minVal, dixy );
d_xy = max( minVal, d_xy );

end

% image sizes not yet fitting for my purpose
return;


if ~exist('maxDiff', 'var')
  maxDiff = 50;
end

if ~exist('varWindow', 'var')
  maxDiff = 35;
end

useScale = 1;
upperScale = 1.0; % the higher the more distinctive - 1.5 ok

smoothVar = 0.01;
%minVal    = 0.015; % worked okish, but seems a bit too small
minVal    = 0.05;

% some consts
%b = 0.95;a = 0.01 * 255^b;
%b = 0.75;a = 0.025 * 255^b;
%b = 0.55;a=100/255; % unger

% set empirically:
b = 1; % b=1:L1, b=2:L2, etc. the metric used 
a = 10;% a is important in the sense that this is the ratio by which the edges are operating - too small: prefers squares, too high: 

b = 1;
a = a^b;

% dy: y - y+1
dy = abs(-u(2:end, :) + u(1:end-1, :));
%dx: x - x+1
dx = abs(-u(:, 2:end) + u(:, 1:end-1));

d_xy = abs(-u(2:end, 2:end) + u(1:end-1, 1:end-1));
dixy = abs(-u(1:end-1, 2:end) + u(2:end, 1:end-1));


dx = 255*min(dx, maxDiff);
dy = 255*min(dy, maxDiff);
d_xy = 255*min(d_xy, maxDiff);
dixy = 255*min(dixy, maxDiff);


if varWindow>0
  var_I = getLocalVariance(u, varWindow);
  var_I = double(var_I);
else
  var_I = ones(size (u));
end

var_I = max(0, var_I);
var_I = sqrt(var_I); % not really important ?

if exist('var_I', 'var')
  var_Iy   = 0.5 * (var_I(1:end-1,:) + var_I(2:end,:));
  var_Ix   = 0.5 * (var_I(:,1:end-1) + var_I(:,2:end));
  var_I_xy = 0.5 * (var_I(1:end-1,1:end-1) + var_I(2:end,2:end));
  var_Iixy = 0.5 * (var_I(2:end,1:end-1) + var_I(1:end-1,2:end));

  % should be scaled based on both dx and dy and set to 0.98 percentiles or
  % so
  
  if useScale
    dx = scale_image(dx./(var_Ix+smoothVar), 0.0, upperScale);
    dy = scale_image(dy./(var_Iy+smoothVar), 0.0, upperScale);
    d_xy = scale_image(d_xy./(var_I_xy+smoothVar), 0.0, upperScale);
    dixy = scale_image(dixy./(var_Iixy+smoothVar), 0.0, upperScale);
  else
    dx = (dx./(var_Ix+smoothVar));
    dy = (dy./(var_Iy+smoothVar));
    d_xy = (d_xy./(var_I_xy+smoothVar));
    dixy = (dixy./(var_Iixy+smoothVar));
  end
else
  dy = scale_image(dy, 0.0, 1.0);
  dx = scale_image(dx, 0.0, 1.0);
end

dx = exp( - a*dx.^b );
dy = exp( - a*dy.^b );
d_xy = exp( - a*d_xy.^b );
dixy = exp( - a*dixy.^b );

tmp = cat(1, dx(:), dy(:));
scale = 1 / mean(tmp(:));
dx = max( minVal, dx .* scale);
dy = max( minVal, dy .* scale);

tmp = cat(1, d_xy(:), dixy(:));
scale = 1 / mean(tmp(:));
dixy = max( minVal, dixy .* scale);
d_xy = max( minVal, d_xy .* scale);
