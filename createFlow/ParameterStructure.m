function parameters = ParameterStructure( theta, lambda, warps, ...
  pyramid_factor, sFolder, infoFileName, imgNr, reduzFactor, inneritsD, ...
  inneritsS, use_edges, stt_factor, pyr_Smooth, useGradientImages, doStereo )
 
if ~exist ('stt_factor', 'var')
  stt_factor = 0.75;
end
if ~exist ('pyr_Smooth', 'var')
  pyr_Smooth = 0;
end
if ~exist ('useGradientImages', 'var')
  useGradientImages = 1;
end
if ~exist ('use_edges', 'var')
  use_edges = 0;
end
if ~exist ('inneritsD', 'var')
  inneritsD = 3;
end
if ~exist ('inneritsS', 'var')
  inneritsS = 25;
end
if ~exist ('reduzFactor', 'var')
  reduzFactor = 1;
end

if ~exist ('doStereo', 'var')
  doStereo = 0;
end

p.doStereo = doStereo; % stereo or flow

p.maxits = 1; % always use that: warp instead: costs nothing (only on gpu)

% construct and use a gradient image for the data term as well (add also weight later)
p.useGradientImages = useGradientImages;
p.wGradient         = 1;  % weight of gradient image

p.reduzFactor = reduzFactor;

% preprocessing parameter:
p.use_structure_texture    = 1;   % on  by default, off for Hoguet ?!
p.structure_texture_sigma  = 1.0; % 1.0 default
p.structure_texture_factor = stt_factor;% combination of texture and structure term
%
p.theta       = 3.0; % weight coupling the auxxiliary variables
p.lambda      = 10;    % weight of the data term
p.lambdaF     = 0;     % not used
p.warps       = 8;     % # of warps per level
p.innerItsD   = 3;     % # iterations of data term minimizer per warp
p.innerItsS   = 25;    % # iterations of smooth term minimizer per warp
p.use_edges   = use_edges; % use gradients of reference view as smoothing prior

p.pyramid_factor = 0.75;
p.doMedian       = 1;     % median filtering at each warp

% default
% boxes off others: on?
p.stockI = 0; % build pyramid images from previous level, not from finest
p.gaussI = 0; % smooth images in pyramid with gaussian filter
if pyr_Smooth == 1
  p.stockI = 1; % build pyramid images from previous level, not from finest
  p.gaussI = 0; % smooth images in pyramid with gaussian filter
elseif pyr_Smooth ==2
  p.stockI = 1; % build pyramid images from previous level, not from finest
  p.gaussI = 1; % smooth images in pyramid with gaussian filter
end

p.storeFolder              = './scratch/test'; 
p.sFolder                  = './scratch/test';
p.infoFileName             = 'Test.inf';
p.imgNr                    = 11;

p.writeIntermediateResults = 0; % write intermediate flow results (2d and 3d)
p.writeIntermediateImages  = 1; % write intermediate results : images

% TV norm with quadratic behavior between 0 and epsilon
p.epsilon = 0.1;
p.lambdaF = 0;

if exist( 'sFolder', 'var')
  while(sFolder(end) == '/')
    sFolder = sFolder(1:end-1);
  end
  p.sFolder = sFolder;
  p.storeFolder = sFolder;
end
if exist( 'infoFileName', 'var') 
  idx = findstr(infoFileName, '.inf');
  if ~isempty(idx) == 1
    infoFileName = infoFileName(1:idx-1);
  end
  p.infoFileName = infoFileName;
end
if exist( 'imgNr', 'var')
  p.imgNr = imgNr;
end
if exist( 'pyramid_factor', 'var')
  p.pyramid_factor = pyramid_factor;
end
if exist( 'use_edges', 'var')
  p.use_edges = use_edges;
  p.use_structure_texture = use_edges;
end
if exist( 'inneritsS', 'var')
  p.innerItsS= inneritsS;
end
if exist( 'inneritsD', 'var')
  p.innerItsD= inneritsD;
end
if exist( 'warps', 'var')
  p.warps = warps;
end
% if exist( 'lambdaF', 'var')
%   p.lambdaF = lambdaF;
% end
if exist( 'lambda', 'var')
  p.lambda = lambda;
end
if exist( 'theta', 'var')
  p.theta = theta;
end
if exist( 'storeFolder', 'var')
  while(storeFolder(end) == '/')
    storeFolder = storeFolder(1:end-1);
  end
  p.storeFolder = storeFolder;
end

if exist( 'stt_factor', 'var')
  p.structure_texture_factor = stt_factor;
end

parameters = p;

end


% weighting of the data term
%lambda  = 30;%16.0; %  35 to 50 ok lambda in [10.0, 1000.0]
%lambdaF = 2.000001;   % lambda for the fundamental matrix data term

% quadratic relaxation parameter
%theta = 0.5; % theta in [0.1, 0.5] 0.25 standard
% 0.5 when doing use_structure_texture = 2;
% that is with gauaussian blur
%0.5 strong - 0.1 weak smoothing
% 0.5 sharp edges
% larger theta needs higher innerits ? NO

% TV norm with quadratic behavior between 0 and epsilon
%epsilon = 0.1;

% use edge weighted TV
%use_edges = 0; % 0 = no, 1 = yes
% a better formulation might be to punish the edge weight between the
% different pixels, -> anisotropic smoothing
% high edge weight means large consistency in the flow,
% low edge weight means low consistency
% sum over the neighbours of u_i
% 0 = sum (w_ij * (u_j - u_i)) <=> u_i = 1/sum(w_ij) * sum(w_ij * u_j)
% update fixpoint equation: u_i = u_i + lambda * 1/sum(w_ij) * sum(w_ij * u_j)
% w_ij small where high gradient in image!

% number of TV-L1 iterations per warp
%maxits = 1; % 1
% number of inner TV-L1 iterations (smoothing steps) per warp
%innerits = 25;

% warping parameters
%max_pyramid_levels = 1000; % max number of pyramid levels
%pyramid_factor = 0.95; % in [0.5, 0.95]
%warps = 5; % number of warps per level

% preprocessing
%use_structure_texture = 1;
%structure_texture_sigma = 1.0;
% test for best parameter on the venus set
%0.55 : 0.266533;%0.58 : 0.270450;%0.6 : 0.276;%0.565: 0.269;%0.5 : 0.295;%0.63 : 0.288;
%structure_texture_factor = 0.82; % for enhanced , that is with l1
%neighborhood difference as data term
