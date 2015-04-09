%
%   Duality-based variational TGV Optical Flow (and Stereo)
%
%   Author: Christoph Vogel
%

function [flow_out, flow_out2] = pyramidFlowDouble ( I1, I2, I3, I4, p, startLev, init_uv )

[M, N, C] = size(I1);
I1 = scale_image(double(I1), 0, 1);
I2 = scale_image(double(I2), 0, 1);
I3 = scale_image(double(I3), 0, 1);
I4 = scale_image(double(I4), 0, 1);

if ~exist('startLev', 'var')
  startLev=80;% 28 is fine for kitti
end

width_Pyramid  = cell(1,1);
height_Pyramid = cell(1,1);

%N*p.pyramid_factor^levelX = 28;
levelX = log(startLev/N)/log(p.pyramid_factor);
levelY = log(startLev/M)/log(p.pyramid_factor);
levels = round(min(levelY, levelX));
pyramid_factor_x = exp(log(startLev/N)/levels);
pyramid_factor_y = exp(log(startLev/M)/levels);

% precalculate image sizes
pyramid_levels = levels+1;
width_Pyramid{1} = N;
height_Pyramid{1} = M;
for i = 2:pyramid_levels
  width_Pyramid{i}  = round(pyramid_factor_x*width_Pyramid{i-1});
  height_Pyramid{i} = round(pyramid_factor_y*height_Pyramid{i-1});
  if max(width_Pyramid{i}, height_Pyramid{i}) < startLev
    pyramid_levels = i;
    break;
  end
end

  % set up image pyramides
  I1_Pyramid = cell(pyramid_levels,1);
  I2_Pyramid = cell(pyramid_levels,1);
  I3_Pyramid = cell(pyramid_levels,1);
  I4_Pyramid = cell(pyramid_levels,1);
  
  for i = 1:pyramid_levels
    images = cat ( 3, imresize(I1, [height_Pyramid{i} width_Pyramid{i}], 'bicubic'), ...
                      imresize(I2, [height_Pyramid{i} width_Pyramid{i}], 'bicubic'), ...
                      imresize(I3, [height_Pyramid{i} width_Pyramid{i}], 'bicubic'), ...
                      imresize(I4, [height_Pyramid{i} width_Pyramid{i}], 'bicubic'));

    I1_Pyramid{i} = images(:,:,1);
    I2_Pyramid{i} = images(:,:,2);
    I3_Pyramid{i} = images(:,:,3);
    I4_Pyramid{i} = images(:,:,4);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start coarse to fine processing
for level = pyramid_levels:-1:1;
  
  M = height_Pyramid{level};
  N = width_Pyramid{level};
  
  if level == pyramid_levels
    
    % initialization
    u = zeros(M,N);
    v = zeros(M,N);
    
    if exist('init_uv', 'var') && numel(init_uv) >= N*M*2
      scale = N/size(init_uv, 2);
      u = imresize ( init_uv(:,:,1), [M, N] );
      v = imresize ( init_uv(:,:,2), [M, N] );
      u = u * scale;
      v = v * scale;
    end

    w   = zeros(M,N,8);
    pp  = zeros(M,N,8);
    qq  = zeros(M,N,16);

    u(:,:,2)=u;v(:,:,2)=v;

  else

    rescale_factor_u = width_Pyramid{level+1}/width_Pyramid{level};
    rescale_factor_v = height_Pyramid{level+1}/height_Pyramid{level};

    % prolongate to finer grid
%     u = imresize(u,[M N], 'bilinear')/rescale_factor_u;
%     v = imresize(v,[M N], 'bilinear')/rescale_factor_v;

    pp_tmp = u;
    u = zeros(M,N, size( pp_tmp, 3 ));
    for i=1:size( pp_tmp, 3 )
      u(:,:,i) = imresize(pp_tmp(:,:,i), [M N], 'bilinear')/rescale_factor_u;
    end
    pp_tmp = v;
    v = zeros(M,N, size( pp_tmp, 3 ));
    for i=1:size( pp_tmp, 3 )
      v(:,:,i) = imresize(pp_tmp(:,:,i), [M N], 'bilinear')/rescale_factor_v;
    end
    pp_tmp = pp;
    pp = zeros(M,N, size( pp_tmp, 3 ));
    for i=1:size( pp_tmp, 3 )
      pp(:,:,i) = imresize(pp_tmp(:,:,i), [M N], 'bilinear');
    end
    pp_tmp = qq;
    qq = zeros(M,N, size( pp_tmp, 3 ));
    for i=1:size( pp_tmp, 3 )
      qq(:,:,i) = imresize(pp_tmp(:,:,i), [M N], 'bilinear');
    end
    pp_tmp = w;
    w  = zeros(M,N, size( pp_tmp, 3 ));
    for i=1:size( pp_tmp, 3 )
      w(:,:,i) = imresize(pp_tmp(:,:,i), [M N], 'bilinear');
    end
    clear('pp_tmp');
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [u, v, w, pp, qq] = tgv_level_flowDouble(I1_Pyramid{level}, I2_Pyramid{level}, ...
                                            I3_Pyramid{level}, I4_Pyramid{level}, ...
                                            u, v, w, pp, qq, p );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% store final flow
flow_out = zeros(M,N,2);
flow_out(:,:,1) = u(:,:,1);
flow_out(:,:,2) = v(:,:,1);
flow_out2(:,:,1) = u(:,:,2);
flow_out2(:,:,2) = v(:,:,2);
end
