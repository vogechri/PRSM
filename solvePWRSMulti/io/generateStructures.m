%%% generate camera matrices and further information used later 
function [data_level, ref_level] = generateStructures ( data_cam, ImgL, Rl )

Icam_Pyramid = data_cam(1).I;
[M, N, ~] = size(ImgL{1});

I_struct = struct( 'I', {}, 'I_w', {}, 'I_wx', {}, 'I_wy', {}, ...
  'u', {}, 'du_d', {}, 'du_wx', {}, 'du_wy', {},'du_wz', {}, 'valid', {}, 'oMap', {} );

data_level = struct( 'I', {}, 'R', {}, 'F', {}, 'Ft', {}, 'epi', {}, 'pop' ,{}, 'popRef', {} );

%Viewport = [ (N)/2 0 (N+1)/2 ; 0 (M)/2 (M+1)/2 ; 0  0  1 ];
I_struct(1).I     =  ImgL{1};
I_struct(2).I     =  ImgL{2};
ref_level.I       = I_struct;
ref_level.R       = data_cam.R;%Viewport * Rl;

for j = 1:size(data_cam,2)
  for k=1:numel(ImgL)
    I_struct(k).I =  Icam_Pyramid{k};
  end
  
  data_level(j).I = I_struct;

  data_level(j).R   = data_cam(j).R;
  data_level(j).epi = data_cam(j).epi;
  data_level(j).pop     = data_cam(j).pop;
  data_level(j).popRef  = data_cam(j).popRef;

  data_level(j).Kl = data_cam(j).Kl;
  data_level(j).Kr = data_cam(j).Kr;
  data_level(j).Rl = data_cam(j).Rl;
  data_level(j).Rr = data_cam(j).Rr;
  data_level(j).Tl = data_cam(j).Tl;
  data_level(j).Tr = data_cam(j).Tr;
end

% constants used without the lower levels
[u, lambda_uv] = Compute_lambdaUV(data_level(j), N, M);

for j = 1:size(data_level,1)
  % direction of derivative does not change with d_:
  [~, dur_d] = Compute_ur_all_d(ref_level.R, data_level(j).R, data_level(j).epi, u, lambda_uv, 10);
  q2dx       = dur_d(:,:,1) ./ sqrt(dur_d(:,:,1).^2 + dur_d(:,:,2).^2);
  data_level(j).q2dx = q2dx;
end
