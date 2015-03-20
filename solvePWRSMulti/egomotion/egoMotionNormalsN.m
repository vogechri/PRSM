%%% do gradient descent on whole normal start at N_start,
%%% N_res and RT_res define the desired 2d motion we wnat to fit to
function [N_scale, R_scale] = egoMotionNormalsN(ref, cam, Seg, N_start, N_res, Rt_res, p2d_, Rglob)

nPix = 20;
centers = findPlaneCenter( Seg, 0, N_res );

kt = repmat ( cam.Kl * Rglob(1:3,4), 1, size(N_start, 2));
Mr = cam.Kl * Rglob(1:3,1:3);

h3 = zeros(3, 3, numel(Seg.Ids));
for i = 1:numel(Seg.Ids)

  rt = squeeze(Rt_res(:,:,i));
%  rt = cat( 1, Rglob, [0,0,0,1]);
  rt = [eye(3), centers(:,i); 0,0,0,1] * rt * [eye(3), -centers(:,i); 0,0,0,1];
  
  %  H = rt(1:3,1:3)-(rt(1:3,1:3)'*rt(1:3,4))*N(1:3,i)';
  H = rt(1:3,1:3) -(rt(1:3,4))*N_res(1:3,i)';
  H = cam.Kl * H;% * iKl;
  % H = R{i}-(t{i})*N{i}';
  
%  H = H';% matlab -- cpp
  h3(:,:,i) = H;
end

centers0D =  bsxfun( @rdivide, centers,  centers(3,:));

p2d_x = reshape(p2d_ ,numel(Seg.Img), 3 )';

p2dSAll = zeros(3, size(N_start, 2), nPix);
p2dGAll = zeros(3, size(N_start, 2), nPix);
p2dAAll = zeros(3, size(N_start, 2), nPix);

for i= 1:numel(Seg.Ids)
  ids = Seg.Ids{i}+1;

  num = cat(2, ceil( rand(1,nPix-3)* numel(ids)), 1, numel(ids));
 
  p2dS = cat(2, centers0D(:, i), p2d_x( :,ids(num) ));

  p2dA = Mr * p2dS;

  p2dG = h3(:,:,i) * p2dS;
  p2dG = bsxfun( @rdivide, p2dG, p2dG(3,:) );
  
  p2dSAll(:, i, :) = reshape(p2dS, 3, 1, nPix);
  p2dGAll(:, i, :) = reshape(p2dG, 3, 1, nPix);
  p2dAAll(:, i, :) = reshape(p2dA, 3, 1, nPix);
end

nStart  = N_start(1:3,:);

% highly dependend on amount of steps taken :(
[nRes, fX_] = minimizeVektor(nStart, @derivativeEgoMotionNormalsN, [20 0.1], kt, p2dSAll, p2dAAll, p2dGAll);

N_scale = nRes;

centers = findPlaneCenter( Seg, 0, N_scale );
R_scale = rot_global2local(repmat(cat( 1, Rglob, [0,0,0,1]), [1,1,size(N_scale,2)]), centers);

%getKittiErr3dSF ( Seg, ref,cam, N_scale, R_scale );

return;
