%%% computes the pixel values for the homogra[hy defined by the normals,
%%% rotation and translation of the references camera
function cam = Homography_cam_pixel_Center(ref, cam, k, N, Rt, Seg, time, center )

Kr = cam(k).Kr;
Mc = cam(k).Rr;
mc = cam(k).Tr;
Kl = cam(1).Kl;
iKl = inv(Kl);

[m,n] = size(cam(k).q2dx);

p2d = ones(m,n,3);
p2d(:,:,1) = repmat( [1:n],  m, 1 );
p2d(:,:,2) = repmat( [1:m]', 1, n );

pix = zeros(size(p2d));
NM = size(p2d,1) * size(p2d,2);

% S(1).Img, S(1).NeighIds, S(1).Ids, S(1).NeighIIds
if (size(N,2) > 1) || (time ~= 1)
  
for i = 1:min( numel(Seg.Ids), size(Rt,3))
  if time > 1
    rt = squeeze(Rt(:,:,i));
    rt = [eye(3), center(:,i); 0,0,0,1] * rt * [eye(3), -center(:,i); 0,0,0,1];

%    H = Mc*rt(1:3,1:3)-(mc+Mc*rt(1:3,1:3)'*rt(1:3,4))*N(1:3,i)';
    H = Mc*rt(1:3,1:3)-(mc+Mc*rt(1:3,4))*N(1:3,i)';
    H = Kr * H * iKl;
    %    H = Mc*R{i}-(mc+Mc*t{i})*N{i}';
  else
    H = eye(3) - mc*N(1:3,i)';
    H = Kr * H * iKl;
    %    H = Mc-mc*N{i}';
  end
  
  p2dx = p2d(1+Seg.Ids{i});
  p2dy = p2d(1+Seg.Ids{i}+   NM);
  p2dz = p2d(1+Seg.Ids{i}+ 2*NM);
  
  pix(1+Seg.Ids{i})       = H(1,1) * p2dx + H(1,2) * p2dy + H(1,3) * p2dz;
  pix(1+Seg.Ids{i}+ NM)   = H(2,1) * p2dx + H(2,2) * p2dy + H(2,3) * p2dz;
  pix(1+Seg.Ids{i}+ 2*NM) = H(3,1) * p2dx + H(3,2) * p2dy + H(3,3) * p2dz;
 
end
else

  H = eye(3) - mc*N(1:3,1)';
  H = Kr * H * iKl;

  p2dx = p2d(:,:,1);
  p2dy = p2d(:,:,2);
  p2dz = p2d(:,:,3);
  
  pix(:,:,1) = H(1,1) * p2dx + H(1,2) * p2dy + H(1,3) * p2dz;
  pix(:,:,2) = H(2,1) * p2dx + H(2,2) * p2dy + H(2,3) * p2dz;
  pix(:,:,3) = H(3,1) * p2dx + H(3,2) * p2dy + H(3,3) * p2dz;
  
end

c3  = squeeze(pix(:,:,3));
pix = bsxfun (@rdivide, pix, pix(:,:,3));
cam(k).I(time).u = pix;
cam(k).I(time).u(:,:,3) = c3;
