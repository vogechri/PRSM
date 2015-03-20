%%% computes the pixel values for the homogra[hy defined by the normals, 
%%% rotation and translation of the references camera
function ref = Homography_ref_pixel_Center(ref, cam, time, N, Rt, Seg, center )

%here : Mc = Id, mc = 0;
%Kr = cam(1).Kr;
%Mc = cam(1).Rr;
%mc = cam(1).Tr;

Kl = cam(1).Kl;
iKl = inv(Kl);

[m,n] = size(cam(1).q2dx);

p2d = ones(m,n,3);
p2d(:,:,1) = repmat( [1:n],  m, 1 );
p2d(:,:,2) = repmat( [1:m]', 1, n );

if time == 1
  ref.I(time).u = p2d;
  return;
end

pix = zeros(size(p2d));
NM = size(p2d,1) * size(p2d,2);

% the rotation has to be around the center of the patch !
% 
% (Mc|mc) * [(I|c)(R|T) (I|-c)] = (R|c-R*c+T)-> new T, namely 
%
% (Mc|mc) * [(I|c)(R|T) (I|-c)] = (R|c-R*c+T)
% -> new T, namely c-R*c+T; new R as usual
% use the flow for finding the mean of the patch

% S(1).Img, S(1).NeighIds, S(1).Ids, S(1).NeighIIds
for i = 1:min( size(Rt,3), numel(Seg.Ids))
  rt = squeeze(Rt(:,:,i));
  rt = [eye(3), center(:,i); 0,0,0,1] * rt * [eye(3), -center(:,i); 0,0,0,1];

  % H = R{i}-(t{i})*N{i}';
  H = rt(1:3,1:3) -(rt(1:3,4))*N(1:3,i)';
  H = Kl * H * iKl;

  pix(1+Seg.Ids{i})       = H(1,1) * p2d(1+Seg.Ids{i}) + H(1,2) * p2d(1+Seg.Ids{i}+  NM) + H(1,3) * p2d(1+Seg.Ids{i}+ 2*NM);
  pix(1+Seg.Ids{i}+ NM)   = H(2,1) * p2d(1+Seg.Ids{i}) + H(2,2) * p2d(1+Seg.Ids{i}+  NM) + H(2,3) * p2d(1+Seg.Ids{i}+ 2*NM);
  pix(1+Seg.Ids{i}+ 2*NM) = H(3,1) * p2d(1+Seg.Ids{i}) + H(3,2) * p2d(1+Seg.Ids{i}+  NM) + H(3,3) * p2d(1+Seg.Ids{i}+ 2*NM);
end

c3 = pix(:,:,3);
pix = bsxfun (@rdivide, pix ,pix(:,:,3));
ref.I(time).u = pix;
ref.I(time).u(:,:,3) = c3;