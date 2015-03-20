% computes the pixel values for the homogra[hy defined by the normals,
% rotation and translation of the references camera
function Homos = get_Homographies_cam(ref, cam, k, N, Rt, time) %, Seg, center)

Kr = cam(k).Kr;
Mc = cam(k).Rr;
mc = cam(k).Tr;
Kl = cam(1).Kl;
iKl = inv(Kl);

Homos = zeros(3,3, size(N,2));
  
for i = 1:size(Rt,3)
  
  if time >= 1
    rt = squeeze(Rt(:,:,i));
    H = Mc*rt(1:3,1:3)-(mc+Mc*rt(1:3,4))*N(1:3,i)';
    H = Kr * H * iKl;
    %    H = Mc*R{i}-(mc+Mc*t{i})*N{i}';
  else
    H = eye(3) - mc*N(1:3,i)';
    H = Kr * H * iKl;
    %    H = Mc-mc*N{i}';
  end
  Homos(:,:,i) = H;
end
