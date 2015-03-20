%%% computes the pixel values for the homography defined by the normals,
%%% rotation and translation of the references camera
function [homos, viewNormals] = compute_Homos_vNormals(ref, cam, N, Rt, Seg, center)

%here : Mc = Id, mc = 0;
Kl = cam(1).Kl;
Kr = cam(1).Kr;
Mc = cam(1).Rr;
mc = cam(1).Tr;
Kl = cam(1).Kl;
iKl = inv(Kl);

Nhoms = numel(Seg.Ids);

% normal homographies to apply on pixel positions
homos = zeros(9, 3 * Nhoms);

% these normals deliver the depth at a pixel in a different image than the
% reference image
viewNormals = zeros(3, 3 * Nhoms);

% first left,t+1:
% the rotation has to be around the center of the patch !
%
% (Mc|mc) * [(I|c)(R|T) (I|-c)] = (R|c-R*c+T)-> new T, namely
%
% (Mc|mc) * [(I|c)(R|T) (I|-c)] = (R|c-R*c+T)
% -> new T, namely c-R*c+T; new R as usual
% use the flow for finding the mean of the patch

% S(1).Img, S(1).NeighIds, S(1).Ids, S(1).NeighIIds
for i = 1:numel(Seg.Ids)
  rt = squeeze(Rt(:,:,i));
  
  rt = [eye(3), center(:,i); 0,0,0,1] * rt * [eye(3), -center(:,i); 0,0,0,1];
  
  %  H = rt(1:3,1:3)-(rt(1:3,1:3)'*rt(1:3,4))*N(1:3,i)';
  H = rt(1:3,1:3) -(rt(1:3,4))*N(1:3,i)';
  H = Kl * H;% * iKl;
  % H = R{i}-(t{i})*N{i}';
  
%  H = H';% matlab -- cpp
  homos(:,i+Nhoms) = H(:);
  
  Rn = (rt(1:3,1:3) * N(1:3,i));
  Rnt = center(:,i)'*N(1:3,i) - (rt(1:3,4)+center(:,i))' * Rn +1;
  viewNormals(:,i+Nhoms) = Rn ./ Rnt;
end


% right at t=0
for i = 1:numel(Seg.Ids)
  H = eye(3) - mc*N(1:3,i)';
  H = Kr * H;% * iKl;
  %    H = Mc-mc*N{i}';

%  H = H';% matlab -- cpp
  homos(:,i) = H(:);
  
  Rn = (Mc * N(1:3,i));
  Rnt = 1-(mc' * Rn);
  viewNormals(:,i) = Rn ./ Rnt;
end

% right at t=1
for i = 1:numel(Seg.Ids)
  rt = squeeze(Rt(:,:,i));
  rt = [eye(3), center(:,i); 0,0,0,1] * rt * [eye(3), -center(:,i); 0,0,0,1];
  
  %    H = Mc*rt(1:3,1:3)-(mc+Mc*rt(1:3,1:3)'*rt(1:3,4))*N(1:3,i)';
  H = Mc*rt(1:3,1:3)-(mc+Mc*rt(1:3,4))*N(1:3,i)';
  H = Kr * H;% * iKl;
  %    H = Mc*R{i}-(mc+Mc*t{i})*N{i}';

%  H = H';% matlab -- cpp
  homos(:,i+2*Nhoms) = H(:);
  
  Rn = (Mc * viewNormals(1:3,i));
  Rnt = 1-(mc' * Rn);
  viewNormals(:,i+2*Nhoms) = Rn ./ Rnt;
end
