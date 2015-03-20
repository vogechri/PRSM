function Rt_lin = fullRotationFit(r_start, SegIds, cam, Rt_start, q2d_2, p2d_1, iDepth, sigmas, RT_its, r2d_2 )

if ~exist('RT_its','var')
  RT_its = 10;
end

Rt_lin = zeros(4,4,numel(SegIds));

if exist('r2d_2','var')
  % this is correct, but am if vis 2d not visible it sucks:
  [r_end, fX_] = minimizeVektor(r_start, @fitHomoAllR, [RT_its 0.1], SegIds, cam.Kl, cam.Kr*cam.Rr, cam.Kr*cam.Tr, Rt_start, q2d_2(1:2,:), p2d_1, iDepth, 0, 0.5*sigmas, r2d_2(1:2,:));
else
  [r_end, fX_] = minimizeVektor(r_start, @fitHomoAllR, [RT_its 0.1], SegIds, cam.Kl, cam.Kr*cam.Rr, cam.Kr*cam.Tr, Rt_start, q2d_2(1:2,:), p2d_1, iDepth, 0, 0.5*sigmas);
end

for k=1:size(r_start, 2)
  R = Rt_start(:,:,k);
  Rt = cat(1, r_end(4:6,k), r_end(1:3,k));
  R_cross = [0, -Rt(6), Rt(5); Rt(6), 0, -Rt(4); -Rt(5), Rt(4), 0];

  r = Rt(4:6);
  sinA = norm(r);
  alpha = asin(sinA);
  cosA = cos(alpha);

  R1 = zeros(3);
  R1(1,1:3:end) = cosA;
  R1(2,2:3:end) = cosA;
  R1(3,3:3:end) = cosA;

  R2 = R_cross;
  R3=(1-cosA) * (r*r')./max( eps, sinA^2);

  Rot = R1 + R2 + R3;
  Rt_ = [R*Rot,Rt(1:3)];
  Rt_lin(:,:,k) = cat(1, Rt_, [0,0,0,1]);
end
