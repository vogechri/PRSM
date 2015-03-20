%%% algebraci fitting with subsequent gradient descent. cost function
%%% lorentzian applied on the reprojection error
function [Rt_, l1_error]= algebraicL2_Rot(RT, K, M, m, R, q1, q2, p, iDepth, l1_error, sigmaS, gdits)

if ~exist('L1', 'var')
  L1 = 1;
end

if ~exist('sigmaS', 'var')
  sigmaS = 0.5;
end

if ~exist('gdits', 'var')
  gdits = 10;
end

usealgebraic = 0;

if usealgebraic
kp1 = cat(1, zeros(1, size( p,2 )),  p(3,:), -p(2,:) );
kp2 = cat(1, -p(3,:), zeros(1, size( p,2 )),  p(1,:) );
kp3 = cat(1,  p(2,:), -p(1,:),  zeros(1, size( p,2 )) );

KPx = zeros (3, 3*size(kp1,2));
KPx(:,1:3:end) = kp1;
KPx(:,2:3:end) = kp2;
KPx(:,3:3:end) = kp3;

r = RT(1:3);
t = RT(4:6);
R_cross = [0, -r(3), r(2); r(3), 0, -r(1); -r(2), r(1), 0];

sinA = norm(r);
alpha = asin(sinA);
cosA = cos(alpha);

R1 = zeros(3);
R1(1,1:3:end) = cosA;
R1(2,2:3:end) = cosA;
R1(3,3:3:end) = cosA;

R2 = R_cross;
r = r / max(sinA, eps);
R3=(1-cosA) * (r*r');

Rot = R1 + R2 + R3;
R   = R*Rot;


r = zeros(3,1);
%ph1 = -krkp + KRKPx * r + K*t*iD
%min (q1-ph1/ph3)^2
%or |q1*ph3 - ph1| = |q1*[-krkp3 + (KRKPx * r)3 + (K*t)3*iD] - (-krkp1+ (KRKPx * r)1 + (K*t)1*iD)|
%or |-q1*krkp3+krkp1 + q1*(KRKPx*r)3 + q1*(K*t)3*iD - (KRKPx*r)1 -(K*t)1*iD)|
% is linear in r and t - no constraint stable enough

RT(1:3) = r; % to make ph1 == ph1_
%%%%%%%%%%%%%%%%%%%%%%%%%
ip   = p;
krkp = K*R*ip;
dt   = bsxfun(@times, iDepth, t );
Kdt  = K * dt;

KRKPx = K*R*KPx;

KRKPxt_ = cat(1, reshape(KRKPx(:, 1:3:end), 1, numel(KRKPx)/3), reshape(KRKPx(:, 2:3:end), 1, numel(KRKPx)/3, 1), reshape(KRKPx(:, 3:3:end), 1, numel(KRKPx)/3));

am1t = cat( 1, KRKPxt_(:,1:3:end), bsxfun(@times, K(1:3:end)', iDepth) );
am2t = cat( 1, KRKPxt_(:,2:3:end), bsxfun(@times, K(2:3:end)', iDepth) );
am3t = cat( 1, KRKPxt_(:,3:3:end), bsxfun(@times, K(3:3:end)', iDepth) );

if ~isempty(q2)
  mrkp  = (M*R)*ip-bsxfun(@times, m, iDepth);
  MRKPx = (M*R)*KPx;
  MRKPxt = cat(2, reshape(MRKPx(:, 1:3:end), numel(MRKPx)/3, 1), reshape(MRKPx(:, 2:3:end), numel(MRKPx)/3, 1), reshape(MRKPx(:, 3:3:end), numel(MRKPx)/3, 1));
  bm1t = cat( 1, MRKPxt(1:3:end,:)', bsxfun(@times, M(1:3:end)', iDepth) );
  bm2t = cat( 1, MRKPxt(1:3:end,:)', bsxfun(@times, M(2:3:end)', iDepth) );
  bm3t = cat( 1, MRKPxt(1:3:end,:)', bsxfun(@times, M(3:3:end)', iDepth) );
  Bm = cat(3, bm1t,bm2t,bm3t );
  Mdt = M * dt;
  ph2 = -mrkp + reshape(MRKPxt*r, 3, numel(iDepth)) + Mdt;
end

%or |-q1*krkp3+krkp1 + q1*(KRKPx*r)3 + q1*(K*t)3*iD - (KRKPx*r)1 -(K*t)1*iD)|
rhs = -bsxfun(@times, q1(1:2,:), krkp(3,:)) + krkp(1:2,:);
A_r1 = bsxfun(@times, q1(1,:), KRKPxt_(:,3:3:end) ) - KRKPxt_(:,1:3:end);
A_r2 = bsxfun(@times, q1(2,:), KRKPxt_(:,3:3:end) ) - KRKPxt_(:,2:3:end);

A_r  = reshape( cat(1, A_r1, A_r2), 3, size(A_r1,2)*2 );

A_t1 = bsxfun(@times, K(3:3:end)', iDepth.* q1(1,:) ) - bsxfun(@times, K(1:3:end)', iDepth );
A_t2 = bsxfun(@times, K(3:3:end)', iDepth.* q1(2,:) ) - bsxfun(@times, K(2:3:end)', iDepth );

A_t  = reshape( cat(1, A_t1, A_t2), 3, size(A_t1,2)*2 );

A = cat(2, A_r', A_t');

if exist('l1_error','var')
 
  x_out = l1_error(1:end/2) > median(l1_error(1:end/2));
  l1_error ( l1_error > median(l1_error) ) = 1000*1000;
  
  D = diag(1./max(0.001, sqrt(l1_error)));
  B=A'*D;
  rt = (B*B') \ (B*D*rhs(:));
  
  if norm(rt(1:3)) > 0.2
    %fix rt to 0:
    C=B(4:6,:);
    rt = (C*C') \ C*rhs(:);rt = cat(1, zeros(3,1), rt);
  end
 
else
  rt = (A'*A) \ A'*rhs(:);
  
  if norm(rt(1:3)) > 0.2
    %fix rt to 0:
    C=A(:,4:6);
    rt = (C'*C) \ C'*rhs(:);rt = cat(1, zeros(3,1), rt);
  end

end
l1_error = abs(A*rt-rhs(:));

[F, dF] = fitHomo_MexWindows(-rt, K, M, m, R, q1(1:2,:), p, iDepth);

%or |-q1*krkp3+krkp1 + q1*(KRKPx*r)3 + q1*(K*t)3*iD - (KRKPx*r)1 -(K*t)1*iD)|
% note that:
%ph1_ = -krkp + cat(2, Am(:,:,1)' * RT, Am(:,:,2)' * RT, Am(:,:,3)' * RT )';
%ph2_ = -mrkp + cat(2, Bm(:,:,1)' * RT, Bm(:,:,2)' * RT, Bm(:,:,3)' * RT )';

err=A*rt-rhs(:);rx=err(1:2:end);ry=err(2:2:end);
avS = sum( sqrt(rx.^2 + ry.^2) ) / numel(rx);
else
  rt = RT;
  avS=sum(sqrt(l1_error(1:2:end).^2+l1_error(2:2:end).^2)) / numel(iDepth);
end

sigmaM = sigmaS * avS;

% lorentzian
step = 1; % highly dependent on this - sigh
if isempty(q2)
  [x, fX_] = minimize(-rt, @fitHomoR, [gdits 0.05], K, M, m, R, q1(1:2,1:step:end), p(:,1:step:end), iDepth(1:step:end), 0, sigmaM);%, 0, q2(1:2,1:step:end));
else
  [x, fX_] = minimize(-rt, @fitHomoR, [gdits 0.05], K, M, m, R, q1(1:2,1:step:end), p(:,1:step:end), iDepth(1:step:end), 0, sigmaM, q2(1:2,1:step:end));
end

Rt = cat(1, x(4:6), x(1:3));
%fprintf('Score: %f %f, - translation %.3f,%.3f,%.3f - rotation %.3f,%.3f,%.3f\n', fX_(1), fX_(end), Rt(1:6) );

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

if ~isreal(Rt_)
  breakhere = 1;
end
