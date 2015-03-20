%%% fits solution with a simple iterated (2 its) least squares (L1 distance)
%%% minimizing the reprojection error in 2d. Here left and right cam flow
function [Rt_, l1_error]= algebraicL2_Rot_C(RT, K, M, m, R, q1, q2, p, iDepth, visq1, visq2, l1_error, mprior)

% visq2 : visible points of q2 is non-useless points

% scale ideally respects configuration and xyz stuff aka p2d/dxyz
multMe = (min(1024,numel(p)/3))/256;
if exist('mprior', 'var')
  multMe = mprior * multMe;
end

ss = 0.1;
mid = round(numel(p)/6);
pd = K*p(:,mid);

pxyz1 = pd ./ iDepth(mid);
pdx = (pxyz1+[ss,0,0]');pdx = pdx./pdx(3);
pdy = (pxyz1+[0,ss,0]');pdy = pdy./pdy(3);
pdz = (pxyz1+[0,0,ss]');pdz = pdz./pdz(3);
scalexyz = [norm(pdx-pd), norm(pdy-pd), norm(pdz-pd)];
mat= multMe * diag( cat(2, [1,1,1], (1/ss)./scalexyz));

kp1 = cat(1, zeros(1, size( p,2 )),  p(3,:), -p(2,:) );
kp2 = cat(1, -p(3,:), zeros(1, size( p,2 )),  p(1,:) );
kp3 = cat(1,  p(2,:), -p(1,:),  zeros(1, size( p,2 )) );

KPx = zeros (3, 3*size(kp1,2));
KPx(:,1:3:end) = kp1;
KPx(:,2:3:end) = kp2;
KPx(:,3:3:end) = kp3;

%save('derivativeStuff.mat', 'RT','K','M','m','R','q1','q2','p','iDepth','KPx','L1');

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
% is linear in r and t - no constraint - is it stable ?

RT(1:3) = r; % to make ph1 == ph1_
%%%%%%%%%%%%%%%%%%%%%%%%%
ip   = p;
krkp = K*R*ip;
KRKPx = K*R*KPx;

if sum(visq1)+sum(visq2) > 10 %sum(visq2) >= sum(visq1) &&
  q1 = q1(:,visq1);
  iDepth_ = iDepth(visq1);
  krkp=krkp(:,visq1);
  KRKPx=KRKPx(:,reshape( repmat(visq1', 3,1), [size(KRKPx, 2), 1]));
else
  iDepth_ = iDepth;
end

KRKPxt_ = cat(1, reshape(KRKPx(:, 1:3:end), 1, numel(KRKPx)/3), reshape(KRKPx(:, 2:3:end), 1, numel(KRKPx)/3, 1), reshape(KRKPx(:, 3:3:end), 1, numel(KRKPx)/3));

%or |-q1*krkp3+krkp1 + q1*(KRKPx*r)3 + q1*(K*t)3*iD - (KRKPx*r)1 -(K*t)1*iD)|
rhs = -bsxfun(@times, q1(1:2,:), krkp(3,:)) + krkp(1:2,:);
A_r1 = bsxfun(@times, q1(1,:), KRKPxt_(:,3:3:end) ) - KRKPxt_(:,1:3:end);
A_r2 = bsxfun(@times, q1(2,:), KRKPxt_(:,3:3:end) ) - KRKPxt_(:,2:3:end);

A_r  = reshape( cat(1, A_r1, A_r2), 3, size(A_r1,2)*2 );

A_t1 = bsxfun(@times, K(3:3:end)', iDepth_.* q1(1,:) ) - bsxfun(@times, K(1:3:end)', iDepth_ );
A_t2 = bsxfun(@times, K(3:3:end)', iDepth_.* q1(2,:) ) - bsxfun(@times, K(2:3:end)', iDepth_ );

A_t  = reshape( cat(1, A_t1, A_t2), 3, size(A_t1,2)*2 );

A = cat(2, A_r', A_t');


if any(visq2) || ~any(visq1)

  mrkp  = (M*R)*ip-bsxfun(@times, m, iDepth);
  MRKPx = (M*R)*KPx;

if sum(visq1)+sum(visq2) > 10 % sum(visq2) > sum(visq1) && 
  q2 = q2(:,visq2);
  iDepth_ = iDepth(visq2);
  mrkp=mrkp(:,visq2);
  MRKPx=MRKPx(:,reshape( repmat(visq2', 3,1), [size(MRKPx, 2), 1]));
else
  iDepth_ = iDepth; 
end

  MRKPxt_ = cat(1, reshape(MRKPx(:, 1:3:end), 1, numel(MRKPx)/3), reshape(MRKPx(:, 2:3:end), 1, numel(MRKPx)/3, 1), reshape(MRKPx(:, 3:3:end), 1, numel(MRKPx)/3));

  rhsb =-bsxfun(@times, q2(1:2,:), mrkp(3,:)) + mrkp(1:2,:);
  B_r1 = bsxfun(@times, q2(1,:), MRKPxt_(:,3:3:end) ) - MRKPxt_(:,1:3:end);
  B_r2 = bsxfun(@times, q2(2,:), MRKPxt_(:,3:3:end) ) - MRKPxt_(:,2:3:end);

  B_r  = reshape( cat(1, B_r1, B_r2), 3, size(B_r1,2)*2 );

  B_t1 = bsxfun(@times, M(3:3:end)', iDepth_.* q2(1,:) ) - bsxfun(@times, M(1:3:end)', iDepth_ );
  B_t2 = bsxfun(@times, M(3:3:end)', iDepth_.* q2(2,:) ) - bsxfun(@times, M(2:3:end)', iDepth_ );

  B_t  = reshape( cat(1, B_t1, B_t2), 3, size(B_t1,2)*2 );
  B = cat(2, B_r', B_t');
  
  rhs = cat(1, rhs(:), rhsb(:));
  A = cat(1, A, B);
  
  diagOutW2 =  cat(1, squeeze(mrkp(3,:)), squeeze(mrkp(3,:)) );
else
  diagOutW2 = [];
end

diagOutW1 = -cat(1, squeeze(krkp(3,:)), squeeze(krkp(3,:)) );
diagOutW = 1./cat(1, diagOutW1(:), diagOutW2(:) );%.^2

if exist('l1_error','var')

  fullL1 = max( l1_error(1:2:end), l1_error(2:2:end));

  % half are outlier:
  outlier = ( fullL1 > median(fullL1) );outlier = cat( 1, outlier, outlier );

  diagOut = 1./max(0.1, sqrt(l1_error));
  diagOut(outlier) = 0.001;
  % or (L2)
  diagOut = ones(numel(diagOut),1);
  diagOut(outlier) = 0;  
  B = bsxfun(@times, A, diagOut);

  B   = cat(1,B,   mat );
  rhs = cat(1, diagOut.*rhs(:), zeros(6,1) );
  rt = B \ (rhs(:));
  A=B;% well now its the same -- influences sigma though but it did before as well

  if norm(rt(1:3)) > 0.2
    %fix rt to 0:
    C=B(:,4:6);
    rt = C \ (rhs(:));rt = cat(1, zeros(3,1), rt);
  end
  
else
  A   = cat(1,A,   mat );
  rhs = cat(1, rhs(:), zeros(6,1) );
  rt = A\rhs(:);
  
  if norm(rt(1:3)) > 0.2
    %fix rt to 0:
    C=A(:,4:6);
    rt = C \ rhs(:);rt = cat(1, zeros(3,1), rt);
  end
end

l1_error = abs(A*rt-rhs(:));
l1_error = l1_error(1:end-6);

Rt = -cat(1, rt(4:6), rt(1:3));
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
  breakhere = 1; % fitting produced non-sense
end
