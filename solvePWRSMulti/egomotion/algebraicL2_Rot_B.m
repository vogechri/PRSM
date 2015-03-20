%%% fits solution with a simple iterated (2 its) least squares (L1 distance)
%%% minimizing the reprojection error in 2d. 
function [Rt_, l1_error]= algebraicL2_Rot_B(RT, K, M, m, R, q1, q2, p, iDepth, l1_error, evalOnly, mprior)

% way to strict at larger segs
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

if ~exist('evalOnly', 'var')
  evalOnly = 0;
end

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
% is linear in r and t - no constraint - is it stable apparently

RT(1:3) = r; % to make ph1 == ph1_
%%%%%%%%%%%%%%%%%%%%%%%%%
ip    = p;
krkp  = K*R*ip;
dt    = bsxfun(@times, iDepth, t );
%Kdt   = K * dt;
KRKPx = K*R*KPx;

% the transposed krkpx:
KRKPxt_ = cat(1, reshape(KRKPx(:, 1:3:end), 1, numel(KRKPx)/3), reshape(KRKPx(:, 2:3:end), 1, numel(KRKPx)/3, 1), reshape(KRKPx(:, 3:3:end), 1, numel(KRKPx)/3));

if ~isempty(q2)
  mrkp  = (M*R)*ip-bsxfun(@times, m, iDepth);
  MRKPx = (M*R)*KPx;
  MRKPxt = cat(2, reshape(MRKPx(:, 1:3:end), numel(MRKPx)/3, 1), reshape(MRKPx(:, 2:3:end), numel(MRKPx)/3, 1), reshape(MRKPx(:, 3:3:end), numel(MRKPx)/3, 1));
  bm1t = cat( 1, MRKPxt(1:3:end,:)', bsxfun(@times, M(1:3:end)', iDepth) );
  bm2t = cat( 1, MRKPxt(1:3:end,:)', bsxfun(@times, M(2:3:end)', iDepth) );
  bm3t = cat( 1, MRKPxt(1:3:end,:)', bsxfun(@times, M(3:3:end)', iDepth) );
  %Bm = cat(3, bm1,bm2,bm3 );
  Bm = cat(3, bm1t,bm2t,bm3t );
  Mdt = M * dt;
  ph2 = -mrkp + reshape(r'*MRKPxt', 3, numel(iDepth)) + Mdt;
end

%ph1 = -krkp + reshape(KRKPx'*r, 3, numel(iDepth)) + Kdt;
%ph2 = -mrkp + reshape(MRKPx'*r, 3, numel(iDepth)) + Mdt;

%ph1 = -krkp + reshape(r'*KRKPxt_, 3, numel(iDepth)) + Kdt;
% r set to 0 always thus
%%%%%%%%%%%%
ph1 = -krkp;% + Kdt; % this defines the weight per pixel (3rd line of it) !!
% same sign ? -- see below:
%%%%%%%%%%%%

%or |-q1*krkp3+krkp1 + q1*(KRKPx*r)3 + q1*(K*t)3*iD - (KRKPx*r)1 -(K*t)1*iD)|
% krkp(3,:) is 1+Kpx*r = K*R*p
rhs = bsxfun(@times, q1(1:2,:), ph1(3,:)) + krkp(1:2,:); % -q1*krkp3+krkp1
%rhs = -bsxfun(@times, q1(1:2,:), krkp(3,:)) + krkp(1:2,:); % -q1*krkp3+krkp1
A_r1 = bsxfun(@times, q1(1,:), KRKPxt_(:,3:3:end) ) - KRKPxt_(:,1:3:end); % q1*(KRKPx*r)3 - (KRKPx*r)1
A_r2 = bsxfun(@times, q1(2,:), KRKPxt_(:,3:3:end) ) - KRKPxt_(:,2:3:end);

A_r  = reshape( cat(1, A_r1, A_r2), 3, size(A_r1,2)*2 );

A_t1 = bsxfun(@times, K(3:3:end)', iDepth.* q1(1,:) ) - bsxfun(@times, K(1:3:end)', iDepth ); % q1*(K*t)3*iD - (K*t)1*iD
A_t2 = bsxfun(@times, K(3:3:end)', iDepth.* q1(2,:) ) - bsxfun(@times, K(2:3:end)', iDepth );

A_t  = reshape( cat(1, A_t1, A_t2), 3, size(A_t1,2)*2 );
A = cat(1, A_r, A_t)';

if exist('l1_error','var') && numel(l1_error) >1
 
%  fullL1 = l1_error(1:2:end)  +   l1_error(2:2:end);
  fullL1 = max( l1_error(1:2:end), l1_error(2:2:end));
 
  % half are outlier:
  outlier = ( fullL1 > median(fullL1) );outlier = cat( 1, outlier, outlier );
  diagOut = 1./max(0.1, sqrt(l1_error));
  diagOut(outlier) = 0.001;
  % or (L2)
  diagOut = ones(numel(diagOut),1);
%  diagOut = 1./cat(1, squeeze(-ph1(3,:)), squeeze(-ph1(3,:)) );diagOut=diagOut(:);%.^2
  diagOut(outlier) = 0;

  B = bsxfun(@times, A, diagOut);

%  A   = cat(1,A,   mat );  
  B   = cat(1,B,   mat );
  rhs = cat(1, diagOut.*rhs(:), zeros(6,1) );
  rt = B \ (rhs(:));
  A=B;% well now its the same -- influences sigma though but it did before as well  
 
  if norm(rt(1:3)) > 0.2
    C=B(:,4:6);
    rt = C \ rhs(:);rt = cat(1, zeros(3,1), rt);
  end
  
else

  % new also wo l1 reweighting 
  diagOut = 1./cat(1, squeeze(-ph1(3,:)), squeeze(-ph1(3,:)) );diagOut=diagOut(:);%.^2;
  diagOut = ones(size(diagOut));

  A = bsxfun(@times, A, diagOut);
  B   = cat(1,A,   mat );
  rhs = cat(1, diagOut.*rhs(:), zeros(6,1) );
  rt = B \ ((rhs(:)));
  A=B;
  if norm(rt(1:3)) > 0.2
    %fix rt to 0:
    C=A(:,4:6);
    rt = (C'*C) \ C'*rhs(:);rt = cat(1, zeros(3,1), rt);
  end
end

% only eval error
if evalOnly
  rt = RT;
end

l1_error = abs(A*rt-rhs(:));
l1_error = l1_error(1:end-6);

Rt = -cat(1, rt(4:6), rt(1:3));
R_cross = [0, -Rt(6), Rt(5); Rt(6), 0, -Rt(4); -Rt(5), Rt(4), 0];

r = Rt(4:6);
sinA = norm(r);
alpha = asin(sinA);
cosA = cos(alpha);
%center = zeros( 3, 1 );

R1 = zeros(3);
R1(1,1:3:end) = cosA;
R1(2,2:3:end) = cosA;
R1(3,3:3:end) = cosA;

R2 = R_cross;
R3=(1-cosA) * (r*r')./max( eps, sinA^2);

Rot = R1 + R2 + R3;
Rt_ = [R*Rot,Rt(1:3)];

if ~isreal(Rt_)
  breakhere = 1; % bad fitr check
end
