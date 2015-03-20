%%% find regions with similar motion close to each other -- leads to more
%%% stable proposals 
function [newIds,iSegIds] = ...
  algebraicMatrixQ9( cam, R, q1, q, p, iDepth, Idsk, Seg, validPixel, maxSegs )

if ~exist('maxSegs','var')
  maxSegs = 1000;
end

K=cam.Kl;
kp1 = cat(1, zeros(1, size( p,2 )),  p(3,:), -p(2,:) );
kp2 = cat(1, -p(3,:), zeros(1, size( p,2 )),  p(1,:) );
kp3 = cat(1,  p(2,:), -p(1,:),  zeros(1, size( p,2 )) );

KPx = zeros (3, 3*size(kp1,2));
KPx(:,1:3:end) = kp1;
KPx(:,2:3:end) = kp2;
KPx(:,3:3:end) = kp3;

%%%%%%%%%%%%%%%%%%%%%%%%%
ip   = p;
krkp = K*R*ip;
KRKPx = K*R*KPx;

% the transposed krkpx:
KRKPxt_ = cat(1, reshape(KRKPx(:, 1:3:end), 1, numel(KRKPx)/3), reshape(KRKPx(:, 2:3:end), 1, numel(KRKPx)/3, 1), reshape(KRKPx(:, 3:3:end), 1, numel(KRKPx)/3));

% r set to 0 always
%ph1 = -krkp + Kdt;

%or |-q1*krkp3+krkp1 + q1*(KRKPx*r)3 + q1*(K*t)3*iD - (KRKPx*r)1 -(K*t)1*iD)|
rhs = -bsxfun(@times, q1(1:2,:), krkp(3,:)) + krkp(1:2,:);
A_r1 = bsxfun(@times, q1(1,:), KRKPxt_(:,3:3:end) ) - KRKPxt_(:,1:3:end);
A_r2 = bsxfun(@times, q1(2,:), KRKPxt_(:,3:3:end) ) - KRKPxt_(:,2:3:end);
A_r  = reshape( cat(1, A_r1, A_r2), 3, size(A_r1,2)*2 );

% works with ids
A_t1 = bsxfun(@times, K(3:3:end)', iDepth(:)' .* q1(1,:) ) - bsxfun(@times, K(1:3:end)', iDepth(:)' );
A_t2 = bsxfun(@times, K(3:3:end)', iDepth(:)' .* q1(2,:) ) - bsxfun(@times, K(2:3:end)', iDepth(:)' );
A_t  = reshape( cat(1, A_t1, A_t2), 3, size(A_t1,2)*2 );

A = cat(2, A_r', A_t');

rhs=rhs';
add = numel(A)/12;
Aout = cell(numel(Idsk),1);
bout = cell(numel(Idsk),1);

% normals:
e  = cam.Kr*cam.Tr; % not like this??? - the same mostly
Ax = cam.Kr*cam.Rr;

Ap = Ax*p;

b1 = (bsxfun(@times, e,p(1,:)));
b2 = (bsxfun(@times, e,p(2,:)));
b3 = (bsxfun(@times, e,p(3,:)));
Bb = -cat( 1, b1(:)', b2(:)', b3(:)' );
  
bn = q-Ap;
%%%%%%%%%%%%%

MMs = zeros(9,9,numel(Idsk));
bbs = zeros(1,numel(Idsk));
Mbs = zeros(9,numel(Idsk));

for i=1:numel(Idsk)
  
  id4s = Idsk{i};
  bx  = bn(:,id4s+1);bx=bx(:);
  id3s = cat(2,1+3*id4s, 2+3*id4s, 3+3*id4s )';
  Mxx = Bb (:,id3s);

  MMxn  = Mxx*Mxx';
  Mbxn  = Mxx*bx;
  bbxn  = bx'*bx;
  %%%
  
  ids = cat(1,1+2*Idsk{i}', 2+2*Idsk{i}' );%ids = ids(:)';
  Mxx = A(ids,:);  
  bx  = cat(2, rhs(1+Idsk{i}), rhs(1+Idsk{i}+add) );%+ numel(rhs)/2 ));%, zeros(6,1) );
  bx = bx';bx=bx(:);
    
  %%% small motion prior:
  multMe = numel(Idsk{i})/256;
  ss = 0.1;% * (numel( Idsk{i})/256);
  mid = round(numel( Idsk{i} )/2);
  mid = Idsk{i}(mid);
  pd = K*p(:,mid);
  pxyz1 = pd ./ iDepth(mid);
  pdx = (pxyz1+[ss,0,0]');pdx = pdx./pdx(3);
  pdy = (pxyz1+[0,ss,0]');pdy = pdy./pdy(3);
  pdz = (pxyz1+[0,0,ss]');pdz = pdz./pdz(3);
  scalexyz = [norm(pdx-pd), norm(pdy-pd), norm(pdz-pd)];

  mat= multMe * diag( cat(2, [1,1,1], (1/ss)./scalexyz));
  Mxx  = cat(1, Mxx, mat );
  bx   = cat(1, bx, zeros(6,1) );
  %%%

  Aout{i} = Mxx;
  bout{i} = bx;

  MMx  = [MMxn,zeros(3,6);zeros(6,3), Mxx'*Mxx];
  Mbx  = cat(1, Mbxn, Mxx'*bx);
  bbx  = bx'*bx + bbxn;
  
  MMs(:,:,i) = MMx;
  Mbs(:,i)   = Mbx;
  bbs(:,i)   = bbx;
end

p1=MMs;
p2=Mbs;
p3=bbs;
p4=numel(Seg.Ids)+1;
%p5=Seg.Ids;
p5=Seg.NeighIds;
p6 = maxSegs;  % final segs

SegPix = zeros(1,numel(Seg.Ids));
for i=1:numel(Seg.Ids)
  SegPix(i) = numel(Seg.Ids{i});
end
%p7= SegPix; % # pix per seg
%p8 = 2.0;   % max factor worsen 

%save( sprintf( './testQ.mat'), 'p1','p2','p3','p4','p5', 'p6', 'p7','p8');
aa = gradGrowMatlabM9(p1,p2,p3,p4,p5, p6); % no division -- ?

%{
img = Seg.Img;
for i=1:numel(Seg.Ids)
  img(Seg.Ids{i}+1) = aa(i);
end
figure(9), imagesc(img)
%}

newIds  = cell(1,numel(unique(aa)));
iSegIds = cell(1,numel(unique(aa)));

for i=1:numel(aa)
  iSegIds{ aa(i)+1 } = cat(1, iSegIds{ aa(i)+1 }, i);
  newIds{ aa(i)+1 }  = cat(1,  newIds{ aa(i)+1}, Seg.Ids{i} );
end
