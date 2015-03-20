%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compute the internal and external camera parameters from the projection
%%% m
%%%%%%%
%%% P: projection matrix 
%%%
%%% K: calibration matrix 
%%% R: camera orientation in World coordinates
%%% C: camera center
%%% p: the principal point in image coordinates (pixel)
function [K, R, C, p] = cameraParameters ( P )

%% TO FILL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% apply rq decomposition on the first 3x3 matrix of P
[K, RC] = rq( P(1:3, 1:3) );

% Camera system is right handed by definition (see lecture slides)! 
% the rq-decomposition does not guarantee to deliver a rotation matrix but
% only a orthonormal matrix
if (det (RC) < 0)
  RC = -RC;
  K = -K;
end

%% extract C, the camera center (in world coordinates)
[U,D,V] = svd (P);
CC = V(:,end);
C = CC / CC(end); % set homogeneous coord to 1

% P is invariant to scale and we want inner calibration in pixel
scale = K(3,3);
K = K / scale;

% we want the values to be positive on the diagonal,
% since focal distance and pixel size are positive numbers
% NO MATLAB has to be negative entry on the 2nd diagonal entry
D = eye(3);
for idx = 1:size(K, 1)
  if mod(idx,2)
    if K(idx, idx) < 0
      D(idx, idx) = -1;
    end
  else
    if K(idx, idx) > 0
      D(idx, idx) = -1;
    end    
  end
end

K  = K * D;
RC = D * RC;

% the principal point
p = [ K(1,3) K(2,3)];

%camera orientation wanted here
R = inv(RC);

% Camera system is right handed by definition (see lecture slides)! 
% the rq-decomposition does not guarantee to deliver a rotation matrix but
% only a orthonormal matrix
if (det (R) < 0)
  R = -R;
%  K = -K;
end

% camera center
C = C(1:3);

% this matrix is identical to P except scale, but due to the conversion to
% cartesian coordinates, the scale does not harm the projection
%testP = K*[inv(R) -inv(R)*C(1:3)];
