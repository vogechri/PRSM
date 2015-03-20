% computes the terms lambda_uv and u as defined in the reference 
%
function [u, lambda_uv] = Compute_lambdaUV(cam, N, M)

%ru = ones(M,N,3);
u  = ones(M,N,3);
u(:,:,1) = repmat( [1:N],  M, 1 );
u(:,:,2) = repmat( [1:M]', 1, N );

if isstruct(cam)
  viewDir = cam.Kl(1:3,3);
  R_l     = cam.Kl*cam.Rl;
else
  R_l = cam;
  [K Rc] = rq(R_l);
  Rc = Rc ./ K(3,3);
  K = K ./ K(3,3);
  viewDir = K(1:3,3);
end
% WOW GL encodes the image coordiantes by -255/256 .. 255/256
% when having 256 pixels in each dimension or -(N-1/N) .. N-1/N
% and -(M-1/M) .. M-1/M with a stepsize of 2/N of course (N pixel).
% This leads to weird recalculations everywhere!

%u(:,:,1) = repmat( [-(N-1)/N:2/N:(N-1)/N],  M, 1 );
%u(:,:,2) = repmat( [-(M-1)/M:2/M:(M-1)/M]', 1, N );

%u(:,:,1) = repmat( [-255/256:1/128:255/256],  M, 1 );
%u(:,:,2) = repmat( [-255/256:1/128:255/256]', 1, N );

% test
%{
Viewport = [ (N-1)/2 0 (N+1)/2 ; 0 (M-1)/2 (M+1)/2 ; 0  0  1];
v = inv(Viewport);
test = zeros(3,1);
for i=1:N
  for j=1:M
    test(1)  = u(i,j,1);
    test(3)  = u(i,j,3);
    test(2)  = u(i,j,2);
    u(i,j,:) = v*test;
  end
end
%}
iR = inv(R_l);

% this applies only if we are using normalized (GL coordinates)

% SCHIEF
%scale = iR * [0 0 1]';

% GERADE
%scale = iR * [(N+1)/2 (M+1)/2 1]';
scale = iR * viewDir;
% richtig: ?

scale = scale / norm(scale);
scale = iR' * scale;

lambda_uv = scale(1) * u(:,:,1) + scale(2) * u(:,:,2) + scale(3) * u(:,:,3);
lambda_uv = 1 ./ lambda_uv;
%u(:,:,3) = - u(:,:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if we are dealing with the car people, Andreas Wedel, etc ground truth
%lambda_uv = ones(M,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% das geht besser wenn man iR in scale und rotation zerlegt
