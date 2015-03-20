function DispImg_=getDisparitySGM(cam, ref, minDisp, maxDisp, windowsize, timestep, par )
  
[M, N ] = size(cam.I(1).I);
Mnew = M;
Nnew = N;

inftyDisp = 0.05;

if ~exist('timestep','var')
  timestep = 1;
end

img1 = ref.I(timestep).I;
img2 = cam.I(timestep).I;
% obsolete ?! 
%{
if ~exist('par','var')
  par.use_structure_texture    = 1;
  par.structure_texture_sigma  = 1;
  par.structure_texture_factor = 0.8;
end

img1=structure_texture_decomposition(img1, par);
img2=structure_texture_decomposition(img2, par);
%}

useRect = 0;
img1 =  uint8(255*img1);
img2 =  uint8(255*img2);

Kr  = cam.Kr(1,3);
Kl  = cam.Kl(1,3);
Dl = [0 0 0 0];
Dr = [0 0 0 0];
R = eye(3);
T = cam.Tr;

i1 = [zeros( Mnew, maxDisp) img1 ];
i2 = [zeros( Mnew, maxDisp-1) img2 zeros( Mnew, 1)];
Kl(1,3) = cam.Kl(1,3) + maxDisp;
Kr(1,3) = cam.Kr(1,3) + maxDisp;

% original p1=windowsize^2*channels*4 == 36 if windowsize==3 must be dividable
% by 4 p2 = 8*p1
% before kitti version - 
%[DispImg__,~]  = Stereo_SGM( Mnew, Nnew+maxDisp, i1, i2, Kl, Kr, Dl, Dr, R, T, minDisp+1, maxDisp, windowsize, useRect, 16, 8*36);

[DispImg__,~]  = Stereo_SGM( Mnew, Nnew+maxDisp, i1, i2, Kl, Kr, Dl, Dr, R, T, minDisp+1, maxDisp, windowsize, useRect, 50, 800);

% old stuff
%[DispImg__ Q H1 H2] = Stereo_SGM( Mnew, Nnew+maxDisp, i1, i2, Kl, Kr, Dl, Dr, R, T, minDisp+1, maxDisp, windowsize, useRect);

% from kitti web: - better  - well seems so 
%[DispImg__,~]  = Stereo_SGM_kitti( Mnew, Nnew+maxDisp, i1, i2, Kl, Kr, Dl, Dr, R, T, minDisp+1, maxDisp, windowsize, useRect, 50, 800);

DispImg_ = DispImg__(:, maxDisp+1:end) - (minDisp+1);

% unreliable:
%DispImg_(:,1:20) = -1;
occMask = zeros(size(DispImg_));
occMask(DispImg_ < 0) = 1;
DispImg_(DispImg_ == 0) = min( inftyDisp, min(DispImg_(DispImg_ > 0)) / 2);
