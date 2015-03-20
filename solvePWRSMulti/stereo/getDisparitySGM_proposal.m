%%% get opencv stereo estiamte
function DispImg__x1 = getDisparitySGM_proposal(cam, ref, stereoT_2d_new, minDisp, maxDisp, windowsize, sits )

if ~exist('sits','var')
  sits = 10;
end

[M,N,~] = size(cam.I(1).I);
u  = ones(M,N,3);
u(:,:,1) = repmat( [1:N],  M, 1 );
u(:,:,2) = repmat( [1:M]', 1, N );

DispImg__ = getDisparitySGM(cam, ref, minDisp, maxDisp, windowsize, 1 );
oobs = (u(:,:,1)-stereoT_2d_new(:,:,1) < 0);
DispImg__(oobs) = -1;

DispImg__fill = DispImg__;
DispImg__fill(DispImg__<0) = stereoT_2d_new(DispImg__<0);

DispImg__fix = DispImg__;
DispImg__fix(1:end,1) = -1;
DispImg__x = TvL1Matrix(sits, size(DispImg__,1), size(DispImg__,2), int32(DispImg__fix<0), DispImg__fill, DispImg__fill);

% smooth once to kill staircase ?
for i=1:3
  DispImg__x = medfilt2(DispImg__x,[3 3],'symmetric');
end
% this part looses a lot of infiormation:
DispImg__x1 = TvL1Matrix(1, size(DispImg__x,1), size(DispImg__x,2), int32(ones(size(DispImg__fix))), DispImg__x, DispImg__x);
