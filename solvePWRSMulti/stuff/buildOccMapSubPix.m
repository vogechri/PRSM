% build the occlusion map consiting of 1: free and 0 : occluded
% input pixC  : the coordinates in pixel (also )
% input depth : the depth at the pixel, at the position pixC
%
% the whole process is problematic since here either 1.0 -> 3-6 and
% 1.75-eps -> 5-8, leading to an overlap of 2, while 1.25 and 1.75 also
% have that overlap, but 1.25-eps and 1.5 also only have an overlap of 2
% so requiring an overlap of 2 or 3 for occlusion leads to bad results
% either way
function oMap = buildOccMapSubPix ( pixC, depth, SegImg )

global occlusionDepth;

if numel(SegImg) > 1 % since i prevent segments to occlude themselves ok
  storeOD = occlusionDepth;
  occlusionDepth = 0.00001;
end
%idx = uint16 (round( pixC (:,:,1) ));
%idy = uint16 (round( pixC (:,:,2) ));

idxMin  = int16((pixC(:,:,1) - 0.375) * 4)-3; % values from 0 to .. 4N-1
idxMin2 = int16((pixC(:,:,1) - 0.25 ) * 4)-3;
idyMin  = int16((pixC(:,:,2) - 0.375) * 4)-3; % values from 0 to .. 4M-1
idyMin2 = int16((pixC(:,:,2) - 0.25 ) * 4)-3;

%idxMin = int16((pixC(:,:,1)) * 4)-4; % values from 0 to .. 4N-1
%idyMin = int16((pixC(:,:,2)) * 4)-4;% values from 0 to .. 4M-1

% test : 
%for i = 0.975:0.0125:3 fprintf('%f: %d - %d and %d - %d\n', i, int16((i - 0.375) * 4),int16((i + 0.375) * 4), int16((i - 0.25) * 4),int16((i + 0.5) * 4)); end
% inside : zBuf, oMap, 
% NO: idiotic
%oMap = OcclusionMapSubPix_MexWindows( idxMin, idyMin, idxMin2, idyMin2, depth, occlusionDepth );
if numel(SegImg) > 1
  oMap = OcclusionMapSubPixSeg_MexWindows( idyMin, idxMin, idyMin2, idxMin2, depth, occlusionDepth, SegImg );
else
  oMap = OcclusionMapSubPix_MexWindows( idyMin, idxMin, idyMin2, idxMin2, depth, occlusionDepth );
end
% debug:
%oMap = OcclusionMapSubPix_MexWindows_write( idyMin, idxMin, idyMin2, idxMin2, depth, occlusionDepth );

%oMap = OcclusionMapSubPix_MexWindows_nonStrict( idxMin, idyMin, idxMin2, idyMin2, depth, offset ); 
% min reicht: min:miin+3 : handcoded - no for
%for allPix, getDepth all 16 pixel do putValue if smallest
% test: parallel
% for allPix getDepth checkAllPixel, if sum occluded >=8 : occluded

%%%%%%%%% TEST %%%%%%%%%
%{
[N M] = size (depth);
oMap2 = ones(M,N);
  
% init map:
zBuf  = 100000 * ones(M,N);

m = (idx(:,:) <= N) & (idx(:,:) >= 1) & (idy(:,:) <= M) & (idy(:,:) >= 1);

% this has to go in a mex file

% build zBuffer
for i = 1:M
  for j = 1:N
    if m(i,j)
      id_x = idx(i,j);
      id_y = idy(i,j);
      zBuf( id_y, id_x ) = min ( zBuf( id_y, id_x ), depth(i,j) );
    end
  end
end

% build occlusion Map, with some sort of tolerance
for i = 1:M
  for j = 1:N
    if m(i,j)

      id_x = idx(i,j);
      id_y = idy(i,j);
      
      if zBuf( id_y, id_x ) + offset < depth(i,j);
        oMap2( i, j ) = 0; % occluded
      end
      
    end
  end
end

if any(any(oMap - oMap2))
  a = 2;
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test:
% oMap = OcclusionMapSubPixSeg_MexWindows( idyMin, idxMin, idyMin2, idxMin2, depth, occlusionDepth, SegImg );
% oMapTiny = OcclusionMapSubPixSeg_MexWindows( idyMin(40:80, 80:120), idxMin(40:80, 80:120), idyMin2(40:80, 80:120), idxMin2(40:80, 80:120), depth(40:80, 80:120), occlusionDepth, SegImg );

% {
se = strel('diamond',1);
%strel('square',3);
%se = strel('ball',5,5);
oMapO = imopen(oMap,se);
oMapC = imclose(oMapO,se);
oMapC = imclose(oMapC,se);

%oMapE = imerode(oMap,se);
%oMapD = imdilate(oMapE,se);
%oMapD = imdilate(oMapD,se);

oMap = oMapC;

if numel(SegImg) > 1
  occlusionDepth = storeOD;
end

% }