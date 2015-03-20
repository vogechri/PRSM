function [cam, ref, imageName, flow2DGt, flow2DGt_noc] = loadChallangingFlow(strFolder, nr, p)

% loading the gt or some other stuff - no need for that only for painting

ImgL      = cell(1,2*p.frames+2);
ImgR      = cell(1,2*p.frames+2);

strNumber0 = sprintf( '%06d', nr );
strNumber1 = sprintf( '%06d',nr+1 );
imageName = sprintf('image%s', strNumber0);

%if strcmp(strFolder , '../../data/SceneFlow-data/RoundAbout' );
  iName = '';
%else
%  iName = 'image';
%end

ImgL{1} = double(sub_imread(sprintf('%s/%s%s_0.pgm', strFolder, iName, strNumber0 )))/4096;
ImgR{1} = double(sub_imread(sprintf('%s/%s%s_1.pgm', strFolder, iName, strNumber0 )))/4096;
ImgL{2} = double(sub_imread(sprintf('%s/%s%s_0.pgm', strFolder, iName, strNumber1 )))/4096;
ImgR{2} = double(sub_imread(sprintf('%s/%s%s_1.pgm', strFolder, iName, strNumber1 )))/4096;

if nr > 0
  strNumber2 = sprintf( '%06d',nr-1 );
  ImgLold = double(sub_imread(sprintf('%s/%s%s_0.pgm', strFolder, iName, strNumber2 )))/4096;
  ImgRold = double(sub_imread(sprintf('%s/%s%s_1.pgm', strFolder, iName, strNumber2 )))/4096;
end

strNumber3 = sprintf( '%06d',nr+2 );
if exist( sprintf('%s/%s%s_0.pgm', strFolder, iName, strNumber3 ) , 'file')
  ImgLnew = double(sub_imread(sprintf('%s/%s%s_0.pgm', strFolder, iName, strNumber3 )))/4096;
  ImgRnew = double(sub_imread(sprintf('%s/%s%s_1.pgm', strFolder, iName, strNumber3 )))/4096;
end


%data_supp = loadCalibrationKITTI_new( calibFile );
%data_supp.I = ImgR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sx = 1517.29;
sy = 1517.29*1.00746;
x0 = 392.266;
y0 = 263.092;
b  = 0.388933;
minX = 0;
minY = 0;
maxY = 542;
maxX = 657;

ImgR{1}  = ImgR{1} ( minY+1:maxY-1, minX+1:maxX-1 );
ImgR{2}  = ImgR{2} ( minY+1:maxY-1, minX+1:maxX-1 );
ImgL{1}  = ImgL{1} ( minY+1:maxY-1, minX+1:maxX-1 );
ImgL{2}  = ImgL{2} ( minY+1:maxY-1, minX+1:maxX-1 );

if nr > 0
ImgLold = ImgLold  ( minY+1:maxY-1, minX+1:maxX-1 );
ImgRold = ImgRold  ( minY+1:maxY-1, minX+1:maxX-1 );
end

% for texture
% imwrite(ImgL{1}, sprintf('%s/%s/%s%s_c0_unwarped.png', strFolder, 'gray-left',  iName, strNumber0 ));
% imwrite(ImgR{1}, sprintf('%s/%s/%s%s_c1_unwarped.png', strFolder, 'gray-left',  iName, strNumber0 ));
% imwrite(ImgL{2}, sprintf('%s/%s/%s%s_c0_unwarped.png', strFolder, 'gray-left',  iName, strNumber1 ));
% imwrite(ImgR{2}, sprintf('%s/%s/%s%s_c1_unwarped.png', strFolder, 'gray-left',  iName, strNumber1 ));

x0 = x0 - minX;
y0 = y0 - minY;
[M, N] = size (ImgL{1});

Viewport = [ (N)/2 0 (N+1)/2 ; 0 (M)/2 (M+1)/2 ; 0  0  1 ];
FlipY = [ 1 0 0 ; 0 -1 M+1; 0 0 1];
%FlipY = [ 1 0 0 ; 0 1 0; 0 0 1];

K   = [sx 0 x0; 0 sy y0; 0 0 1];
R_l = [1 0 0 0; 0 1 0 0; 0 0 1 0];
R_r = [1 0 0 -b; 0 1 0 0; 0 0 1 0];

% rectified image
T = [-b 0 0]';
Rot = eye(3);

Pl = FlipY * K * R_l;
Pr = FlipY * K * R_r;


% for compatibility to openGL version
%Pl = inv (Viewport) * Pl;
%Pr = inv (Viewport) * Pr;

% drawImages:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A, B, C] = svd (Pl); % ! C is nor C' : F = A*B*C' ! -> epi2 is col of C
cl = C(:, end);
el = Pr * cl;

F = Pr * ( Pl'*inv(Pl*Pl') );

ex = [ [0 , -el(3), el(2)]; [ el(3), 0, -el(1)]; [-el(2) , el(1), 0]];
F = ex * F;

F = F / norm(F,'fro');

epi = el;
Rr = Pr(:,1:3);
popl = Pl(:,end);
popr = Pr(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_supp(1).I = ImgR;
data_supp(1).R = Rr;
data_supp(1).F = F;
data_supp(1).Ft = zeros(3);
data_supp(1).epi = epi;
data_supp(1).popRef = popl;
data_supp(1).pop = popr;
data_supp(1).Rot = Rot;
data_supp(1).Tra = T;

[data_supp(1).Kl, data_supp(1).Rl, data_supp(1).Tl, pp] = cameraParameters ( Pl );
[data_supp(1).Kr, data_supp(1).Rr, data_supp(1).Tr, pp] = cameraParameters ( Pr );
data_supp(1).Tl = - inv(data_supp(1).Rl) * data_supp(1).Tl;
data_supp(1).Tr = - inv(data_supp(1).Rr) * data_supp(1).Tr;
data_supp(1).Rl = inv(data_supp(1).Rl);
data_supp(1).Rr = inv(data_supp(1).Rr);

flow2DGt     = repmat( zeros(size(ImgR{1})), [1,1,4] );
flow2DGt_noc = repmat( zeros(size(ImgR{1})), [1,1,4] );
[cam, ref] = generateStructures ( data_supp, ImgL, data_supp.R );

if nr > 0
  ref.Iold = ImgLold;
  cam.Iold = ImgRold;
end

if exist( 'ImgLnew' , 'var')
  ref.Inew = ImgLnew;
  cam.Inew = ImgRnew;
end