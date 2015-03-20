function [camc, refc, imageName, flow2DGt, flow2DGt_noc] = loadZealandFlow( strFolder, nr, p, iName )

% loading the gt or some other stuff - no need for that only for painting

ImgL      = cell(1,2*p.frames+2);
ImgR      = cell(1,2*p.frames+2);

strNumber0 = sprintf( '%04d', nr );
strNumber1 = sprintf( '%04d',nr+1 );
imageName = sprintf('image%s', strNumber0);

%if strcmp(strFolder , '../../data/SceneFlow-data/RoundAbout' );
%  iName = '';
%else
%  iName = 'image';
%end

% 12 bit, so 
%temp1 = double(sub_imread(sprintf('%s/%s%s_0.pgm', strFolder, iName, strNumber0 )))/4096;

ImgL{1} = double(sub_imread(sprintf('%s/left_%s/img_%s_c0.pgm', strFolder, iName, strNumber0 )))/65535;
ImgR{1} = double(sub_imread(sprintf('%s/right_%s/img_%s_c2.pgm', strFolder, iName, strNumber0 )))/65535;
ImgL{2} = double(sub_imread(sprintf('%s/left_%s/img_%s_c0.pgm', strFolder, iName, strNumber1 )))/65535;
ImgR{2} = double(sub_imread(sprintf('%s/right_%s/img_%s_c2.pgm', strFolder, iName, strNumber1 )))/65535;
ImgC{1} = double(sub_imread(sprintf('%s/center_%s/img_%s_c1.pgm', strFolder, iName, strNumber0 )))/65535;
ImgC{2} = double(sub_imread(sprintf('%s/center_%s/img_%s_c1.pgm', strFolder, iName, strNumber1 )))/65535;

if nr > 0
  strNumber2 = sprintf( '%04d',nr-1 );
  ImgLold = double(sub_imread(sprintf('%s/left_%s/img_%s_c0.pgm', strFolder, iName, strNumber2 )))/65535;
  ImgRold = double(sub_imread(sprintf('%s/right_%s/img_%s_c2.pgm', strFolder, iName, strNumber2 )))/65535;
  ImgCold = double(sub_imread(sprintf('%s/center_%s/img_%s_c1.pgm', strFolder, iName, strNumber2 )))/65535;
end

strNumber3 = sprintf( '%04d',nr+2 );
if exist( sprintf('%s/left_%s/img_%s_c0.pgm', strFolder, iName, strNumber3), 'file')
  ImgLnew = double(sub_imread(sprintf('%s/left_%s/img_%s_c0.pgm', strFolder, iName, strNumber3 )))/65535;
  ImgRnew = double(sub_imread(sprintf('%s/right_%s/img_%s_c2.pgm', strFolder, iName, strNumber3 )))/65535;
  ImgCnew = double(sub_imread(sprintf('%s/center_%s/img_%s_c1.pgm', strFolder, iName, strNumber3 )))/65535;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sx       = 1031.02;
sy       = 0.998045*1031.02;
x0       = 302.454;
y0       = 285.46;

b1  = 0.80537 ; % left - right
b2  = 0.505707; % left - center

minX = 0;
minY = 0;
maxY = 480;
maxX = 640;
%{
ImgR{1}  = ImgR{1} ( minY+1:maxY-1, minX+1:maxX-1 );
ImgR{2}  = ImgR{2} ( minY+1:maxY-1, minX+1:maxX-1 );
ImgL{1}  = ImgL{1} ( minY+1:maxY-1, minX+1:maxX-1 );
ImgL{2}  = ImgL{2} ( minY+1:maxY-1, minX+1:maxX-1 );

if nr > 0
  ImgLold = ImgLold  ( minY+1:maxY-1, minX+1:maxX-1 );
  ImgRold = ImgRold  ( minY+1:maxY-1, minX+1:maxX-1 );
end

if exist( 'ImgLnew' , 'var')
  ImgLnew = ImgLnew ( minY+1:maxY-1, minX+1:maxX-1 );
  ImgRnew = ImgRnew ( minY+1:maxY-1, minX+1:maxX-1 );
end
%}

x0 = x0 - minX;
y0 = y0 - minY;
[M, N] = size (ImgL{1});

%Viewport = [ (N)/2 0 (N+1)/2 ; 0 (M)/2 (M+1)/2 ; 0  0  1 ];
FlipY = [ 1 0 0 ; 0 -1 M+1; 0 0 1];
%FlipY = [ 1 0 0 ; 0 1 0; 0 0 1];

K   = [sx 0 x0; 0 sy y0; 0 0 1];
R_l = [1 0 0 0; 0 1 0 0; 0 0 1 0];
R_r = [1 0 0 -b1; 0 1 0 0; 0 0 1 0];
R_c = [1 0 0 -b2; 0 1 0 0; 0 0 1 0];

% rectified image
T1 = [-b1 0 0]';
T2 = [-b2 0 0]';
Rot = eye(3);

Pl = FlipY * K * R_l;
Pr = FlipY * K * R_r;
Pc = FlipY * K * R_c;

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
Rc = Pc(:,1:3);
popc = Pc(:,end);


ec = Pc * cl;
Fc = Pc * ( Pl'*inv(Pl*Pl') );
exc = [ [0 , -ec(3), ec(2)]; [ ec(3), 0, -ec(1)]; [-ec(2) , ec(1), 0]];
Fc = exc * Fc;
Fc = Fc / norm(Fc,'fro');
epic = ec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_supp(1).I = ImgR;
data_supp(1).R = Rr;
data_supp(1).F = F;
data_supp(1).Ft = zeros(3);
data_supp(1).epi = epi;
data_supp(1).popRef = popl;
data_supp(1).pop = popr;
data_supp(1).Rot = Rot;
data_supp(1).Tra = T1;

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
%  ref.Iold = ImgLold;
  cam.Iold = {ImgLold, ImgRold};
end

if exist( 'ImgLnew' , 'var')
%  ref.Inew = ImgLnew;
  cam.Inew = {ImgLnew, ImgRnew};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_supp(1).I   = ImgC;
data_supp(1).R   = Rc;
data_supp(1).F   = Fc;
data_supp(1).Ft  = zeros(3);
data_supp(1).epi = epic;
data_supp(1).popRef = popl;
data_supp(1).pop = popc;
data_supp(1).Rot = Rot;
data_supp(1).Tra = T2;

[data_supp(1).Kl, data_supp(1).Rl, data_supp(1).Tl, pp] = cameraParameters ( Pl );
[data_supp(1).Kr, data_supp(1).Rr, data_supp(1).Tr, pp] = cameraParameters ( Pc );
data_supp(1).Tl = - inv(data_supp(1).Rl) * data_supp(1).Tl;
data_supp(1).Tr = - inv(data_supp(1).Rr) * data_supp(1).Tr;
data_supp(1).Rl = inv(data_supp(1).Rl);
data_supp(1).Rr = inv(data_supp(1).Rr);

flow2DGt     = repmat( zeros(size(ImgC{1})), [1,1,4] );
flow2DGt_noc = repmat( zeros(size(ImgC{1})), [1,1,4] );
[camc, refc] = generateStructures ( data_supp, ImgL, data_supp.R );

if nr > 0
%  ref.Iold = ImgLold;
  camc.Iold = ImgCold;
end

if exist( 'ImgCnew' , 'var')
%  ref.Inew = ImgLnew;
  camc.Inew = ImgCnew;
end

if nr > 1
%  ref.Iold = ImgLold;
  camc.Iold = {ImgLold, ImgCold};
end

if exist( 'ImgCnew' , 'var')
%  ref.Inew = ImgLnew;
  camc.Inew = {ImgLnew, ImgCnew};
end


%cam(2) = camc;
