function [data_cam] = loadCalibrationKITTI_new( filename )


% get the camera intrinsic and extrinsic calibration
%calib = loadCalibrationCamToCam(fullfile(dir,'calib_cam_to_cam.txt'));

%    calib.S_rect{cam} = S_rect_;
%    calib.R_rect{cam} = R_rect_;

% open file
fid = fopen(filename,'r');

if fid<0
  calib = [];
  return;
end
res = fscanf(fid, 'P0: %f %f %f %f %f %f %f %f %f %f %f %f \nP1: %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

calib.P_rect{1} = reshape(res(1:12), [4,3])';
calib.P_rect{2} = reshape(res(13:24), [4,3])';

T = inv(calib.P_rect{1}(1:3,1:3)) * calib.P_rect{2}(:,4);
K = calib.P_rect{1}(1:3,1:3);

data_cam.Kl = calib.P_rect{1}(1:3,1:3);
data_cam.Kr = calib.P_rect{1}(1:3,1:3);
data_cam.Rl = eye(3);
data_cam.Rr = eye(3);
data_cam.Tl = inv(calib.P_rect{1}(1:3,1:3)) * calib.P_rect{1}(:,4);
data_cam.Tr = inv(calib.P_rect{2}(1:3,1:3)) * calib.P_rect{2}(:,4);

data_cam.Rot = eye(3);
data_cam.Tra = data_cam.Tr;

data_cam.pop    = zeros(3,1);
data_cam.popRef = (data_cam.Kr)*data_cam.Tr;
data_cam.epi    = (data_cam.Kr)*data_cam.Tr;

data_cam.F  = [0,0,0;0,0,0.707106781186548;0,-0.707106781186548,0];
data_cam.Ft = zeros(3);
data_cam.R  = data_cam.Kl;

%N = calib.S_rect{1}(1);M = calib.S_rect{1}(2);
%Viewport = [ (N)/2 0 (N+1)/2 ; 0 (M)/2 (M+1)/2 ; 0  0  1 ];
%%FlipY = [ 1 0 0 ; 0 -1 M+1; 0 0 1];
%FlipY = [ 1 0 0 ; 0 1 0; 0 0 1];


%[data_supp(1).Kl, data_supp(1).Rl, data_supp(1).Tl, pp] = cameraParameters ( Pl );
%[data_supp(1).Kr, data_supp(1).Rr, data_supp(1).Tr, pp] = cameraParameters ( Pr );

