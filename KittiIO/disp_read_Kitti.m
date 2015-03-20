function D = disp_read_Kitti (filename)
% loads disparity map D from png file
% for details see readme.txt

%I =  imread(filename);
I =  pngStuff(filename);
D = -double(I)/256;
