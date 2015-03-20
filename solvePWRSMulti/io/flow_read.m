function F = flow_read (filename)
% loads flow field F from png file
% for details see readme.txt

I = double(imread(filename));

F_u = (I(:,:,1)-2^15)/64;
F_v = (I(:,:,2)-2^15)/64;
F_valid = min(I(:,:,3),1);
%F_u(F_valid==0) = 0;
%F_v(F_valid==0) = 0;
F(:,:,1) = F_u;
F(:,:,2) = F_v;
F(:,:,3) = F_valid;

