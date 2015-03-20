function flow_write (F,filename)
% saves flow field F to png file
% for details see readme.txt

F = double(F);

I(:,:,1) = uint16(max(min(shiftdim(F(:,:,1))*64+2^15,2^16-1),0));
I(:,:,2) = uint16(max(min(shiftdim(F(:,:,2))*64+2^15,2^16-1),0));
I(:,:,3) = uint16(max(min(shiftdim(F(:,:,3)),1),0));
%imwrite(I,filename);
pngStuff(filename,I);

