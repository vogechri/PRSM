function plotDepthMapsFromViews (dR0, dL1, dR1, I1, I2, I3, str, par)

temp = 1./dR1;
mind = min(temp(:));maxd = max(temp(:));
temp = 1./dL1;
mind = min(mind, min(temp(:)));maxd = max( maxd, max(temp(:)));
temp = 1./dR0;
mind = min(mind, min(temp(:)));maxd = max( maxd, max(temp(:)));
%mind = min(temp(:));maxd = min(0, max(temp(:)));
maxd = min(0, maxd); mind = max(-4, mind);
f=figure(7); imshow(temp,[mind,maxd], 'Border','tight'), colormap(hsv),axis image,axis off;
tmp2 = getframe(f);tmp2=tmp2.cdata;
close(f);

i1 = repmat(uint8(I1*255), [1,1,3]);
figure(7); 
I1 = tmp2*0.5 + i1 * 0.5;
imshow(I1, 'Border','tight'),axis image, axis off;
truesize(7);
export_fig( sprintf('%s/dispR0_%d_%s.png', par.sFolder, par.imgNr, str), '-m1');
close(7);
%%%%%%%%%%%%%%%
temp = 1./dL1;
%mind = min(temp(:));maxd = min(0, max(temp(:)));

f=figure(7); imshow(temp,[mind,maxd], 'Border','tight'), colormap(hsv),axis image,axis off;
tmp2 = getframe(f);tmp2=tmp2.cdata;
close(f);

i1 = repmat(uint8(I2*255), [1,1,3]);
figure(7);
I1 = tmp2*0.5 + i1 * 0.5;
imshow(I1, 'Border','tight'),axis image, axis off;
truesize(7);
export_fig( sprintf('%s/dispL1_%d_%s.png', par.sFolder, par.imgNr, str), '-m1');
close(7);
%%%%%%%%%%%%%%%%
temp = 1./dR1;
%mind = min(temp(:));maxd = min(0, max(temp(:)));

f=figure(7); imshow(temp,[mind,maxd], 'Border','tight'), colormap(hsv),axis image,axis off;
tmp2 = getframe(f);tmp2=tmp2.cdata;
close(f);

i1 = repmat(uint8(I3*255), [1,1,3]);
figure(7);
I1 = tmp2*0.5 + i1 * 0.5;
imshow(I1, 'Border','tight'),axis image, axis off;
truesize(7);
export_fig( sprintf('%s/dispR1_%d_%s.png', par.sFolder, par.imgNr, str), '-m1');
close(7);
