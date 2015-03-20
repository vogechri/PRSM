function plotDepthMapsVC (dS, Is, str, par, myMinMaxD, myMaxMinD)

if ~exist('myMinMaxD','var')
  myMinMaxD = 0;
end
if ~exist('myMaxMinD','var')
  myMaxMinD = -2;
end

  temp = 1./dS(:,:,1);
  mind = min(temp(:));maxd = min(0, max(temp(:)));
for i=2:size(dS,3)
  mind = min(mind, min(temp(:)));maxd = max( maxd, max(temp(:)));
end
%temp = 1./dL1;
%mind = min(mind, min(temp(:)));maxd = max( maxd, max(temp(:)));
%temp = 1./dR0;
%mind = min(mind, min(temp(:)));maxd = max( maxd, max(temp(:)));
%mind = min(temp(:));maxd = min(0, max(temp(:)));
maxd = min(myMinMaxD, maxd); mind = max(myMaxMinD, mind);

for i=1:size(dS,3)
  temp = 1./dS(:,:,i);  
  temp = max( temp, mind);
  temp = min( temp, maxd);
  f=figure(7);set(f, 'visible','off');sc(temp, 'hicontrast');
  %f=figure(7); imshow(temp,[mind,maxd], 'Border','tight'), colormap(hsv),axis image,axis off;
  tmp2 = export_fig(f, '-nocrop', '-a1');
  %tmp2 = getframe(f);tmp2=tmp2.cdata;
  close(f);

  i1 = repmat(uint8(Is(:,:,i)*255), [1,1,3]);
  f=figure(7);set(f, 'visible','off');
  I1 = tmp2*0.5 + i1 * 0.5;
  imshow(I1, 'Border','tight'),axis image, axis off,set(f, 'visible','off');
  truesize(7);
  export_fig( sprintf('%s/%03d_%s_cam%d.png', par.sFolder, par.imgNr, str,i ), '-m1');
  close(7);
end
