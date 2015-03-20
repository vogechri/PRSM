function plotDataS( Seg, img1, img2, dm1X, dm2X, fNr, par, fname )

f=figure(fNr+1000);set(f, 'visible','off');
imshow(dm1X, [0, Seg.PatchSize^2 ], 'Border','tight'), colormap(jet),axis image,axis off;
%tmp1 = getframe(f);tmp1=tmp1.cdata;
tmp1 = export_fig(f, '-nocrop', '-a1');
close(f);

f=figure(fNr+1000);set(f, 'visible','off');
imshow(dm2X, [0, Seg.PatchSize^2 ], 'Border','tight'), colormap(jet),axis image,axis off;
%tmp2 = getframe(f);tmp2=tmp2.cdata;
tmp2 = export_fig(f, '-nocrop', '-a1');
close(f);

Img1 = repmat(uint8(img1*255), [1,1,3]);
Img2 = repmat(uint8(img2*255), [1,1,3]);

f=figure(fNr);set(f, 'visible','off');
II = cat ( 1, tmp1*0.5 + Img1 * 0.5, tmp2*0.5 + Img2 * 0.5);
imshow(II, 'Border','tight'),axis image, axis off;
%subplot(2,1,1);
%imshow(tmp1*0.5 + Img1 * 0.5, 'Border','tight'),axis image;
%subplot(2,1,2); imagesc(dm2X), colorbar;
%imshow(tmp2*0.5 + Img2 * 0.5, 'Border','tight'),axis image;
truesize(fNr);
export_fig( sprintf('%s/DataScores_%d_%s.png', par.sFolder, par.imgNr, fname), '-m1');
close(fNr);
