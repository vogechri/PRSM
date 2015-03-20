%%% A helpful plotting function
function plotAnalysis(ref, cam, N_res, Rt_res, Seg, u, fNr, par, use_centered_view,...
                      ending, occIds1, occIds2, occIds3, noOccl)

% nice hack to plot flow only: use header below and disable computation of
% flow_ (line 43)
%function plotAnalysis(ref, cam, flow_, Seg, u, fNr, par, ending, occIds1, occIds2, occIds3, noOccl)

if ~exist('use_centered_view', 'var')
  use_centered_view = 0;
end
% plot wx,wy,wz
plot3DFlow = 0;

if ~exist('fNr', 'var')
  fNr = 1;
end

if ~exist('ending', 'var')
  ending = 'pA2';
end

[m,n,~] = size(u);

maxMot = 1;
maxD   = 30;%15;

writeImg = 1;
plotdata = 1;
plotmatlab = 0;

if writeImg
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plot nice figures:
  I_key = ref.I(1).I;
end

p2d_ = permute( reshape ( inv(cam(1).Kl) * reshape(permute (u , [3,1,2]), 3, n*m), [3,m,n]), [2,3,1]);
if use_centered_view
  center = findPlaneCenter( Seg, p2d_, N_res);
else
  center = zeros(size(N_res(1:3,:)));  
end

flow_ =  reconstruc3DFlowHom( N_res, Rt_res, Seg, p2d_, use_centered_view );
origD = squeeze(flow_(:,:,1));

% to look at in matlab
if plotmatlab
%  maxMot = 0.125;
  figure(fNr  ), imagesc(flow_(:,:,1), [0,maxD] ), colorbar
  figure(fNr+1), imagesc(flow_(:,:,2), [-maxMot,maxMot]), colorbar
  figure(fNr+2), imagesc(flow_(:,:,3), [-maxMot,maxMot]), colorbar
  figure(fNr+3), imagesc(flow_(:,:,4), [-maxMot,maxMot]), colorbar
  return;
end

if plotdata
  cam = Homography_cam_pixel_Center(ref, cam, 1, N_res, Rt_res, Seg, 1, center);
  cam = Homography_cam_pixel_Center(ref, cam, 1, N_res, Rt_res, Seg, 2, center);
  ref = Homography_ref_pixel_Center(ref, cam, 2, N_res, Rt_res, Seg, center);
  
  [ ref, cam ] = computeWarpedImagesCam( ref, cam );
  
  f=figure(fNr+4);set(f, 'visible','off'),imshow(ref.I(1).I, 'Border','tight'), colormap(gray);
  if writeImg
    export_fig( sprintf('%s/originalL.png', par.sFolder), '-m1');
    close force all
  end
  f=figure(fNr+5);set(f, 'visible','off'),imshow(cam(1).I(1).I_w, 'Border','tight'), colormap(gray);
  if writeImg
    export_fig( sprintf('%s/%s_rt_HQ.png', par.sFolder, ending), '-m1');
    close force all
  end
  if writeImg
    f = figure(fNr+6);set(f, 'visible','off'),imshow(cam.I(1).I, 'Border','tight'), colormap(gray);
    export_fig( sprintf('%s/originalR.png', par.sFolder), '-m1');
    close(f);
    close force all
    f = figure(fNr+6);set(f, 'visible','off'),imshow(ref.I(2).I, 'Border','tight'), colormap(gray);
    export_fig( sprintf('%s/originalLt.png', par.sFolder), '-m1');
    close(f);
    close force all
    f = figure(fNr+6);set(f, 'visible','off'),imshow(cam.I(2).I, 'Border','tight'), colormap(gray);
    export_fig( sprintf('%s/originalRt.png', par.sFolder), '-m1');
    close(f);
    close force all
  end
  f=figure(fNr+7);set(f, 'visible','off'),imshow(ref(1).I(2).I_w, 'Border','tight'), colormap(gray);
  if writeImg
    export_fig( sprintf('%s/%s_lt1_HQ.png', par.sFolder, ending), '-m1');
    close force all
  end
  f=figure(fNr+8);set(f, 'visible','off'),imshow(cam(1).I(2).I_w, 'Border','tight'), colormap(gray);
  if writeImg
    export_fig( sprintf('%s/%s_rt1_HQ.png', par.sFolder, ending), '-m1');
    close force all
  end
end

if writeImg
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if size(I_key,3) > 1
    I_key1 = mean(I_key,3);
  else
    I_key1 = I_key;
  end
  I_key1 = uint8(I_key1*255);
  I_key3 = repmat(I_key1, [1,1,3]);

  tmp = flow_(:,:,1);
  tmp(tmp> maxD) = maxD;
  f = figure(fNr+17);set(f, 'visible','off');
  imshow(tmp, [min(tmp(:)),max(tmp(:))], 'Border','tight'), colormap(jet),axis image,axis off;
  tmp = export_fig(f, '-nocrop', '-a1');
  
  tmp = 0.5*tmp+0.5*I_key3;
  close(f);f = figure(fNr+17);set(f, 'visible','off');
  imshow(tmp), axis image, axis off;
  export_fig( sprintf('%s/%s_d_HQ.png', par.sFolder, ending), '-m1');


  tmp = flow_(:,:,1);tmp( tmp > maxD/2) = maxD/2;
  close(f);f = figure(fNr+17);set(f, 'visible','off');
  imshow(tmp, [min(tmp(:)),max(tmp(:))], 'Border','tight'), colormap(jet),axis image,axis off;
  tmp = export_fig(f, '-nocrop', '-a1');
  
  tmp = 0.5*tmp+0.5*I_key3;
  close(f);f = figure(fNr+17);set(f, 'visible','off');
  imshow(tmp), axis image, axis off;
  export_fig( sprintf('%s/%s_d_HQ_d2.png', par.sFolder, ending), '-m1');

  tmp = flow_(:,:,1);tmp( tmp > maxD/3) = maxD/3;
  close(f);f = figure(fNr+17);set(f, 'visible','off');
  imshow(tmp, [min(tmp(:)),max(tmp(:))], 'Border','tight'), colormap(jet),axis image,axis off;
  tmp = export_fig(f, '-nocrop', '-a1');
  
  tmp = 0.5*tmp+0.5*I_key3;
  close(f);f = figure(fNr+17);set(f, 'visible','off');
  imshow(tmp), axis image, axis off;
  export_fig( sprintf('%s/%s_d_HQ_d3.png', par.sFolder, ending), '-m1');


  tmp = origD;tmp( tmp > maxD*1.5) = maxD*1.5;
  close(f);f = figure(fNr+17);set(f, 'visible','off');
  imshow(tmp, [min(tmp(:)),max(tmp(:))], 'Border','tight'), colormap(jet),axis image,axis off;
  tmp = export_fig(f, '-nocrop', '-a1');
  
  tmp = 0.5*tmp+0.5*I_key3;
  close(f);f = figure(fNr+17);set(f, 'visible','off');
  imshow(tmp), axis image, axis off;
  export_fig( sprintf('%s/%s_d_HQ_d4.png', par.sFolder, ending), '-m1');

  tmp = origD;tmp( tmp > maxD*2.5) = maxD*2.5;
  close(f);f = figure(fNr+17);set(f, 'visible','off');
  imshow(tmp, [min(tmp(:)),max(tmp(:))], 'Border','tight'), colormap(jet),axis image,axis off;
  tmp = export_fig(f, '-nocrop', '-a1');
  
  tmp = 0.5*tmp+0.5*I_key3;
  close(f);f = figure(fNr+17);set(f, 'visible','off');
  imshow(tmp), axis image, axis off;
  export_fig( sprintf('%s/%s_d_HQ_d5.png', par.sFolder, ending), '-m1');

  %-%%%%%%%%%%%%%%%%%%%%%%%%%%
  [disp, uDisp, vDisp, dDisp] = convert3Dto2D(ref, cam, flow_(:,:,1), flow_(:,:,2), flow_(:,:,3), flow_(:,:,4));

  %-% disparity image
  tmp = -disp;
  f = figure(fNr+18);set(f, 'visible','off');
  sc(tmp, 'hicontrast');
  tmp = export_fig(f, '-nocrop', '-a1');  
  tmp = 0.5*tmp+0.5*I_key3;
  close(f);f = figure(fNr+18);set(f, 'visible','off');
  imshow(tmp), axis image, axis off;
  export_fig( sprintf('%s/%s_disp_HQ.png', par.sFolder, ending), '-m1');

  tmp = dDisp;
  f = figure(fNr+19);set(f, 'visible','off');
  sc(tmp, 'hicontrast');
  tmp = export_fig(f, '-nocrop', '-a1');
  tmp = 0.5*tmp+0.5*I_key3;
  close(f);f = figure(fNr+18);set(f, 'visible','off');
  imshow(tmp), axis image, axis off;
  export_fig( sprintf('%s/%s_dispDiff_HQ.png', par.sFolder, ending), '-m1');


  if exist('occIds1', 'var' )
    % occ
    dd = origD;
    occludedArea = false(size(Seg.Img));
    if numel(occIds1) == numel(Seg.Img)
      occludedArea( occIds1 ~=0 ) = true;
    else
      for i=1:numel(occIds1)
        occludedArea (Seg.Img == occIds1(i)) = true;
      end
    end
    tmp = dd;tmp( tmp > maxD) = maxD;
    f = figure(fNr+18);set(f, 'visible','off');
    imshow(tmp, [min(tmp(:)),max(tmp(:))], 'Border','tight'), colormap(jet),axis image,axis off;
    tmp = export_fig(f, '-nocrop', '-a1');
    tmp( repmat( occludedArea, [1,1,3]) ) = 255;% white, black also ok
    
    tmp = 0.5*tmp+0.5*I_key3;
    close(f);f = figure(fNr+18);set(f, 'visible','off');
    imshow(tmp), axis image, axis off;
    export_fig( sprintf('%s/%s_d_HQ_Occ.png', par.sFolder, ending), '-m1');
  end

    max_flow = max([abs(uDisp(:)); abs(vDisp(:))]);
    F_mag = sqrt(uDisp.*uDisp+vDisp.*vDisp);
    F_dir = atan2(vDisp,uDisp);
    I_flow2d = flow_map(F_mag(:),F_dir(:),ones(1,numel(uDisp)),max_flow,8);
    I_flow2d = reshape(I_flow2d, [size(uDisp, 1), size(uDisp, 2), 3]);


    if exist('occIds2', 'var' )
      occludedArea = false(size(Seg.Img));
      if numel(occIds2) == numel(Seg.Img)
        occludedArea( occIds2 ~=0 ) = true;
      else
        for i=1:numel(occIds2)
          occludedArea (Seg.Img == occIds2(i)) = true;
        end
      end
    end
    f = figure(fNr+18);set(f, 'visible','off');
    imshow(I_flow2d, 'Border','tight'),axis image,axis off;
    tmp = export_fig(f, '-nocrop', '-a1');
  
    if exist('occIds2', 'var' )
      tmp( repmat( occludedArea, [1,1,3]) ) = 0;% white, black also ok
    end
    
    tmp = 0.5*tmp+0.5*I_key3;
    close(f);f = figure(fNr+18);set(f, 'visible','off');
    imshow(tmp), axis image, axis off;
    export_fig( sprintf('%s/%s_d_HQ_OccLT.png', par.sFolder, ending), '-m1');
  
  if exist('occIds3', 'var' )
    % occ
    tmp = cat( 3, flow_(:,:,2), flow_(:,:,3), flow_(:,:,4));tmp=tmp./max(tmp(:));
    occludedArea = false(size(Seg.Img));
    if numel(occIds3) == numel(Seg.Img)
      occludedArea( occIds3 ~=0 ) = true;
    else
      for i=1:numel(occIds3)
        occludedArea (Seg.Img == occIds3(i)) = true;
      end
    end
    f = figure(fNr+18);set(f, 'visible','off');
    imshow(tmp, 'Border','tight'),axis image,axis off;
    tmp = export_fig(f, '-nocrop', '-a1');
    tmp( repmat( occludedArea, [1,1,3]) ) = 255;% white, black also ok
    
    tmp = 0.5*tmp+0.5*I_key3;
    close(f);f = figure(fNr+18);set(f, 'visible','off');
    imshow(tmp), axis image, axis off;
    export_fig( sprintf('%s/%s_d_HQ_OccRT.png', par.sFolder, ending), '-m1');
  end

if plot3DFlow
  close(f);f = figure(fNr+18);set(f, 'visible','off');
  tmp = flow_(:,:,2);
  if (max(tmp(:)) < -maxMot)
    maxMot = abs ( max(tmp(:)) );
  end
  if (min(tmp(:)) >  maxMot)
    maxMot = abs ( min(tmp(:)) );
  end
  imshow(tmp, [max(-maxMot, min(tmp(:))),min(maxMot, max(tmp(:)))], 'Border','tight'), colormap(jet),axis image,axis off;
  tmp = export_fig(f, '-nocrop', '-a1');  
  tmp = 0.5*tmp+0.5*I_key3;
  close(f);f = figure(fNr+17);set(f, 'visible','off');
  imshow(tmp), axis off;
  export_fig( sprintf('%s/%s_wX_HQ.png', par.sFolder, ending), '-m1');
  
  close(f);f = figure(fNr+19);set(f, 'visible','off');
  tmp = flow_(:,:,3);
  if (max(tmp(:)) < -maxMot)
    maxMot = abs ( max(tmp(:)) );
  end
  if (min(tmp(:)) >  maxMot)
    maxMot = abs ( min(tmp(:)) );
  end
  imshow(tmp, [max(-maxMot, min(tmp(:))),min(maxMot, max(tmp(:)))], 'Border','tight'), colormap(jet),axis image,axis off;
  tmp = export_fig(f, '-nocrop', '-a1');  
  tmp = 0.5*tmp+0.5*I_key3;
  close(f);f = figure(fNr+17);set(f, 'visible','off');
  imshow(tmp), axis off;
  export_fig( sprintf('%s/%s_wY_HQ.png', par.sFolder, ending), '-m1');
  
  close(f);f = figure(fNr+20);set(f, 'visible','off');
  tmp = flow_(:,:,4);
  if (max(tmp(:)) < -maxMot)
    maxMot = abs ( max(tmp(:)) );
  end
  if (min(tmp(:)) >  maxMot)
    maxMot = abs ( min(tmp(:)) );
  end
  imshow(tmp, [max(-maxMot, min(tmp(:))),eps+min(maxMot, max(tmp(:)))], 'Border','tight'), colormap(jet),axis image,axis off;
  tmp = export_fig(f, '-nocrop', '-a1');  
  tmp = 0.5*tmp+0.5*I_key3;
  close(f);f = figure(fNr+17);set(f, 'visible','off');
  imshow(tmp), axis off;
  export_fig( sprintf('%s/%s_wZ_HQ.png', par.sFolder, ending), '-m1');
  close(f);
  % }
end

  if exist('noOccl', 'var' )
    noOccImg = zeros(size(Seg.Img)); 
    for i=1:numel(Seg.Ids) 
      noOccImg(Seg.Img==i-1)=noOccl(i);
    end
    f = figure(fNr+18);set(f, 'visible','off');
    imshow(log(noOccImg)./max(log(noOccImg(:))), 'Border','tight'), colormap(jet),axis image,axis off;
    tmp = export_fig(f, '-nocrop', '-a1');
    tmp = 0.5*tmp+0.5*I_key3;
    close(f);f = figure(fNr+18);set(f, 'visible','off');
    imshow(tmp), axis image, axis off;
    export_fig( sprintf('%s/%s_d_HQ_noSol.png', par.sFolder, ending), '-m1');
  end

  pause(1);% less UIJ_AreThereWindowShowsPending errors ? since then windows really are closed
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = flow_map(F_mag,F_dir,F_val,max_flow,n)

I(:,1) = mod(F_dir/(2*pi),1);
I(:,2) = F_mag * n / max_flow;
I(:,3) = n - I(:,2);
I(:,[2 3]) = min(max(I(:,[2 3]),0),1);
I = hsv2rgb(reshape(I,[],1,3));

I(:,:,1) = I(:,:,1); %.* F_val;
I(:,:,2) = I(:,:,2); %.* F_val;
I(:,:,3) = I(:,:,3); %.* F_val;
