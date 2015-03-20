function [ref, cam ] = computeWarpedImagesCam( ref, cam )

%Bicubic_MexWindows
if size (ref.I(1).I, 3) < 2
  for j=1:size(cam,2)
    for k=1:numel(cam(j).I)
      [cam(j).I(k).I_w, cam(j).I(k).I_wx, cam(j).I(k).I_wy] = ...
        Interpol_mex( cam(j).I(k).I, cam(j).I(k).u(:,:,1), cam(j).I(k).u(:,:,2));
    end
  end
else
  for j=1:size(cam,2)
    for k=1:numel(cam(j).I)
      for l = 1:size (ref.I(1).I, 3)
        [cam(j).I(k).I_w(:,:,l), cam(j).I(k).I_wx(:,:,l), cam(j).I(k).I_wy(:,:,l)] = ...
          Interpol_mex(cam(j).I(k).I(:,:,l), cam(j).I(k).u(:,:,1), cam(j).I(k).u(:,:,2));
      end
    end
  end
end

if size (ref.I(1).I, 3) < 2
  for k=2:numel(ref.I)
    [ref.I(k).I_w, ref.I(k).I_wx, ref.I(k).I_wy] = ...
      Interpol_mex(ref.I(k).I, ref.I(k).u(:,:,1), ref.I(k).u(:,:,2));
  end
else
  for k=2:numel(ref.I)
    for l = 1:size (ref.I(1).I, 3)
      [ref.I(k).I_w(:,:,l), ref.I(k).I_wx(:,:,l), ref.I(k).I_wy(:,:,l)] = ...
        Interpol_mex(ref.I(k).I(:,:,l), ref.I(k).u(:,:,1), ref.I(k).u(:,:,2));
    end
  end
end
