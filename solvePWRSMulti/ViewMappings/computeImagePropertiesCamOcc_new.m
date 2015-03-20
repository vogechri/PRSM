function [ref, cam] = computeImagePropertiesCamOcc_new(ref, u, ...
  lambda_uv, cam, x, doOccMap)

R_l = ref.R;

% for each image one map (not for the reference cam) 1 is t1 of reference
%occlusions = struct( 'map', {} );

addBorder = 0.5;

[M, N] = size (x.d_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% right at t0

for j=1:size( cam, 2 )

    [cam(j).I(1).u, cam(j).I(1).du_d, oMap] = Compute_ur_all_d_weight(R_l, ...
      cam(j).R, cam(j).epi, u, lambda_uv, x.d_, doOccMap);
  
  if doOccMap
    cam(j).I(1).oMap = oMap;
  else
    if ~isempty(cam(j).I(1).oMap)
      oMap = cam(j).I(1).oMap;
    end
  end

  % boundary handling % promote to epi-constraints as well ???
    m = (cam(j).I(1).u(:,:,1) > N+addBorder) | (cam(j).I(1).u(:,:,1) < 1-addBorder) | ...
        (cam(j).I(1).u(:,:,2) > M+addBorder) | (cam(j).I(1).u(:,:,2) < 1-addBorder) | (oMap(:,:) == 0);

  cam(j).I(1).valid = m;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% right at t+1

for k=1:numel(x.w_)
  ik = k+1;
  w_ = x.w_(k);
  for kk= k-2:-2:1
    w_.x = w_.x + x.w_(kk).x;
    w_.y = w_.y + x.w_(kk).y;
    w_.z = w_.z + x.w_(kk).z;
  end
  
  for j=1:size( cam, 2 )

      [cam(j).I(ik).u, cam(j).I(ik).du_d, cam(j).I(ik).du_wx, cam(j).I(ik).du_wy, cam(j).I(ik).du_wz, ...
        oMap] = Compute_urx_all_d_weight(R_l, cam(j).R, u,...
        cam(j).epi, lambda_uv, x.d_, w_.x, w_.y, w_.z, doOccMap );
    
    if doOccMap
      cam(j).I(ik).oMap = oMap;
    else
      if ~isempty(cam(j).I(ik).oMap)
        oMap = cam(j).I(ik).oMap;
      end
    end

    % boundary handling % promote to epi-constraints as well ???
      m = (cam(j).I(ik).u(:,:,1) > N+addBorder) | (cam(j).I(ik).u(:,:,1) < 1-addBorder) | ...
          (cam(j).I(ik).u(:,:,2) > M+addBorder) | (cam(j).I(ik).u(:,:,2) < 1-addBorder) | (oMap(:,:) == 0);
    cam(j).I(ik).valid = m;%repmat(m, [1,1,5]);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% left at t+1
for k=1:numel(x.w_)
  ik = k+1;
  w_ = x.w_(k);
  for kk= k-2:-2:1
    w_.x = w_.x + x.w_(kk).x;
    w_.y = w_.y + x.w_(kk).y;
    w_.z = w_.z + x.w_(kk).z;
  end
  
  [ref.I(ik).u, ref.I(ik).du_d, ref.I(ik).du_wx, ref.I(ik).du_wy, ref.I(ik).du_wz, oMap] = ...
    Compute_ulx_all_d_weight(R_l, u, lambda_uv, x.d_, w_.x, w_.y, w_.z, doOccMap);
  % boundary handling % promote to epi-constraints as well ???
  
  if doOccMap
    ref.I(ik).oMap = oMap;
  else
    if ~isempty(ref.I(ik).oMap)
      oMap = cam(j).I(1).oMap;
    end
  end
  
    m = (ref.I(ik).u(:,:,1) > N+addBorder) | (ref.I(ik).u(:,:,1) < 1-addBorder) | ...
        (ref.I(ik).u(:,:,2) > M+addBorder) | (ref.I(ik).u(:,:,2) < 1-addBorder) | (oMap(:,:) == 0);
  ref.I(ik).valid = m;
  
end
