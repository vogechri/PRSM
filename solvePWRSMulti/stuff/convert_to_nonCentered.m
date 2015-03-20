%%% helper function convert rigid motion per segment center to origin
%%% centered motion
function [N_lin_prop, Rt_lin_prop, centers2D_prop] = convert_to_nonCentered( Seg, N_prop, RT_prop )

centers2D = cat(2, Seg.Centers, ones (size(Seg.Centers,1),1))';
centers2D_prop = [];
N_lin_prop = [];
Rt_lin_prop = [];
k=0;
while size(N_lin_prop,2) < size(N_prop,2)
  Rt_lin = RT_prop(:,:,k*numel(Seg.Ids)+1:(k+1)*numel(Seg.Ids) );
  N_lin  = N_prop(1:3, k*numel(Seg.Ids)+1:(k+1)*numel(Seg.Ids) );
  centers = findPlaneCenter( Seg, 0, N_lin );
  for i = 1:size(centers, 2)
    Rt_lin(1:3,4,i) =  Rt_lin(1:3,4,i) - Rt_lin(1:3,1:3,i) * centers(:,i) + centers(:,i);
  end
  % rest:
  N_lin_prop  = cat(2, N_lin_prop, N_lin);
  Rt_lin_prop = cat(3, Rt_lin_prop, Rt_lin);
  centers2D_prop = cat(2, centers2D_prop, centers2D);
  k=k+1;
end
