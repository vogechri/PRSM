%   TGV variational optical flow (stereo)
%
%   Author: Christoph Vogel

function [u, v, w, p, q] = tgv_run_flow(I1, I2, u, v, w, p, q, par )

alpha0 = 5;  % 5 ?? to be evaluated
alpha1 = 1;  % 1 

if ~par.doStereo
  alpha0 = 1.5; % yes, or 1.5
end

[M, N, C] = size(I1);

doMexFile =1;

if ~doMexFile
if par.doStereo
  g = edge_Tensor(I1);
  G = System_Matrix( alpha1, alpha0, M,N, par, g);% add edge based weighting later
else
  G = System_Matrix( alpha1, alpha0, M,N, par );% add edge based weighting later
end
[u,v,w,p,q] = solve_TGV( I1, I2, G, u,v, w, p,q, par);

else 
  if par.doStereo
    g = edge_Tensor(I1);
    [yids, xids, vids] = Matrix_data( alpha1, alpha0, M,N, par, g );
    [u,v,w,p,q] = fullTGV_Step( I1,I2, par.ring, u,v, par.lambda, yids, xids, vids, max(yids), max(xids), cat(1, u(:), v(:), w(:) ), cat(1, p(:), q(:)), par.warps, par.maxits, par.cEps, 1, par.dataTerm  ); 
  else
    [yids, xids, vids] = Matrix_data( alpha1, alpha0, M,N, par );
    [u,v,w,p,q] = fullTGV_Step( I1,I2, par.ring, u,v, par.lambda, yids, xids, vids, max(yids), max(xids), cat(1, u(:), v(:), w(:) ), cat(1, p(:), q(:)), par.warps, par.maxits, par.cEps, 0, par.dataTerm  );
  end
end

if par.doStereo
  % do that only late if ever - results get worse indeed
  u = min(u, 0); % for some reason (smoothing) can push to > 0
end
