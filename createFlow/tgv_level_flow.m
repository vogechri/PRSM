%   TGV variational optical flow (stereo)
%
%   Author: Christoph Vogel

function [u, v, w, p, q] = tgv_level_flow(I1, I2, u, v, w, p, q, par )

alpha0 = 5;  % 5 ?? to be evaluated tripple check for kitti 
alpha1 = 1;  % 1 

if ~par.doStereo
  alpha0 = 1.5; % 1 or 1.5 ? 1: non-occ, 1.5 occ better
end

[M, N, ~] = size(I1);
if par.doStereo
  g = edge_Tensor(I1);

  % debug
%  if alpha0 == 0 %par.doTV
%    [yids, xids, vids] = Matrix_data_TV( alpha1, M,N, par );%, level_constraints );
%  else  
    [yids, xids, vids] = Matrix_data( alpha1, alpha0, M,N, par, g );
%  end
  
  
%    p1=I1;p2=I2;p3=par.ring;p4=u;p5=v;p6=par.lambda;p7=yids;p8=xids;
%    p9=vids;p10=max(yids);p11=max(xids);p12=cat(1, u(:), v(:), w(:) );
%    p13=cat(1, p(:), q(:));p14=par.warps;p15=par.maxits;p16=par.cEps;
%    p17=1;p18=par.dataTerm;
%    save( sprintf( 'of_censPic.mat' ),'p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11','p12','p13','p14','p15','p16','p17','p18');
  
%  [u,v,w,p,q] = fullTGV_Step( I1,I2, par.ring, u,v, par.lambda, yids, xids, vids, max(yids), max(xids), cat(1, u(:), v(:), w(:) ), cat(1, p(:), q(:)), par.warps, par.maxits, par.cEps, 1 );

%  [u,v,w,p,q] = fullTGV_Step_PF24( I1,I2, par.ring, u,v, par.lambda, yids, xids, vids, max(yids), max(xids), cat(1, u(:), v(:), w(:) ), cat(1, p(:), q(:)), par.warps, par.maxits, par.cEps, 1, par.dataTerm );
  [u,v,w,p,q] = fullStepTGV_PF24_s( I1,I2, par.ring, u,v, par.lambda, yids, xids, vids, max(yids), max(xids), cat(1, u(:), v(:), w(:) ), cat(1, p(:), q(:)), par.warps, par.maxits, par.cEps, 1, par.dataTerm );  

else
  [yids, xids, vids] = Matrix_data( alpha1, alpha0, M,N, par );
  
  % these 2 work differently ! for CSAD and CENSUS
%  [u,v,w,p,q] = fullTGV_Step_PF24( I1,I2, par.ring, u,v, par.lambda, yids, xids, vids, max(yids), max(xids), cat(1, u(:), v(:), w(:) ), cat(1, p(:), q(:)), par.warps, par.maxits, par.cEps, 0, par.dataTerm );
  [u,v,w,p,q] = fullStepTGV_PF24_s( I1,I2, par.ring, u,v, par.lambda, yids, xids, vids, max(yids), max(xids), cat(1, u(:), v(:), w(:) ), cat(1, p(:), q(:)), par.warps, par.maxits, par.cEps, 0, par.dataTerm );
end

if par.doStereo
  % do that only late if ever - results get worse indeed
  u = min(u, 0); % for some reason (smoothing) can push to > 0
end
