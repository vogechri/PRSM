%   TGV variational optical flow (stereo)
%
%   Author: Christoph Vogel

function [u, v, w, p, q] = tgv_level_flowDouble(I1, I2, I3, I4, u, v, w, p, q, par )

alpha0 = 5;  % 5 ?? to be evaluated tripple check for kitti 
alpha1 = 1;  % 1 

if ~par.doStereo
  alpha0 = 1.5; % 1 or 1.5 ? 1: non-occ, 1.5 occ better
end

[M, N, ~] = size(I1);
if par.doStereo
  g = edge_Tensor(I1);

%  if alpha0 == 0 %par.doTV
%    [yids, xids, vids] = Matrix_data_TV( alpha1, M,N, par );%, level_constraints );
%  else  
    [yids, xids, vids] = Matrix_data( alpha1, alpha0, M,N, par, g );
%  end
  
  
% p1=I1;p2=I2;p3=par.ring;p4=u;p5=v;p6=par.lambda;p7=yids;p8=xids;
% p9=vids;p10=max(yids);p11=max(xids);p12=cat(1, u(1:N*M)', v(1:N*M)', w(1:N*M*4)' );
% p13=cat(1, p(1:N*M*4)', q(1:N*M*8)');p14=par.warps;p15=par.maxits;p16=par.cEps;
% p17=1;p18=par.dataTerm;p19=0;p20=I3;p21=I4;
% p22=squeeze(u(:,:,2));p23=squeeze(v(:,:,2));
% p24=cat(1, u(N*M+1:2*N*M)', v(N*M+1:2*N*M)', w(1+N*M*4:end)' );
% p25=cat(1, p(1+N*M*4:end)', q(1+N*M*8:end)');
%    save( sprintf( 'of_censPicDouble.mat' ),'p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11','p12','p13','p14','p15','p16','p17','p18', 'p19','p20','p21','p22','p23','p24','p25');
  
  [u,v,w,p,q] = ofmexD( I1,I2, par.ring, squeeze(u(:,:,1)),squeeze(v(:,:,1)), ...
                        par.lambda, yids, xids, vids, max(yids), max(xids), ...
                        cat(1, u(1:N*M), v(1:N*M), w(1:N*M*4) ), cat(1, p(1:N*M*4), q(1:N*M*8)), ...
                        par.warps, par.maxits, par.cEps, 1, par.dataTerm, 0, ...
                        I3, I4, squeeze(u(:,:,2)),squeeze(v(:,:,2)), ...
                        cat(1, u(N*M+1:2*N*M), v(N*M+1:2*N*M), w(1+N*M*4:end) ), ...
                        cat(1, p(1+N*M*4:end), q(1+N*M*8:end)) );

else
  [yids, xids, vids] = Matrix_data( alpha1, alpha0, M,N, par );

%  [u1,v1,w1,p1,q1, u2,v2,w2,p2,q2] 
  [u,v,w,p,q] = ofmexD( I1,I2, par.ring, squeeze(u(:,:,1)),squeeze(v(:,:,1)), ...
                        par.lambda, yids, xids, vids, max(yids), max(xids), ...
                        cat(1, u(1:N*M)', v(1:N*M)', w(1:N*M*4)' ), cat(1, p(1:N*M*4)', q(1:N*M*8)'), ...
                        par.warps, par.maxits, par.cEps, 0, par.dataTerm, 0, ...
                        I3, I4, squeeze(u(:,:,2)),squeeze(v(:,:,2)), ...
                        cat(1, u(N*M+1:2*N*M)', v(N*M+1:2*N*M)', w(1+N*M*4:end)' ), ...
                        cat(1, p(1+N*M*4:end)', q(1+N*M*8:end)') );  
end

% u = cat(3, u1,u2);
% v = cat(3, v1,v2);
% w = cat(3, w1,w2);
% p = cat(3, p1,p2);
% q = cat(3, q1,q2);

if par.doStereo
  % do that only late if ever - results get worse indeed
  u = min(u, 0); % for some reason (smoothing) can push to > 0
end
