% compute the position in the right image at timestep t
%
% see 3D scene flow proposal
function [ur, dur_d] = Compute_ur_all_d(R_l, R_r, epi, u, lambda_uv, d )

[M N K] = size(u);

qu = zeros(M,N,K);
zz = zeros(M,N,K);
ur = zeros(M,N,2);
%{
% test points in left view:
Q = inv (R_l);

qu(:,:,1) = lambda_uv .*(Q(1,1) * u(:,:,1) + Q(1,2) * u(:,:,2) + Q(1,3) * u(:,:,3));
qu(:,:,2) = lambda_uv .*(Q(2,1) * u(:,:,1) + Q(2,2) * u(:,:,2) + Q(2,3) * u(:,:,3));
qu(:,:,3) = lambda_uv .*(Q(3,1) * u(:,:,1) + Q(3,2) * u(:,:,2) + Q(3,3) * u(:,:,3));

zz(:,:,1) = d .* qu(:,:,1);
zz(:,:,2) = d .* qu(:,:,2);
zz(:,:,3) = d .* qu(:,:,3);
%255,100: 3d-0.211084 0.978661 -3.929999 OK till here
ur (:,:,1) = zz(:,:,1) ./ zz(:,:,3);
ur (:,:,2) = zz(:,:,2) ./ zz(:,:,3);
% end test
%}
Q = R_r * inv (R_l);

qu(:,:,1) = lambda_uv .*(Q(1,1) * u(:,:,1) + Q(1,2) * u(:,:,2) + Q(1,3) * u(:,:,3));
qu(:,:,2) = lambda_uv .*(Q(2,1) * u(:,:,1) + Q(2,2) * u(:,:,2) + Q(2,3) * u(:,:,3));
qu(:,:,3) = lambda_uv .*(Q(3,1) * u(:,:,1) + Q(3,2) * u(:,:,2) + Q(3,3) * u(:,:,3));

zz(:,:,1) = epi(1) + d .* qu(:,:,1);
zz(:,:,2) = epi(2) + d .* qu(:,:,2);
zz(:,:,3) = epi(3) + d .* qu(:,:,3);
%255 128-> 1 129: 0.015351 -2.537391 2.256235 w-comp all OK ?
ur (:,:,1) = zz(:,:,1) ./ zz(:,:,3);
ur (:,:,2) = zz(:,:,2) ./ zz(:,:,3);

nn = zz(:,:,3).^2;

dur_d (:,:,1) = qu(:,:,1) * epi(3) - qu(:,:,3) * epi(1); 
dur_d (:,:,2) = qu(:,:,2) * epi(3) - qu(:,:,3) * epi(2); 

dur_d (:,:,1) = dur_d (:,:,1) ./ nn;
dur_d (:,:,2) = dur_d (:,:,2) ./ nn;
