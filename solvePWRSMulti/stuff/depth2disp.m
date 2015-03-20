function disparity = depth2disp(u, d, epi, R_l, R_r, lambda_uv, q2dx)

[M N K] = size(u);

qu = zeros(M,N,K);
zz = zeros(M,N,K);
urt = zeros(M,N,2);

% R_r = inv(ViewP) * K * Rot
% R_l = inv(ViewP) * K
% -> inv(R_l)*R_r = Rot
Q = R_r * inv (R_l);

qu(:,:,1) = lambda_uv .*(Q(1,1) * u(:,:,1) + Q(1,2) * u(:,:,2) + Q(1,3) * u(:,:,3));
%qu(:,:,2) = lambda_uv .*(Q(2,1) * u(:,:,1) + Q(2,2) * u(:,:,2) + Q(2,3) * u(:,:,3));
qu(:,:,3) = lambda_uv .*(Q(3,1) * u(:,:,1) + Q(3,2) * u(:,:,2) + Q(3,3) * u(:,:,3));

zz(:,:,1) = epi(1) + d .* qu(:,:,1);
%zz(:,:,2) = epi(2) + d .* qu(:,:,2);
zz(:,:,3) = epi(3) + d .* qu(:,:,3);

urt (:,:,1) = zz(:,:,1) ./ zz(:,:,3);
%urt (:,:,2) = zz(:,:,2) ./ zz(:,:,3);

disp2D = urt(:,:,1) - u(:,:,1);
disparity = disp2D ./ q2dx;
% uses sin a = ankathete / hypothenuse, a linear relationship, (sin a != 0)

% disp2D = ur - u(:,:,1:2);
% F = disp2D(:,:,1).^2 + disp2D(:,:,2).^2;
% F = sqrt(F);
% 
% dF = disp2D(:,:,1) .* dur(:,:,1) + disp2D(:,:,2) .* dur(:,:,2);
% dF = dF ./F;



% Ansatz in 2D: e + disp * q = ur
% disp2D(:,:,1) = ur(:,:,1) - e(1)'/e(3);
% disp2D(:,:,2) = ur(:,:,2) - e(2)'/e(3);
% F = disp2D(:,:,1).^2 + disp2D(:,:,2).^2;
% F = sqrt(F);
% 
% dF = disp2D(:,:,1) .* dur(:,:,1) + disp2D(:,:,2) .* dur(:,:,2);
% dF = dF ./F;



% [N M K] = size(ur);
% one = ones(N,M);
% 
% % Ansatz in 2D: u + disp * q = ur
% disp2D = ur - u(:,:,1:2);
% 
% disp = disp2D(:,:,1).^2 + disp2D(:,:,2).^2;
% disp = sqrt(disp);
% 
% % nix anderes als (disp2D / |disp2D|) * dur = |dur|
% dF = disp2D(:,:,1) .* dur_d(:,:,1) + disp2D(:,:,2) .* dur_d(:,:,2);
% one(dF < 0) = -1;
% disp = disp .* one;
% 
% %dF = dF ./F;

