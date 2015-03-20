function [depth dD_dDisp dD2_dDisp2] = disp2depth(disp, u, e, R_l, R_r, lambda_uv, q2dx) % dur_d
%epsilon = 0.0001;

% inverse:
%function [ weight ] = Compute_ur_weight(R_l, R_r, epi, u, lambda_uv, d, dd )
[M N K] = size(u);
qu = zeros(M,N,K);
Q = R_r * inv (R_l);

qu(:,:,1) = lambda_uv .*(Q(1,1) * u(:,:,1) + Q(1,2) * u(:,:,2) + Q(1,3) * u(:,:,3));
%qu(:,:,2) = lambda_uv .*(Q(2,1) * u(:,:,1) + Q(2,2) * u(:,:,2) + Q(2,3) * u(:,:,3));
qu(:,:,3) = lambda_uv .*(Q(3,1) * u(:,:,1) + Q(3,2) * u(:,:,2) + Q(3,3) * u(:,:,3));

%q2dx = dur(:,:,1) ./ sqrt(dur(:,:,1).^2 + dur(:,:,2).^2);
temp = u(:,:,1) + disp .* q2dx;
DZ   =       e(1) -      e(3)  * temp;
DN   = -qu(:,:,1) + qu(:,:,3) .* temp;

depth = DZ ./ DN;

%test = dCheck - depth; % yep

dD_dDispZ = - e(3) * q2dx;
dD_dDispN = qu(:,:,3) .* q2dx;

dD_dDisp = -DZ.*dD_dDispN + DN .*dD_dDispZ;


%dD2_dDisp2 = dD_dDispN * (dD_dDisp); % Hmmm why not .* check with
%dD2_dDisp2 = dD2_dDisp2 ./ (DN.^3);
%checkDeriv

DN2 = DN.^2;

dD_dDisp = dD_dDisp ./ DN2;

dD2_dDisp2 = -2*dD_dDispN ./ DN  .* dD_dDisp;
%dD2_dDisp2 = dD2_dDisp2 ./ (DN2);
end


% disp2D = ur - u(:,:,1:2);
% F = disp2D(:,:,1).^2 + disp2D(:,:,2).^2;
% F = sqrt(F);
% 
% dF = disp2D(:,:,1) .* dur(:,:,1) + disp2D(:,:,2) .* dur(:,:,2);
% dF = dF ./F;


% % Ansatz in 2D: e + disp * q = ur
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
% F = disp2D(:,:,1).^2 + disp2D(:,:,2).^2;
% F = sqrt(F);
% 
% % nix anderes als (disp2D / |disp2D|) * dur = |dur|
% dF = disp2D(:,:,1) .* dur(:,:,1) + disp2D(:,:,2) .* dur(:,:,2);
% one(dF < 0) = -1;
% F = F .* one;
% 
% dF = dF ./F;

