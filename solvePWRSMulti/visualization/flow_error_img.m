function E = flow_error_img(F_gt,F_est)

F_gt_du  = shiftdim(F_gt(:,:,1));
F_gt_dv  = shiftdim(F_gt(:,:,2));
F_gt_val = shiftdim(F_gt(:,:,3));

F_est_du = shiftdim(F_est(:,:,1));
F_est_dv = shiftdim(F_est(:,:,2));

E_du = F_gt_du-F_est_du;
E_dv = F_gt_dv-F_est_dv;
E    = sqrt(E_du.*E_du+E_dv.*E_dv);
E(F_gt_val==0) = 0;