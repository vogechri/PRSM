function err = disp_error(D_gt,D_est,tau)

E = abs(D_gt-D_est);
E(D_gt<=0) = 0;
err = length(find(E>tau))/length(find(D_gt>0));
