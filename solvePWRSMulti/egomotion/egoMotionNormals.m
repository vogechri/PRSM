%%% assume egomotion, fix motion to ego-motion, but adjust scale of normals 
%%% to reproduce the former 2d posiiton in the other image.
%%% this is the case if the normal is wrong but the 2d motion is correct
%%% then the adjustment affects the normal - not the motion part
function [N_scale] = egoMotionNormals(ref, cam, Seg, N_lin, Rt_lin, N_res, Rt_res)

tLocOld = reshape( Rt_res(1:3,4, :), 3, size(Rt_res,3));
tLoc    = reshape( Rt_lin(1:3,4, :), 3, size(Rt_lin,3));
kt      = cam.Kl * tLoc;
centers = findPlaneCenter( Seg, 0, N_lin );
dStart  = centers(3,:);

centers2D = cam.Kl * centers;
centers2D = bsxfun( @rdivide, centers2D,  centers2D(3,:));

centers = findPlaneCenter( Seg, 0, N_res );
centers2D_goal = cam.Kl * (centers+tLocOld);
centers2D_goal = bsxfun( @rdivide, centers2D_goal,  centers2D_goal(3,:));

% x contains the scales of the normals
[x, fX_] = minimizeVektor(1./dStart, @derivativeEgoMotionNormals, [25 0.01], kt, centers2D, centers2D_goal);
N_scale = bsxfun(@times, N_lin(1:3,:), x.*dStart );

% probe: ok
%centers = findPlaneCenter( Seg, 0, N_scale );
%dStart  = centers(3,:);
%[x, fX_] = minimizeVektor(1./dStart, @derivativeEgoMotionNormals, [25 0.01], kt, centers2D, centers2D_goal);