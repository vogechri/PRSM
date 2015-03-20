function [F, dF] = derivativeEgoMotionNormals(d, kt, centers2D, centers2D_goal)

nen = (d.*kt(3,:)+centers2D(3,:));

rx = centers2D_goal(1,:) - (centers2D(1,:)+d.*kt(1,:))./ nen;
ry = centers2D_goal(2,:) - (centers2D(2,:)+d.*kt(2,:))./ nen;

F = rx.^2 + ry.^2;

dFx = 2*rx .* ( ( d.*kt(1,:).*centers2D(3,:)  - d.*kt(3,:).* centers2D(1,:) ) ./ nen.^2 );
dFy = 2*ry .* ( ( d.*kt(2,:).*centers2D(3,:)  - d.*kt(3,:).* centers2D(2,:) ) ./ nen.^2 );

dF = -dFx - dFy;

%F = F';
%dF = dF';