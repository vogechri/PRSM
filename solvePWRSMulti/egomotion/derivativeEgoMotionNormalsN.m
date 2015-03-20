function [F, dF] = derivativeEgoMotionNormalsN(N, kt, p2dSAll, p2dAAll, p2dGAll)

iDepth  = ( sum( bsxfun(@times, p2dSAll, N), 1 ));
q2dRAll = p2dAAll-bsxfun(@times, kt,iDepth);

iNen = 1./q2dRAll(3,:,:);
q2dR = bsxfun(@times, q2dRAll, iNen);

rx = (p2dGAll(1,:,:)-q2dR(1,:,:));
ry = (p2dGAll(2,:,:)-q2dR(2,:,:));

F = sum( squeeze(rx.^2 + ry.^2), 2)';

Z1 = bsxfun( @times, kt(1,:), p2dAAll(3,:,:) ) - bsxfun( @times, kt(3,:), p2dAAll(1,:,:) );
Z2 = bsxfun( @times, kt(2,:), p2dAAll(3,:,:) ) - bsxfun( @times, kt(3,:), p2dAAll(2,:,:) );

dFx = sum( bsxfun(@times,  2*rx .*iNen.^2 .* Z1, p2dSAll), 3 );
dFy = sum( bsxfun(@times,  2*ry .*iNen.^2 .* Z2, p2dSAll), 3 );

dF = dFx + dFy;
