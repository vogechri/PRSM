function center = findPlaneCenter(Seg, p2d_, N_)

center = cat( 2, Seg.Centers, ones(size(Seg.Centers, 1), 1) )';
center = -bsxfun( @rdivide, center, sum( center .* N_(1:3,:), 1) );