%%% compute first approaximation of a image gradinet based tensor
function g = edge_Tensor(f)

q=1;alpha=5;

[M,N] = size(f);
idx = bsxfun(@plus, zeros(M,N),  1:N);
idy = bsxfun(@plus, zeros(M,N), (1:M)');

[~, fx, fy] = Interpol_mex(f, idx, idy);

norm = max(1e-06, sqrt(fx.^2 + fy.^2));
small = (norm <= 1e-06);

fx = fx ./ norm;
fy = fy ./ norm;

fx(small) = 1;
fy(small) = 0;

gg = max(1e-06, exp(-alpha*norm.^q));

g =     cat(3, gg .* fx.*fx, gg .* fx.*fy, gg .* fx.*fy, gg .* fy.*fy);
g = g + cat(3,       fy.*fy,      -fx.*fy,      -fx.*fy,       fx.*fx);
