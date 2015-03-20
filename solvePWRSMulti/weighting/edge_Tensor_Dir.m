%%% compute first approaxiamtion of a image gradinet based tensor
function g = edge_Tensor_Dir(f, addX, addY, direction)

q=1;alpha=5;

[M,N] = size(f);
idx = bsxfun(@plus, zeros(M,N),  1:N);
idy = bsxfun(@plus, zeros(M,N), (1:M)');

idx = idx + addX;
idy = idy + addY;

% cutoff values
remove = (idx>N) | (idx<1);% | (idy>M) | (idy<1);

idx(:,sum(remove,1)>0) = [];
idy(:,sum(remove,1)>0) = [];

remove = (idy>M) | (idy<1);

idx(sum(remove,2)>0,:) = [];
idy(sum(remove,2)>0,:) = [];

[~, fx, fy] = Interpol_mex(f, idx, idy);

fx=reshape(fx, size(idx));
fy=reshape(fy, size(idx));

norm = max(1e-06, sqrt(fx.^2 + fy.^2));
small = (norm <= 1e-06);

fx = fx ./ norm;
fy = fy ./ norm;

fx(small) = 1;
fy(small) = 0;

gg = max(1e-06, exp(-alpha*norm.^q));

g =     cat(3, gg .* fx.*fx, gg .* fx.*fy, gg .* fx.*fy, gg .* fy.*fy);
g = g + cat(3,       fy.*fy,      -fx.*fy,      -fx.*fy,       fx.*fx);

g = ...
direction(1) * (direction(1)*g(:,:,1) + direction(2)*g(:,:,2))+ ...
direction(2) * (direction(1)*g(:,:,3) + direction(2)*g(:,:,4));