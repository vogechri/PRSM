%%% input Normals, Rotations and translation, Sgementation, lines on which
%%% the points lie (K_l ^-1 * pixel)
function flow =  reconstruc3DFlowHom( N_, RT_, Seg, p2d_, centerVersion )

[N, M, D] = size(p2d_);
flow      = zeros(N,M,4);

if ~exist('centerVersion', 'var')
  centerVersion = 0;
end

% only compute flow here
for i=1:min( numel(Seg.Ids), size(N_,2))
  normal = N_(1:3,i);
  points = cat( 2, p2d_(Seg.Ids{i}+1), p2d_(1+Seg.Ids{i} + N*M) , p2d_(1+Seg.Ids{i} + 2*N*M))';
  d = abs(1./(points'*normal));
  flow(Seg.Ids{i}+1) = d;
end

p3d = bsxfun(@times, p2d_, squeeze(flow(:,:,1)));

% center of rotation in center of patch
if centerVersion
  for i=1:numel(Seg.Ids)
    
    center = [Seg.Centers(i,:), 1];
    center3D = [center / abs(center*N_(1:3,i)), 0];
    
    points3D = cat( 2, p3d(1+Seg.Ids{i}), p3d(1+Seg.Ids{i} + N*M) , p3d(1+Seg.Ids{i} + 2*N*M))';
    points3D = cat(1, points3D, ones(1,size(points3D,2)));
    Rt = RT_(1:3,:, i);
    if ~isreal(Rt)
      fprintf('reconstrucFlowHom, Rt not real %d\n', i);
    end
    
    flow3d = Rt * bsxfun( @minus, points3D, center3D') - points3D(1:3,:);
    flow3d = bsxfun( @plus, flow3d, center3D(1:3)');
    flow(1+Seg.Ids{i} +   N*M) = flow3d(1,:);
    flow(1+Seg.Ids{i} + 2*N*M) = flow3d(2,:);
    flow(1+Seg.Ids{i} + 3*N*M) = flow3d(3,:);
  end
else
  % center of rotation not given (i.e. 0)
  for i=1:min( numel(Seg.Ids), size(RT_,3))
    
    points3D = cat( 2, p3d(1+Seg.Ids{i}), p3d(1+Seg.Ids{i} + N*M) , p3d(1+Seg.Ids{i} + 2*N*M))';
    points3D = cat(1, points3D, ones(1,size(points3D,2)));
    Rt = RT_(1:3,:, i);
    if ~isreal(Rt)
      fprintf('reconstrucFlowHom, Rt not real %d\n', i);
    end
    flow3d = Rt * points3D - points3D(1:3,:);
    
    flow(1+Seg.Ids{i} +   N*M) = flow3d(1,:);
    flow(1+Seg.Ids{i} + 2*N*M) = flow3d(2,:);
    flow(1+Seg.Ids{i} + 3*N*M) = flow3d(3,:);
  end
end