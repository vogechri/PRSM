% ugly code
function G = System_Matrix(alpha1, alpha0, M,N, par, g)

% part 1 : setup the system matrix (sparse)
% order ux, uy, wxx, wxy, wyx,wyy, same for v

% missing: insert 0 rows at appropriate positions
test  = reshape(1:N*M, M,N)';
test1 = test(1:end,1:end-1);
test2 = test(1:end,2:end);

% dual positions - leave out first row, and first column:
testDX  = reshape(1:N*M, M,N);
testDX = testDX(1:end,2:end);
testDX = testDX(:)';
testDY  = reshape(1:N*M, M,N);
testDY = testDY(2:end,1:end)';
testDY = testDY(:)';

% gradient: U
guxids_x = cat( 1, [1:M*N-M], [M+1:M*N] );
guxids_y = cat( 1, testDX, testDX );
guxids_v = cat( 1, alpha1 * ones(1,M*N-M), -alpha1 * ones(1,M*N-M));
% use the standard w here, full image, start at 2nd column

gwxids_x = 2*M*N + [M+1:M*N];
gwxids_y = testDX;
gwxids_v = -alpha1 * ones(1,(M*N-M));
ypos = M*N;

guyids_x = cat(1,test1(:)', test2(:)');
guyids_y = ypos + cat( 1, testDY, testDY );
guyids_v = cat(1, alpha1*ones(1,numel(guyids_x)/2), -alpha1*ones(1,numel(guyids_x)/2));
% use the standard w here, full image, start at 2nd column

gwyids_x = 3*M*N + test2(:)';
gwyids_y = ypos+testDY;
gwyids_v = -alpha1 * ones(1,(M*N-N));
%
if ~exist('g','var')
  xids_1 = cat(1, guxids_x, gwxids_x);
  yids_1 = cat(1, guxids_y, gwxids_y);
  vids_1 = cat(1, guxids_v, gwxids_v);
  xids_2 = cat(1, guyids_x, gwyids_x);
  yids_2 = cat(1, guyids_y, gwyids_y);
  vids_2 = cat(1, guyids_v, gwyids_v);
  xids = cat(1, xids_1(:), xids_2(:));
  yids = cat(1, yids_1(:), yids_2(:));
  vids = cat(1, vids_1(:), vids_2(:));
  %    G = sparse(yids, xids, vids);
  if ~par.doStereo
    %%%%%%%%%%%%% Gradient V
    xids_1 = cat(1, N*M+guxids_x,   2*N*M+gwxids_x);
    yids_1 = cat(1, 2*N*M+guxids_y, 2*N*M+gwxids_y);
    vids_1 = cat(1, guxids_v, gwxids_v);
    xids_2 = cat(1, N*M+guyids_x, 2*N*M+gwyids_x);
    yids_2 = cat(1, 2*N*M+guyids_y, 2*N*M+gwyids_y);
    vids_2 = cat(1, guyids_v, gwyids_v);
    
    Uxids = cat(1, xids_1(:), xids_2(:));
    Uyids = cat(1, yids_1(:), yids_2(:));
    Uvids = cat(1, vids_1(:), vids_2(:));
    xids = cat(1, xids(:), Uxids(:));
    yids = cat(1, yids(:), Uyids(:));
    vids = cat(1, vids(:), Uvids(:));
  end
  
else % weighted, sort and combine by y index!
  
  idy = reshape(1:N*M, M,N);
  id_dx1 = [zeros(M,1), idy(:,1:end-1)];
  id_dx2 = [zeros(M,1), idy(:,2:end)];
  id_dy1 = [zeros(1,N); idy(1:end-1,:)];
  id_dy2 = [zeros(1,N); idy(2:end,:)];
  
  g1 = squeeze(g(:,:,1));
  g2 = squeeze(g(:,:,2));
  g3 = squeeze(g(:,:,3));
  g4 = squeeze(g(:,:,4));
  
  xidux = cat(1, id_dx1(:)', id_dx2(:)', id_dy1(:)', id_dy2(:)', 2*N*M+idy(:)');
  xiduy = cat(1, id_dx1(:)', id_dx2(:)', id_dy1(:)', id_dy2(:)', 3*N*M+idy(:)');
  xidvx = N*M+xiduy;%cat(1, id_dx1(:)', id_dx2(:)', id_dy1(:)', id_dy2(:)', 3*N*M+idy(:)');
  xidvy = N*M+cat(1, id_dx1(:)', id_dx2(:)', id_dy1(:)', id_dy2(:)', 4*N*M+idy(:)');
  
  yid = repmat( idy(:)', 5, 1);
  valx = alpha1*cat(1, g1(:)', -g1(:)', g2(:)', -g2(:)', -ones(1,N*M));
  valy = alpha1*cat(1, g3(:)', -g3(:)', g4(:)', -g4(:)', -ones(1,N*M));
  validity = cat(1, id_dx1(:)', id_dx2(:)', id_dy1(:)', id_dy2(:)', id_dy1(:)'+id_dx1(:)');
  validity =validity>0; % also remove 0 entries if 0 in valxes
  
  valvalx = valx.*validity;
  valvaly = valy.*validity;
  % invalidate whole rows:
  valvalx = sum(abs(valvalx(1:4,:)), 1) > 0;
  valvaly = sum(abs(valvaly(1:4,:)), 1) > 0;
  
  validateX = abs(valx)>0 & bsxfun(@and, valvalx, validity);
  validateY = abs(valy)>0 & bsxfun(@and, valvaly, validity);
  xidux = xidux(validateX);
  xiduy = xiduy(validateY);
  xidvx = xidvx(validateX);
  xidvy = xidvy(validateY);
  
  yidx = yid(validateX);
  yidy = yid(validateY);
  valx = valx(validateX);
  valy = valy(validateY);
  
  if ~par.doStereo
    xids = cat(1, xidux(:), xiduy(:), xidvx(:), xidvy(:));
    yids = cat(1, yidx(:), N*M+yidy(:),2*N*M+yidx(:),3*N*M+yidy(:));
    vids = cat(1, valx(:), valy(:), valx(:), valy(:));
  else
    xids = cat(1, xidux(:), xiduy(:));
    yids = cat(1, yidx(:), N*M+yidy(:));
    vids = cat(1, valx(:), valy(:));
  end
end


%    G = sparse(yids, xids, vids);
% add derivatives w.r.t. w now:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xpos = 2*N*M;% first w
ypos = 4*M*N;
% for 4 w fields, dx,dy needed
gwxids_x = xpos + cat( 1, [1:M*N-M], [M+1:M*N] );
gwxids_y = ypos + cat( 1, testDX, testDX );%+ cat( 1, [1:M*N-M], [1:M*N-M] );
gwxids_v = cat( 1, alpha0 * ones(1,M*N-M), -alpha0 * ones(1,M*N-M));
ypos = ypos+M*N;
xids_1 = gwxids_x;
yids_1 = gwxids_y;
vids_1 = gwxids_v;
% now derivative in y direction same variable w though
gwyids_x = xpos + cat(1,test1(:)', test2(:)');
gwyids_y = ypos + cat( 1, testDY, testDY );
gwyids_v = cat(1, alpha0 * ones(1,numel(testDY)), -alpha0*ones(1,numel(testDY)));

if ~par.doStereo
  xids = cat(1, xids, xids_1(:), gwyids_x(:), ...
    N*M+xids_1(:), N*M+gwyids_x(:),...
    2*N*M+xids_1(:), 2*N*M+gwyids_x(:),...
    3*N*M+xids_1(:), 3*N*M+gwyids_x(:));
  
  yids = cat(1, yids, yids_1(:), gwyids_y(:), ...
    2*N*M+yids_1(:), 2*N*M+gwyids_y(:),...
    4*N*M+yids_1(:), 4*N*M+gwyids_y(:),...,
    6*N*M+yids_1(:), 6*N*M+gwyids_y(:));
  
  vids = cat(1, vids, vids_1(:), gwyids_v(:), vids_1(:), gwyids_v(:), ...
    vids_1(:), gwyids_v(:), vids_1(:), gwyids_v(:));
else % shorter
  xids = cat(1, xids, xids_1(:), gwyids_x(:), ...
    N*M+xids_1(:), N*M+gwyids_x(:), 3*N*M+gwyids_x(:,end));
  
  yids = cat(1, yids, yids_1(:), gwyids_y(:), ...
    2*N*M+yids_1(:), 2*N*M+gwyids_y(:), 6*N*M+gwyids_y(:,end));
  
  vids = cat(1, vids, vids_1(:), gwyids_v(:), vids_1(:), gwyids_v(:), ...
    gwyids_v(:,end));
end

G = sparse(yids, xids, vids);

%ret = spTestR(yids, xids, vids, max(yids), max(xids), max(xids), 1:max(xids));
%ret2 = G * [1:max(xids)]';

