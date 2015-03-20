% ugly code
function [yids, xids, vids] = Matrix_data_TV( alpha1, M, N, par, mask, constraints, g)

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

%gwxids_x = 2*M*N + [M+1:M*N];
%gwxids_y = testDX;
%gwxids_v = -alpha1 * ones(1,(M*N-M));
ypos = M*N;

guyids_x = cat(1,test1(:)', test2(:)');
guyids_y = ypos + cat( 1, testDY, testDY );
guyids_v = cat(1, alpha1*ones(1,numel(guyids_x)/2), -alpha1*ones(1,numel(guyids_x)/2));
% use the standard w here, full image, start at 2nd column

%gwyids_x = 3*M*N + test2(:)';
%gwyids_y = ypos+testDY;
%gwyids_v = -alpha1 * ones(1,(M*N-N));
%

if exist('mask', 'var')
  mask_out_x = cat( 1, mask(1:M*N-M)|mask(M+1:M*N), mask(1:M*N-M)|mask(M+1:M*N) );
  mask = mask';
%  mask = false(size(mask));
  mask1 = mask(1:end,1:end-1);
  mask2 = mask(1:end,  2:end);
  mask1 = mask1 | mask2;
  mask_out_y = cat(1,mask1(:)', mask1(:)');
  
  guyids_x(mask_out_y) = [];
  guyids_y(mask_out_y) = [];
  guyids_v(mask_out_y) = [];
  guxids_x(mask_out_x) = [];
  guxids_y(mask_out_x) = [];
  guxids_v(mask_out_x) = [];
%  guxids_v(mask_out_x) = 0;
%  guyids_v(mask_out_y) = 0;
end

if ~exist('g','var')
  xids_1 = guxids_x;
  yids_1 = guxids_y;
  vids_1 = guxids_v;
  xids_2 = guyids_x;
  yids_2 = guyids_y;
  vids_2 = guyids_v;
  xids = cat(1, xids_1(:), xids_2(:));
  yids = cat(1, yids_1(:), yids_2(:));
  vids = cat(1, vids_1(:), vids_2(:));
  %    G = sparse(yids, xids, vids);
  if ~par.doStereo
    %%%%%%%%%%%%% Gradient V
    xids_1 = N*M+guxids_x;
    yids_1 = 2*N*M+guxids_y;
    vids_1 = guxids_v;
    xids_2 = N*M+guyids_x;
    yids_2 = 2*N*M+guyids_y;
    vids_2 = guyids_v;
    
    Uxids = cat(1, xids_1(:), xids_2(:));
    Uyids = cat(1, yids_1(:), yids_2(:));
    Uvids = cat(1, vids_1(:), vids_2(:));
    xids = cat(1, xids(:), Uxids(:));
    yids = cat(1, yids(:), Uyids(:));
    vids = cat(1, vids(:), Uvids(:));
  else
    % last entry:
    xids(end+1) = 2*max(xids);
    yids(end+1) = 2*max(yids);
    vids(end+1) = 0.000001;
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
  
  xidux = cat(1, id_dx1(:)', id_dx2(:)', id_dy1(:)', id_dy2(:)');
  xiduy = cat(1, id_dx1(:)', id_dx2(:)', id_dy1(:)', id_dy2(:)');
  xidvx = N*M+xiduy;%cat(1, id_dx1(:)', id_dx2(:)', id_dy1(:)', id_dy2(:)', 3*N*M+idy(:)');
  xidvy = N*M+cat(1, id_dx1(:)', id_dx2(:)', id_dy1(:)', id_dy2(:)');
  
  yid = repmat( idy(:)', 4, 1);
  valx = alpha1*cat(1, g1(:)', -g1(:)', g2(:)', -g2(:)');
  valy = alpha1*cat(1, g3(:)', -g3(:)', g4(:)', -g4(:)');
  validity = cat(1, id_dx1(:)', id_dx2(:)', id_dy1(:)', id_dy2(:)');
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
    % last entry:
    xids(end+1) = 2*max(xids);
    yids(end+1) = 2*max(yids);
    vids(end+1) = 0.000001;
  end
end

% here need to add 
if max(yids(:)) < 4*N*M
  xids(end+1) = 2*N*M;
  yids(end+1) = 4*N*M;
  vids(end+1) = 0.00000;
end
  
if exist('constraints', 'var') && isfield(constraints, 'weights')
% append constraints:
% find 4 neighs p_i of constriaint position

www = constraints.weights;

ff = cat( 1, floor (constraints.p(:,1))', floor (constraints.p(:,2))');
fc = bsxfun(@plus, ff, [0,1]');%cat( 1, floor (constraints.p(:,1))', ceil (constraints.p(:,2))');
cf = bsxfun(@plus, ff, [1,0]'); %cat( 1, ceil (constraints.p(:,1))', floor (constraints.p(:,2))');
cc = ff+1;%cat( 1, ceil (constraints.p(:,1))', ceil (constraints.p(:,2))');
% and their weights w_i - sum w_i = 1
wff = sum((ff-constraints.p').^2,1);
wfc = sum((fc-constraints.p').^2,1);
wcf = sum((cf-constraints.p').^2,1);
wcc = sum((cc-constraints.p').^2,1);
[~,eliminate] = max(cat(1,wff,wfc,wcf,wcc));

test = cat(3,ff,fc,cf,cc);

%lastrow = max(yids(:)); % different due to masks, but constraints assume this
lastrow = 4*N*M;%max(yids(:));

newxids=zeros(1,3*2*size(constraints.p, 1));
newyids=zeros(1,3*2*size(constraints.p, 1));
newvids=zeros(1,3*2*size(constraints.p, 1));

for i = 1: size(constraints.p, 1)
points = squeeze(test(:,i,:));
points(:, eliminate(i) ) = [];

ids     = sub2ind ([M,N], points(2,:), points(1,:))';
weights = cat(1,points, ones(1,3))\[constraints.p(i,:),1]';

% high constraints also on smoothness at positions:
% these are the posistions the point is present
%{
%[find(xids==ids(1)), find(xids==ids(2)), find(xids==ids(3)]
% the y ids where a point is present, now find all positions - ugh

theYs = unique(yids([find(xids==ids(1)); find(xids==ids(2)); find(xids==ids(3))]));
theids = find(ismember(yids, theYs));
vids( theids ) = vids( theids ) * www(i);
% v flow
theYs = unique(yids([find(xids==ids(1)+N*M); find(xids==ids(2)+N*M); find(xids==ids(3)+N*M)]));
theids = find(ismember(yids, theYs));
vids( theids ) = vids( theids ) * www(i);
%}
% too slow:
% xids = cat (1, xids, ids);
% yids = cat (1, yids, repmat( lastrow+2*i-1, [3,1]));
% vids = cat (1, vids, www(i)*weights );
% 
% xids = cat (1, xids, ids+N*M);
% yids = cat (1, yids, repmat( lastrow+2*i, [3,1]));
% vids = cat (1, vids, www(i)*weights );

newxids(6*(i-1)+1:6*i) = [ids;ids+N*M];
%newyids(6*(i-1)+1:6*i) = [repmat( lastrow+2*i-1, [3,1]); repmat( lastrow+2*i, [3,1])];
newyids(6*(i-1)+1:6*i) = [lastrow+2*i-1;lastrow+2*i-1;lastrow+2*i-1;lastrow+2*i;lastrow+2*i;lastrow+2*i];
newvids(6*(i-1)+1:6*i) = [www(i)*weights;www(i)*weights];

%union( union (find(xids==ids(1)), find(xids==ids(2))), find(xids==ids(3)))

end

xids = cat (1, xids, newxids');
yids = cat (1, yids, newyids');
vids = cat (1, vids, newvids');

% add a row to the matrix:
% sum w_i u_i = u
% weighted sum of their flow vectors is the constraint flow vector

end

%G = sparse(yids, xids, vids);
