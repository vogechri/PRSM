%%% reduces the proposals set by combining proposals which appear very
%%% similarly in the set. Criteria is the reprojection of the image corners
%%% if these are within some distance, propsals are assumed the same
%%% and cut
function [sol, N_linR, Rt_linR, ia] = reduce_NRT_pix (ref,cam, N_lin, Rt_lin)

%{ 
% joint with discretized N,R,t
scala1 = 10000; % 10000
scalat = 1000; % 10000
scalaR = 10000; % 10000
%NRT_uni = cat( 2, round( scala1 * N_lin')/scala1, reshape(round(scala2*(Rt_lin))/scala2,16, size(Rt_lin,3))' );
NRT_uni = cat( 2, round( scala1 * N_lin')/scala1, reshape(round(scalaR*(Rt_lin(1:4,1:end-1,:)))/scalaR,12, size(Rt_lin,3))', reshape(round(scalat*(Rt_lin(1:4,end,:)))/scalat, 4, size(Rt_lin,3))' );
%}
scalaPix = 1/3;%1/3; % ? 1/x, x pixel deviation, values 1-4 pixel ?

% apply on some pixel measure diff in image plane -- still not getting
% the true pixel thus rather conservative pick two border pixel
% check if diff < 1 pixel or so if so merge
Homos = get_Homographies_cam(ref, cam, 1, N_lin, Rt_lin, 1 );%, Seg, center);
p2ds = cat(2, [1;1;1], [size(ref.I(1).I, 2); size(ref.I(1).I, 1); 1], [size(ref.I(1).I, 2); 1; 1], [1; size(ref.I(1).I, 1); 1]);

v1 = squeeze(p2ds(1,1) .* Homos(1:3,1,:) + p2ds(2,1) .* Homos(1:3,2,:) + p2ds(3,1) .* Homos(1:3,3,:));
v2 = squeeze(p2ds(1,2) .* Homos(1:3,1,:) + p2ds(2,2) .* Homos(1:3,2,:) + p2ds(3,2) .* Homos(1:3,3,:));
v1 = bsxfun( @rdivide, v1, v1(3,:) );
v2 = bsxfun( @rdivide, v2, v2(3,:) );

% joint with discretized N,R,t
%NRT_uni = cat( 2, NRT_uni, round( squeeze(v1*scalaPix)' )/scalaPix, round( squeeze(v2*scalaPix)' )/scalaPix);
% no discretization
NRT_uni = cat( 2, round( squeeze(v1*scalaPix)' )/scalaPix, round( squeeze(v2*scalaPix)' )/scalaPix);

% add other corners:
v1 = squeeze( p2ds(1,3) .* Homos(1:3,1,:) + p2ds(2,3) .* Homos(1:3,2,:) + p2ds(3,3) .* Homos(1:3,3,:) );
v2 = squeeze( p2ds(1,4) .* Homos(1:3,1,:) + p2ds(2,4) .* Homos(1:3,2,:) + p2ds(3,4) .* Homos(1:3,3,:) );
v1 = bsxfun( @rdivide, v1, v1(3,:) );
v2 = bsxfun( @rdivide, v2, v2(3,:) );
NRT_uni = cat( 2, NRT_uni, round( squeeze(v1*scalaPix)' )/scalaPix, round( squeeze(v2*scalaPix)' )/scalaPix);

% add second homography -- not needed ?!
%{
Homos = get_Homographies_cam(ref, cam, 1, N_lin, Rt_lin, 0 );%, Seg, center);
v1 = p2ds(1,1) .* Homos(1:3,1,:) + p2ds(2,1) .* Homos(1:3,2,:) + p2ds(3,1) .* Homos(1:3,3,:);
v2 = p2ds(1,2) .* Homos(1:3,1,:) + p2ds(2,2) .* Homos(1:3,2,:) + p2ds(3,2) .* Homos(1:3,3,:);
NRT_uni = cat( 2, NRT_uni, round( squeeze(v1*scalaPix)' )/scalaPix, round( squeeze(v2*scalaPix)' )/scalaPix);
%}

[c,ia,ic] = unique(NRT_uni, 'rows');

%c = c(:,1:19);
c = cat( 2, N_lin', reshape(Rt_lin, 16, size(Rt_lin,3))' );
c = c(ia,:); % could be average of similar ones

N_linR  = c(:,1:3)';
%Rt_linR = reshape( c(:,4:end)', 4,4,size(c,1) );
Rt_linR = reshape( c(:,end-15:end)', 4,4,size(c,1) );
sol = int32(ic-1);

% joint with discretized N,R,t 
%{
n3_row = N_linR(3,:);
n3_rowSource = N_lin(3,ia);
n3_row ( abs(n3_row) < 1/scala1 ) =  sign(n3_rowSource ( abs(n3_row) < 1/scala1 )) / (2*scala1);
N_linR(3,:) = n3_row;
%}