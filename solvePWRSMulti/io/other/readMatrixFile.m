function [Pl, Pr, F, epi, Rl, Rr, Rot, Tra] = readMatrixFile(filename)

% sanity check
if isempty(filename) == 1
    error('readMatrixFile: empty filename');
end;

idx = findstr(filename, '.');
idx = idx(end);

if length(filename(idx:end)) == 1
    error('readMatrixFile: extension required in filename %s', filename);
end;

if strcmp(filename(idx:end), '.sfl') ~= 1    
    error('readMatrixFile: filename %s should have extension ''.sfl''', filename);
end;

fid = fopen(filename, 'r');
if (fid < 0)
    error('readMatrixFile: could not open %s', filename);
end;

Pl = zeros(3,4);
Pr = zeros(3,4);
F  = zeros(3,3);

FlipY = [1 0 0; 0 -1 0 ; 0 0 1];

%tag     = fread(fid, 1, 'float32');
%width   = fread(fid, 1, 'int32');
%height  = fread(fid, 1, 'int32');

% sanity check

%if (tag ~= TAG_FLOAT)
%    error('readFlowFile(%s): wrong tag (possibly due to big-endian machine?)', filename);
%end;

%if (width < 1 || width > 99999)
%    error('readFlowFile(%s): illegal width %d', filename, width);
%end;

%if (height < 1 || height > 99999)
%    error('readFlowFile(%s): illegal height %d', filename, height);
%end;

%nBands = 4;

% arrange into matrix form
tmp = fread(fid, inf, 'double', 'ieee-le');
fclose(fid);

Pl(:) = tmp(1:12);
Pr(:) = tmp(13:24);

Pl(3,3) = - Pl(3,3);
Pr(3,3) = - Pr(3,3);

%Pr
%3.980149 0.000000 -0.099703 -0.099504
%0.000000 -4.000000 0.000000 0.000000
%-0.398015 0.000000 -0.997029 -0.995037
%-1.990074 0.000000 -1.977138 0.024814

% YES INDEED !, But this is done better in the writer!!
%Pl(3,:) = [0 0 -1 0];
%Pr(3,:) = [0 0 -1 0];

% NO this does not apply here, done in the writer
%Pl = FlipY * Pl;
%Pr = FlipY * Pr;
%N = 256;
%M = 256;
%Viewport = [ (N-1)/2 0 (N+1)/2 ; 0 (M-1)/2 (M+1)/2 ; 0  0  1];

%Pl = Viewport * Pl;
%Pr = Viewport * Pr;

[A B C] = svd (Pl); % ! C is nor C' : F = A*B*C' ! -> epi2 is col of C

cl = C(:, end);
el = Pr * cl;
% DO NOT DO THAT !!! (computation of ur and so on)
%el = el/el(3);
%3.934566 -2.530993 2.277198
%test2 = Pr * [ 0.933594 0.996094 -4.000000 1]';
%[0.95463848114013672 -0.98581051826477051 -3.9899985790252686 1]';%[3.974414 2.518199 2.319126 1]';
%test3 = test2 / test2(3); 
%test4 = Pl * [ 0.933594 0.996094 -4.000000 1]';
%test5 = test4 / test4(3); 


F = Pr * ( Pl'*inv(Pl*Pl') );

%A = [[0 , -e(3), e(2)]; [ e(3), 0, -e(1)]; [ e(1) , e(2), e(3)]];
ex = [ [0 , -el(3), el(2)]; [ el(3), 0, -el(1)]; [-el(2) , el(1), 0]];
%F = [ 0  el(3) -el(2); -el(3) 0 el(1); el(2) -el(1) 0] * F;
F = ex * F;

F = F / norm(F,'fro');

[A B C] = svd (F); % ! C is nor C' : F = A*B*C' ! -> epi2 is col of C

test = el'*F;
epi = el;
Rl = Pl(:,1:3);
Rr = Pr(:,1:3);

[Kk, Rot] = rq(Rr);
%[A,B]=qr(inv(Rr));a=inv(A);B=inv(B);
Rot = Rot ./ Kk(3,3);
Kk = Kk ./ Kk(3,3);

Tra = Kk^-1 * Pr(:,end);
