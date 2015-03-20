function img = readFlowFile(filename)

% sanity check
if isempty(filename) == 1
    error('readFlowFile: empty filename');
end;

idx = findstr(filename, '.');
idx = idx(end);

if length(filename(idx:end)) == 1
    error('readFlowFile: extension required in filename %s', filename);
end;

if strcmp(filename(idx:end), '.sfl') ~= 1    
    error('readFlowFile: filename %s should have extension ''.sfl''', filename);
end;

if ~exist(filename, 'file')
  error('readFlowFile: could not find %s\n', filename);
end

fid = fopen(filename, 'r');
if (fid < 0)
    error('readFlowFile: could not open %s', filename);
end;

%tag     = fread(fid, 1, 'float32');
height   = fread(fid, 1, 'int32');
width    = fread(fid, 1, 'int32');

% sanity check

%if (tag ~= TAG_FLOAT)
%    error('readFlowFile(%s): wrong tag (possibly due to big-endian machine?)', filename);
%end;

if (width < 1 || width > 99999)
    error('readFlowFile(%s): illegal width %d', filename, width);
end;

if (height < 1 || height > 99999)
    error('readFlowFile(%s): illegal height %d', filename, height);
end;

nBands = 4;

% arrange into matrix form
tmp = fread(fid, inf, 'double');
%tmp = reshape(tmp, [width*nBands, height]);
%tmp = tmp';
sz = height*width;
img = zeros(height,width,4);
ttmp = zeros(width,height);
ttmp(:) = tmp(1 : sz*(nBands-3));
tttmp = ttmp;
%ttmp( tttmp < 0 ) = 100000.0; % infty
img(:,:,1) = single(ttmp');
ttmp(:) = tmp(sz*(nBands-3)+1 : sz*(nBands-2));
%ttmp( tttmp < 0 ) = 0.0; % no movement
img(:,:,2) = single(ttmp');
ttmp(:) = tmp(sz*(nBands-2)+1 : sz*(nBands-1));
%ttmp( tttmp < 0 ) = 0.0; % no movement
img(:,:,3) = single(ttmp');
ttmp(:) = tmp(sz*(nBands-1)+1 : sz*nBands);
%ttmp( tttmp < 0 ) = 0.0; % no movement
img(:,:,4) = -single(ttmp');

%img(:,:,1) = tmp(:, (1:width)*nBands-3);
%img(:,:,2) = tmp(:, (1:width)*nBands-2);
%img(:,:,3) = tmp(:, (1:width)*nBands-1);
%img(:,:,4) = tmp(:, (1:width)*nBands);
      
fclose(fid);
