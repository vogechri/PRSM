function [flowFileName, errImg] = writeInfoFile(par, imgFolder, imgName, fSize, flow, flowGt, Energies)

storeFolder  = par.sFolder;
% unique name ?!
infoFileName = sprintf( '%s%s', imgName, par.infoFileName);
%flowFileName = sprintf( '%s_%s.ofl', flowFile, par.infoFileName);
flowFileName = sprintf( '%s/%s.ofl', storeFolder, infoFileName);
% sanity check %%%%%%%%%%%%%%%%%%
if isempty(infoFileName) == 1
    warning('writeInfoFile: empty filename');
    return;
end;

idx = strfind(infoFileName, '.inf');
if isempty(idx) == 1
infoFileName = sprintf('%s.inf', infoFileName);
  idx = strfind(infoFileName, '.');
end
  
idx = idx(end);

if length(infoFileName(idx:end)) == 1
    error('writewriteInfoFile: extension required in filename %s', filename);
end;

if strcmp(infoFileName(idx:end), '.inf') ~= 1
  infoFileName(idx:end)
    error('writeFlow: filename %s should have extension ''.inf''', filename);
end;

infoFileName = sprintf('%s/%s', storeFolder, infoFileName);
%fileName = infoFileName;
%%%%%%%%%%%%%%%


if ~exist(infoFileName,'file')

fid = fopen(infoFileName, 'w');
if (fid < 0)
  fprintf(2, 'writeInfoFile: could not open %s', infoFileName);
  return;
end;

%%% writing the paramenter
if isfield(par, 'lambda')
  fwrite(fid, sprintf('lambda %f\n',par.lambda), 'char');
end
if isfield(par, 'theta')
  fwrite(fid, sprintf('theta %f\n',par.theta), 'char');
end
if isfield(par, 'warps')
  fwrite(fid, sprintf('warps %f\n',par.warps), 'char');
end
if isfield(par, 'pyramid_factor')
  fwrite(fid, sprintf('pyramid_factor %f\n',par.pyramid_factor), 'char');
end
if isfield(par, 'useGradientImages')
  if isfield(par, 'wGradient')
    weightGrad = par.useGradientImages * par.wGradient;
  else
    weightGrad = par.useGradientImages;
  end
  fwrite(fid, sprintf('gradImgWeight %f\n', weightGrad), 'char');
end
if isfield(par, 'structure_texture_factor')
  fwrite(fid, sprintf('stt_factor %f\n',par.structure_texture_factor), 'char');
end
if isfield(par, 'innerItsD')
  fwrite(fid, sprintf('innerItsD %f\n',par.innerItsD), 'char');
end
if isfield(par, 'innerItsS')
  fwrite(fid, sprintf('innerItsS %f\n',par.innerItsS), 'char');
end
%if isfield(par, 'use_edges')
%  fwrite(fid, sprintf('use_edges %f\n',par.use_edges), 'char');
%end
if isfield(par, 'doMedian')
  fwrite(fid, sprintf('doMedian %f\n',par.doMedian), 'char');
end
if isfield(par, 'stockI')
  fwrite(fid, sprintf('stockI %f\n',par.stockI), 'char');
end
if isfield(par, 'gaussI')
  fwrite(fid, sprintf('gaussI %f\n',par.gaussI), 'char');
end
if isfield(par, 'epsilon')
  fwrite(fid, sprintf('epsilon %f\n',par.epsilon), 'char');
end
if isfield(par, 'reduz_fac')
  fwrite(fid, sprintf('reduz_fac %d\n', par.reduzFactor ), 'char');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fwrite(fid, sprintf('DESCRIPTION\n'), 'char');
else
  
fid = fopen(infoFileName, 'a');
if (fid < 0)
  fprintf(2, 'writeInfoFile: could not open %s', infoFileName);
  return;
end

end

fwrite(fid, sprintf('imageFolder %s\n',imgFolder), 'char');
fwrite(fid, sprintf('imageName %s\n',imgName), 'char');
fwrite(fid, sprintf('flowFile %s\n',flowFileName), 'char');
fwrite(fid, sprintf('fSize %d %d\n',fSize(1), fSize(2)), 'char');

if (size(flowGt,3) == 3)
  err2 = flow_error( flowGt, flow, 2);
  err3 = flow_error( flowGt, flow, 3);
  err4 = flow_error( flowGt, flow, 4);
  err5 = flow_error( flowGt, flow, 5);
  fwrite(fid, sprintf('EPE< 2/3/4/5: %f/%f/%f/%f\n', err2, err3, err4, err5), 'char');
  fprintf(sprintf('EPE< 2/3/4/5: %f/%f/%f/%f\n', err2, err3, err4, err5), 'char');
  [err, errImg] = getEndPointError(flow, flowGt);
else
if (numel(flowGt) == numel(flow))
  [err, errImg] = getEndPointError(flow, flowGt);
  [aae stdae aepe] = flowAngErr(flowGt(:,:,1), flowGt(:,:,2), flow(:,:,1), flow(:,:,2), 0); % ignore 0 boundary pixels
  fprintf('\nAAE %3.3f average EPE %3.3f \n', aae, aepe);
  
else
  if numel(flowGt) > 1
    [err, errImg] = getEndPointError_sparse(flow, flowGt);
  else
    err = nan;
    errImg = 0;
  end
end
end

fwrite(fid, sprintf('EPE: %f\n',err), 'char');
if exist('aae','var')
fwrite(fid, sprintf('AAE: %f\n',aae), 'char');
end
if exist('aepe','var')
fwrite(fid, sprintf('AEPE: %f\n',aepe), 'char');
end
  
if exist('Energies', 'var')
  fwrite(fid, sprintf('Energies: \n%s\n',Energies), 'char');
end 

%fwrite(fid, sprintf('flowFileGT %f\n',flowFileGT), 'char');

fclose(fid);
