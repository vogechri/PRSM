%%% optical flow / or stereo, see 
%%% An Evaluation of Data Costs for Optical Flow, 
%%% C. Vogel, S. Roth, K. Schindler, 
%%% In: German Conference on Pattern Recognition (GCPR), Saarbrücken, Germany, 2013
function flow = TGV_flow(cEps, lambda, warps, pyramid_factor, ...
  I1, I2, innerIts, ring, dataTerm, doStereo)

p.cEps     = cEps;
p.lambda   = lambda;
p.maxits   = innerIts;
p.ring     = ring;
p.doStereo = doStereo;
p.dataTerm = dataTerm;
p.warps    = warps;
p.pyramid_factor = pyramid_factor;

% print parameters
p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Linux_;
Linux_ = 0;

if ~usejava('jvm')
  Linux_ = 1;
end

global writeJava;
if Linux_ == 0
  writeJava = 1;
else
  writeJava = 0;
end

avail=license('checkout','Image_Toolbox');
while avail == 0
  %     disp('License not available, waiting...');
  pause(1.0);
  avail=license('checkout','Image_Toolbox');
end

  ticID = tic;
  flow = pyramidFlow(I1, I2, p, 16);
  totalTime = toc(ticID);

  fprintf(1, 'Elapsed total time: %f\n', totalTime);
end
