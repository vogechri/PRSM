function [invMatches] = findMatchesToInvalid(invalidRight, disp, valDisp, blockSize)
% [invMatches] = findMatchesToInvalid(invalidRight, disp, valDisp, blockSize)
%
% Given a set of disparities based on the left image of a stereo pair and a
% mask of invalid regions in the right image, finds disparities that
% correspond to matches to those invalid regions. Invalid regions will be
% grown by the given blockSize to also detect matches of blocks that
% include these pixels.
%
% Inputs:
%   invalidRight    Invalid regions in right image
%   disp            Disparities (indexed to left image)
%   valDisp         Validity mask for disp, true for every entry that is
%                   valid
%   blockSize       Block size to grow invalid regions by

if not(exist('valDisp', 'var')), valDisp = disp >= 0; end
if not(exist('blockSize', 'var')), blockSize = 15; end

% Expand artifacts in right image by blocksize to find where you might get
% bad matches
blockMatcher = strel('square', blockSize+2);
invRgt = imdilate(invalidRight, blockMatcher);
% For each match in left image, look up right image pixel it matched to
sz = size(disp);
xValsL = repmat(1:sz(2), [sz(1) 1]);
xValsR = xValsL;

xValsR(valDisp) = xValsL(valDisp) - round(disp(valDisp));
xValsRClipped = xValsR; xValsRClipped(xValsR<1) = 1;

yVals  = repmat((1:sz(1))', [1 sz(2)]);
indR   = sub2ind(sz, yVals, xValsRClipped);
invMatches = invRgt(indR);

