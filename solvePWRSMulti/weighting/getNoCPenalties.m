%%% penalties for pixel without correspondence so oob, occ
function [autos] = getNoCPenalties(ref, cam, par, Seg, Seg2, useThreeFrames )

if numel(cam) ==1
    autoScores    = par.autoScale*ones(1, numel(Seg.Ids)+3*numel(Seg2.Ids));
    autoScoresPix = par.autoScale*ones( size(ref.I(1).I,1), size(ref.I(1).I,2), 4 );
else
    autoScores    = par.autoScale*ones(1, numel(Seg.Ids)+5*numel(Seg2.Ids));
    autoScoresPix = par.autoScale*ones( size(ref.I(1).I,1), size(ref.I(1).I,2), 6 );    
end

if useThreeFrames
    autoScores    = par.autoScale*ones(1, numel(Seg.Ids)+5*numel(Seg2.Ids));
    autoScoresPix = par.autoScale*ones( size(ref.I(1).I,1), size(ref.I(1).I,2), 6 );
end

autos.Seg = autoScores;
autos.Pix = autoScoresPix;
