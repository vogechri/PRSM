% delivers errors on whole dataset reading the output files from a folder1
% and folder2 -- see run_pwrs_red parameter: storeFolder
function fullResult= evalSceneFlowKittiFolderThesis ( )

compareGT = 0;
numList = 0:193;

folder1 = 'C:\Users\vogechri\Desktop\work\results\validate\AllewOn_auto0_halfEvalEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_pwds_PS20_inFrontoff_max300_450_impOff_srand5_scala10000_fitPriorOn_autoI3Check';
% guess ego Prop == proposals from prev frame fitted to current prop-set
folder2 = 'C:\Users\vogechri\Desktop\work\results\autoMean\AllewOn_auto0_halfEvalEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_pwds_PS20_inFrontoff_max300_450_impOff_srand5_scala10000_fitPriorOn_autoI3_reOff_egoPropOn1';
%folder1 = 'C:\Users\vogechri\Desktop\work\results\autoMean\AllewOn_auto0_halfEvalEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_pwds_PS20_inFrontoff_max300_450_impOff_srand5_scala10000_fitPriorOn_autoI3on_reOff_egoPropOn1';
% WTF now pretty change only with auto vs without
folder1 = 'C:\Users\vogechri\Desktop\work\results\autoMean\AllewOn_auto0_halfEvalEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_pwds_PS20_inFrontoff_max300_450_impOff_srand5_scala10000_fitPriorOn_autoI3onNormal_reOff_egoPropOn1';
folder2 = 'C:\Users\vogechri\Desktop\work\results\autoMean\AllewOn_auto0_halfEvalEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_pwds_PS20_inFrontoff_max300_450_impOff_srand5_scala10000_fitPriorOn_autoI3offNormal_reOff_egoPropOn1';
% absurd: mean data energy is SHIT
%folder1 = 'C:\Users\vogechri\Desktop\work\results\autoMean\AllewOn_auto0_halfEvalEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_pwds_PS20_inFrontoff_max300_450_impOff_srand5_scala10000_fitPriorOn_autoI3offWeird_reOff_egoPropOn1';

% what is egoprop? -> prop from last wrt propset
%
% ahh: ? 2fps, auto, just 1 pset, reduced -- possible driven by solution of prev frame:
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_srand5';
%
%THE NORMAL CASE all off - 2 prop sets though REALLY BAD
% and unfortunately refinment:
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx0_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_srand5';
% this is indeed better ? WTF 900 props survive:
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx0_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover';

% 1200: survive -- worse in disp ?
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx0_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1200';
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx0_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover700';
% so only left -> patchCover better; 
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx0_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCoveroff_flowL';
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx0_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx0_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowL';
% here: LR better, also pathcCover does little worse
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx0_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowRL';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx0_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCoveroff_flowRL';

% use different ??
%folder1 = 'C:\Users\vogechri\Desktop\work\results\validate\AllewOn_auto0_halfEvalEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_pwds_PS20_inFrontoff_max300_450_impOff_srand5_scala10000_fitPriorOn_autoI3Check';
% this is refine, no add proposals, 
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\JournalECCV\AllewOn_auto0_halfEvalEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_pwds_PS20_inFrontoff_max300_450_impOff_srand5_scala10000_fitPriorOn_autoI3_reOff_egoPropOn1';


%folder2 = 'C:\Users\vogechri\Desktop\work\results\autoMean\AllewOn_auto0_halfEvalEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_pwds_PS20_inFrontoff_max300_450_impOff_srand5_scala10000_fitPriorOn_autoI3_reOff_egoPropOn1';
%folder1 = 'C:\Users\vogechri\Desktop\work\results\autoMean\AllewOn_auto0_halfEvalEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_pwds_PS20_inFrontoff_max300_450_impOff_srand5_scala10000_fitPriorOn_autoI3on_reOff_egoPropOn1';
% WTF now pretty change only with auto vs without
% THE SAME: auto on - like mine ? 
%folder1 = 'C:\Users\vogechri\Desktop\work\results\autoMean\AllewOn_auto0_halfEvalEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_pwds_PS20_inFrontoff_max300_450_impOff_srand5_scala10000_fitPriorOn_autoI3onNormal_reOff_egoPropOn1';
% auto off ?
%folder2 = 'C:\Users\vogechri\Desktop\work\results\autoMean\AllewOn_auto0_halfEvalEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_pwds_PS20_inFrontoff_max300_450_impOff_srand5_scala10000_fitPriorOn_autoI3offNormal_reOff_egoPropOn1';

folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx0_2fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowL';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_srand5';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx0_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowL';

%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLR';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLR';

%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCoveroff_flowLR';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCoveroff_flowLR';

folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1200_flowLRNew';
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover900_flowLRNew';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover900_flowLRNewpc';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx0_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowRL';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx0_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCoveroff_flowRL';

% CRAP
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover900_flowLRNewpc_weight';
% DOES NOTHING
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover900_flowLRNewpc_midW';

%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover900_flowLRNewpc_midW_weight01';
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover900_flowLRNewpc_midW_bug';
% lower epe higher rest == nothing
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover900_flowLRNewpc_bug_idOutlier';

% BAD
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover900_flowLRNewpc_bug_weightPrior';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover900_flowLRNewpc_bug_weightPriorOff';

folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1200_flowLRNewpc_bug';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug';
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_autoValidOff';
% ok not worse? no square or just superflous at all ? or no L1 ?
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_onealgCStep';
% hm no
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_towAlfCStepnoSquare';
% still ? 
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_towAlfCStep_noW';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_twoStepMine';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_twoStepMine_dw1';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_twoStepMine_noW';
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine';
%{
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine_TEST';
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine_TEST';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCoveroff_flowLR';

folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine_TESTp3f';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLR';

%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine_TESTp3f_loadF';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S1J1E1_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine_TESTp3f_loadF';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine080910_w4';
%     6.0294    3.8021    5.0933    4.7093
%     3.6237    2.3211    3.3526    3.1134
%     2.5722    1.6871    2.5649    2.3830
%     1.9907    1.3320    2.0876    1.9418
% 
% 6.0 & 3.6 & 2.6 & 2.0 & 3.8 & 2.3 & 1.7 & 1.3 & 5.1 & 3.4 & 2.6 & 2.1 & 4.7 & 3.1 & 2.4 & 1.9
% Worst 5 flow occ:
%  181 16.80, 018 13.50, 191 12.00, 062 11.90, 112 11.30, 140 10.80, 049 10.80, 156 10.40
% Worst 5 flow noc:
%  049 7.90, 191 7.80, 140 7.20, 112 7.20, 018 6.90, 092 6.80, 062 6.50, 156 6.40
% Worst 5 disp occ:
%  180 26.00, 071 14.40, 094 13.80, 049 12.30, 036 11.10, 120 10.90, 042 9.80, 161 9.70
% Worst 5 disp noc:
%  180 24.50, 071 14.50, 094 12.50, 049 11.60, 120 10.20, 036 9.80, 042 9.30, 161 8.60
% EPED: 0.717, 0.759 EPEF: 0.625 0.924
%}
% note loaded flow with 5 its ?  or what ?
% whti is it the same - because loading does use 2 different sets as well
% but my sequential code does not load shit thus .. 
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine_2inputs';
% csad: ^ %
% still gains - check also if only 1 input - 
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine_2inputs';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine_2inputs';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E1_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine_2inputs';
%%%

% better init flow but NOT better scene flow
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine_2inputsCen';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine_2inputsCenWon2';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine_2inputsCenWon3';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine_2inputsCenWon4';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine_2inputsCenWon4_s11';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine_2inputsCenWon4_s11a';
% 11.333 f 9.0 d
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_cen632_s12';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_cen632_s11';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_cen632_s10';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_cen128_s10';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_cen128_s11_eps2';
% check against census - lol normal better -> so off !
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_cen128_s11_eps1';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_cennormal_s11_eps1';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_csad_s11_eps1';

%%% testing different init proposals: csad and csad sgm - sgm with 80500
%%% better no kitti sgm
% CSAD vs TCEN as proposal -> CSAD better -- 
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine_2inputs';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_cen_full';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_cen_nosgm';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_csad_nosgm';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm';
% pic 140,79,19 tcensus SUCKS
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_onlysgmS100';
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_onlysgm_80500_S100';
% cant so it : only 74 still 
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_onlysgm_K80500_S100';
% what happened ?
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_csadsgm';

% the same sgmkitti and my old
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_csadsgm_K80500_S100';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_csadoldsgm_K80500_S100';

% not sure what to do here appears on par
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_csadoldsgmonly_K80500_S100';

folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J0E0_sw01_as075_tp015_ots015_csadonlys1j0e0';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_csadonlys0j0e0';
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_csadonly_dste';
%%%%

%%% here i test what i have to do for pwrsfinit -> unclear appears 03
%%% weighting is fine though
% 03 appears better 
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_egoWeight';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_egoWeight06';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_allWeight03';
% both the same WEIRD
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_allWeight03_noG';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_allWeight06_noG';

folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_allWeight03_withG';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_allWeight03_withoutG';
%%%

% here i check what happens if i try to fit 2 proposal sets instead of 1
% for sgm disp init : YES
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmlr_5001kprops';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmlr_200_500_1kprops';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmlr_100_200_500_1kprops';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_onlysgm_80500_S100';

% not so much any more
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_csadsgm_K80500_S100';

% WHAT:
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_csadonly_dste';

% WORSE ??? ok its only 400k
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_casdnewst_400_1kprops';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_casdnewst_400&1kprops';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_csad_nosgm';

% like in journal
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgm';
% guess : that sgm only other is csad and sgm and sgm with csad init 
% prefered -> only sgm ? YES since no Results_Fxx in folder
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm_8910';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm_8910_redo';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmfromcsad_8910_redo';
% what is the point here: ? why so good: 
% 1. loaded proposals ?
% 2. proposals failed ? -> load 3fr. proposals for test
% 3. what ever .. some parameters wrong check dates !
% worked on 20.2. 4am
% egomotionStep - never touched with s0j0e0 
% pwrsfMulti_simpler_v3, getDisparitySGM
% generateProposals_load, generateProposals, initSeg_2dFlowTest and algebraicMatrixQ9LR
% ok generateProposals_load does not work at all (no 08 files)
% thus generateProposals was used -- so ? 
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine080910_w4';

% i though this is it: .. but ? guess: tries to load from non-existent dir
% thus loaded non-present stuff with seg and r! prop
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine080910_w4';

% far worse even so never ran like that 
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_10_from3w';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_10_from4w';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_10_gp';
% ok works
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_10_gp8910';

% no gain ?????????
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S1J1E1_sw01_as075_tp015_ots015_s10';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_s10';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine_2inputs';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E1_sw01_as075_tp015_ots015_patchCover1000_flowLRNewpc_bug_oneStepMine_2inputs';

folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_csadonly_dste';
%{
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_csadonly_dste';

folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_csadonly_dste_8910';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_csadonly_dste';

%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm_8910';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmlr_200_500_1kprops';


folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_csadonly_dste_8910';

folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_csadonly_dste';

folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgm200';
% NO
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgm180';
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgm180_nostt';
%}

% no stt leasd to HORROR pic 74 overall better though 
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm180_nostt';
% pic 140 ALOT BETTER compared to above
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_nostt';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm224_nostt';

% with stt pic 74 ALOTBETTER
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_sttold';

%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_stt';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_stt50';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_sttNew100';

%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmlr_200_500_1kprops';

%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_csadonly_dste';

% stt no stt: nostt ?
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgm180';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgm180_nostt';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgm160';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm224';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm192';
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_nostt';
% {
% stt off worse here ?
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_8910';
% what now -- why? AAAAH ??? 74 here bad -- not sure -- or good ?
% no stt overall worse, no s1r1 pic 74 sucks
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm_8910';
% with stt
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_csadonly_dste_8910';
% consistently worse but worst are better?
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgm160_8910';
% with stt
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_csadonly_dste_8910';




folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_nostt';
folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm_8910';
% 47, 96 fits suck mc helps here
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC';
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp';

% early removal of bad motions - fails 
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp_rmLF';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp_autoSync';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp_control';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp_mcFlipC';

% keep 
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp';
% ok but is there a difference if it is so small ?
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp_mcFlip';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp_mcFlipflip';

% really exactly the same ?
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp_pubVersion_weird';
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp_pubVersion_nosortPix';

% alsoe worse ? wow
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp_pubVersion_srand_96on';
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp_pubVersion_srand1';
% WORSE
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp_pubVersion_srand';

% somehow with srand at seg level worse always
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp_pubVersion_96on';
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp_pubVersion_75on';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp_pubVersion_74on';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp_pubVersion_85rand';

folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp_pubVersion_85rand_removelargeMot';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp_pubVersion_85rand_AlllargeMotIN';

folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_Check';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_Check';
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_Checklm';%?
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_CheckMC';%?
% ok 0.06 but mc on still better for flow pics 8 and 96 !!!
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_CheckMC_pwrs06';%?
% still missing something
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_CheckMCon_pwrs06';
% data constraint to 250 and 300 displacements
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_CheckMCoff_d250m300_pwrs06';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_CheckMCoff_d250m300_pwrs06_inFrontPen';
% pic 141 FUCKED UP 155, 8
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_CheckMCoff_d350m400_pwrs06_frontCamPenNew';
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_CheckMC_pwrs06';%?
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_Check';


% mc off, then pwrs 06: better then stricter 250,300

% pic 188 not good overall worse then mcoff without
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_CheckMCoff_d350m400_pwrs06_frontCamPenNew_WHATs0r0_nom1_md0';
% really does not work -- WHY ? is it highly inferior, here smo 03
% 141 still ridiculously bad also 47
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_CheckMCoff_d350m400_pwrs06_frontCamPenNew_nom1';

% no control ? or pwrs: 0.06 smooth no really only noControl??????????????
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_CheckMCoff_d350m400_pwrs06_frontCamPenNew_noControl_nom1';

%%%% no control but depthcontrol
% ok
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_CheckMCoff_d350m400_pwrs06on_frontCamPenNew_noControl_nom1';
% really not ok: is it pwrs ? WTF maybe MC there ? or .. ??? 
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_CheckMCoff_d350m400_pwrs06on_frontCamPenNew_noControl_nom1';

% my new version sucks: odd behaviour: smo 05 better % epe SUCK pic 141 only
% smo06: % sucks epe ok
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_CheckMCoff_d350m400_pwrs06on_depthCtrl';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_CheckMCoff_d350m400_pwrs05_depthCtrl';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_CheckMCoff_d350m400_pwrs06on_depthCtrl';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_MCoff_d350m400_pwrs05_depthCtrl5';
% no instable
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_MCoff_d350m400_pwrs04_depthCtrl5';
% 141 and pic 8 loose it completely could compile once more with HIGHER pen ?
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_MCoff_d350m400_pwrs03_depthCtrl5';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_MCoff_d350m400_pwrs03_depthCtrl25';
% mc back: ok but too strict ? 
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_MCon_d350m400_pwrs03_depthCtrl5';
% non-sense result -- because of randomness little worse then above ?
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_d350m400_pwrs03_depthCtrl5';

% as = 0.2
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as065_tp015_ots015_sgmOnly_MCon2_d350m400_pwrs03_depthCtrl5_xtra02';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as065_tp015_ots015_sgmOnly_MCon2_d350m400_pwrs03_depthCtrl5_xtra02';

% NOT SO BAD ACTUALLY yet 3frames ? 
% in fact better +0.2 but still 0.75 so higher xtra pen -> better epe worse %
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_d350m400_pwrs03_depthCtrl5_xtra02';
% now old 0.1 xtraPen
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_d350m400_pwrs03_depthCtrl5_xtra01';
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_d350m400_pwrs03_depthCtrl5_xtra01';
% WTF: ? below fine thoguh
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01';
%%% now fine :
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_8910';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_8910';
%%%%%%%%%%%%%


% just mc off
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_Check';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_CheckMC_pwrs06';%?

% standard: MC ? 
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC_cutImp';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgm160';

%%% 3fps should be like this ? 
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_10_gp8910';

%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_csadonly_dste_8910';
% with stt
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_csadonly_dste_8910';
folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm_8910';
% 47, 96 fits suck mc helps here
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_noMC';
%%%

folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s045_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s045_S1J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01';

%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S1J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_fL';
% so only initTest remains as it is better !
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_fLTest';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_fLparts';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_fLR';
% too good ? barely a difference realistic? must be wrong
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_fLTest';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_fR';

folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_Replace23';
% only 0.1-0.05% on average 
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_ReplaceRef23';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_Replace34';
% roughly on par: 22 -- can i repeat ? 22 then 22 ? 
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_Replace22';
% NOT THE ANSWER FOR NO REASON ? still improves but ..
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_3fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_Replace22_rep2';
% NO
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_3fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_Replace22_rep3';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_3fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_Replace11_rep10Dup';
% this is with 3 times replace - lol and does nothing
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_3fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_refMore';

% compare with standard 8.4 .. NOTHING
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_3fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_refMore_norepl';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_d350m400_pwrs03_depthCtrl5_xtra01';

%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S1J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_Replace22';
% here replace but s0r0e0 and time series : so and so 
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_3fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_Replace22_8910';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_csadonly_dste_8910';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm_8910';

folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_3fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_red4stat';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_3fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_red8stat';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_3fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_red3stat';

% check : well .. ? 
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_3fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_redCheck';
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_3fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_redCheck8910';
%Worst 5 flow occ: 140 14.00, 156 13.40, 181 13.00, 062 11.70, 191 11.40, 018 11.10, 083 10.80, 049 10.10
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_3f_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_redCheck8910';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_2fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_redCheck8910';

% no idea seriously the randomness kills it to draw conclusions 
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_2f_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_depthCtrl5_xtra01_checkRepShort4';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_3fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_red4stat';
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_2f_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_depthCtrl5_xtra01_checkRepShort4_B';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_2f_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_depthCtrl5_xtra01_checkRepShort4_reDta';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_2f_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon5_depthCtrl5_xtra01_checkRepShort4_reDta';
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_3fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_red4stat';

folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_2f_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon8_depthCtrl5_xtra01_checkRepShort4';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_2f_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon8_depthCtrl5_xtra01_checkRepShort4B';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_2f_s05_S0J0E0_sw01_as085_tp015_ots015_sgmOnly_MCon8_depthCtrl5_xtra01_checkRepShort4C';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_2f_s05_S0J0E0_sw01_as08_tp015_ots015_sgmOnly_MCon8_depthCtrl5_xtra01_checkRepShort4C';
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_2f_s05_S0J0E0_sw01_as08_tp015_ots015_sgmOnly_MCon8_depthCtrl5_xtra01_checkRepShort4C_bug';
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_2f_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon8_depthCtrl5_xtra01_checkRepShort4_bug2';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_2f_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon8_depthCtrl5_xtra01_checkRepShort4_bug2_md0';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_2f_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCoff_depthCtrl5_xtra01_checkRepShort4_bug2_md0';
% well depth xx basically useless ?!
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_2f_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon_depthCtrl5off_xtra01_checkRepShort4_bug2_md0';
folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_2f_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon_depthCtrl1_xtra01_checkRepShort4_bug2';

folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto0_refx1_2f_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon_depthCtrl2_xtra01_checkRepShort4';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto0_refx1_3fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon_depthCtrl2_xtra01_checkRepShort4';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_3fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_redCheck8910';
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto0_refx1_2f_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon_depthCtrl2_xtra01_repShort4';
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_2f_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon_depthCtrl2_xtra01_repShort4';

%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_csadonly_dste_8910';
%Worst 5 flow occ: 181 12.20, 156 12.00, 062 11.80, 049 11.20, 018 11.00, 112 10.90, 191 10.70, 074 10.60
%folder1 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm_8910';

folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_ref';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_noref';

% not so much 0.2 %
folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\auto1_refx1_2fp_s05_S0J0E0_sw01_as075_tp015_ots015_norefine';

% no mc but smoothing
%folder1 = 'C:\Users\vogechri\Desktop\work\init\VCGIT\Auto1_refx1_2fp_lsaux_s055_S0J0E0_sw01_as075_tp015_ots015_mcOff';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_2f_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCoff_depthCtrl5_xtra01_checkRepShort4_bug2_md0';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\VCGIT\Auto1_refx1_2fp_lsaux_s065_S0J0E0_sw01_as075_tp015_ots015_mcOff';

folder1 = 'C:\Users\vogechri\Desktop\work\init\VCGIT\Auto1_refx1_2fp_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_mcOn';
folder2 = 'C:\Users\vogechri\Desktop\work\init\VCGIT\Auto1_refx1_2fp_lsaux_s045_S0J0E0_sw01_as075_tp015_ots015_mcOn';
%folder1 = 'C:\Users\vogechri\Desktop\work\init\VCGIT\Auto1_refx1_2fp_lsaux_s04_S0J0E0_sw01_as075_tp015_ots015_mcOn';

folder2 = 'C:\Users\vogechri\Desktop\work\init\VCGIT\Auto1_refx1_2fp_lsaux_s045_S0J0E0_sw01_as075_tp015_ots015_mcOn_r5';
% much better 0.045
folder1 = 'C:\Users\vogechri\Desktop\work\init\VCGIT\Auto1_refx1_3fps_lsaux_s045_S0J0E0_sw01_as075_tp015_ots015_mcOn_r5';

% lost it: ??
%folder1 = 'C:\Users\vogechri\Desktop\work\init\VCGIT\Auto1_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_mcOn_r5';
folder2 = 'C:\Users\vogechri\Desktop\work\init\VCGIT\Auto1_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_mcOn_r4';
folder2 = 'C:\Users\vogechri\Desktop\work\init\VCGIT\Auto1_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_mcOn8_r4';

folder1 = 'C:\Users\vogechri\Desktop\work\init\VCGIT\Auto1_refx1_2fp_lsaux_s045_S0J0E0_sw01_as075_tp015_ots015_mcOn_r5';
folder2 = 'C:\Users\vogechri\Desktop\work\init\VCGIT\Auto1_refx1_2f_s045_S0J0E0_sw01_as075_tp015_ots015_mcOn8_r4';


%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_3fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_redCheck8910';
% good as well : EPE !!! not sure though what this is : 
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto0_refx1_3fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon_depthCtrl2_xtra01_checkRepShort4';
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\Auto1_refx1_3fps_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_redCheck8910';

%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_csadonly_dste_8910';
%Worst 5 flow occ: 181 12.20, 156 12.00, 062 11.80, 049 11.20, 018 11.00, 112 10.90, 191 10.70, 074 10.60
%folder2 = 'C:\Users\vogechri\Desktop\work\init\Thesis\AllewOn_auto1_halfEEnergy_refx1_3fps_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm_8910';



%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgmOnly_MCon2_pwrs03_depthCtrl5_xtra01_fLR';

% mc off, pwrs06 ?
%folder1 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_CheckMC_pwrs06';%?
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2fp_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgmOnly_CheckMC';%?
% standard
%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S1J1E0_sw01_as075_tp015_ots015_sgm160';

%folder2 = 'C:\Users\vogechri\Desktop\work\init\vcclean\AllewOn_auto1_halfEEnergy_refx1_2f_lsaux_s05_S0J0E0_sw01_as075_tp015_ots015_sgm160_nostt';
% what is missing for vc ? replacement strategy -- or centered
% representation ? 

% what is different if centered? 
% test needs new smoothness and new computation of homos and normal per
% segment in future nad past frames
%
%replacement needs new optimization & data for whole frame

% test strekalowskies model with piecewise rigidity, argue coupling
% as a linearization to get rid of non-linearity, 
% like |u-f|^2 + G(f) , f = 2d flow
% u paramterized version as U(x,y) = H*(x,y,1)
% fixed u -> data flow step
% fixed f -> 2d fit but jointly as ? taylor data ? OR ??
% what happens here? -> fit and reassign to segments jointly
% affine model -> closed form ? 
% u is piecewise constant only if large alphas; depth could be done
% directly
%
% 
% penalize both
% motion == r,t and depth = normal ? seperately as depth can vary motion
% not
%
% 2. u is depth - let it very smoothly ? in segment ?
% leads to hard jumps between segments -- no since cutoff L2 regularizer
% so i need a square here ? 
% the problem is that |u-f|^2 + |nabla u|^2 u are parameter, then 
% parameter vary smooth, but |u-f|^2 needs F(u) instead F:u-> R^2
% displacements
% then first order F -> proxmap has solution ?
%
% now however u could be 1/depth; then 
% homo == K (R K^-1 p + t* 1/d * K^-1 p)
% thus can change t or d jointly -- linerization: appears fine 
% use 'disp' instead : 1/d replaced by d
% wilder model but regularizing disp difference appears better
% even linear function is ok'ish
% 
% tougher are R and t
% other idea is to regularize F(u) and couple 


% }
doDifference = 0;

% that image is flawed with gt
%  [a,b] = find(numList == 181);numList(b) = [];
%  [a,b] = find(numList == 140);numList(b) = [];
%  [a,b] = find(numList == 156);numList(b) = [];
%  [a,b] = find(numList == 74);numList(b) = [];
%  [a,b] = find(numList ==  1);numList(b) = [];
%  [a,b] = find(numList == 8);numList(b) = [];
%  [a,b] = find(numList == 96);numList(b) = [];
%  [a,b] = find(numList == 141);numList(b) = [];  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% at 3 pixel occ/nocev
allScores=[];

epes   = [0,0,0,0];
epes1  = [0,0,0,0];
epeOcc = [0,0];
epeNoc = [0,0];
epeOHit = [0,0];
epeNHit = [0,0];
allEpe1=[];

bestCombi = zeros(4);
listPic = [];
elms = 0;
sum = zeros(4);
elms2 = 0;
sum2 = zeros(4);
failImage = 0;
failImages = [];
for i=1:numel(numList)
  innerF = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  scoreFile  = dir(fullfile(folder1, innerF, sprintf('/RESULTS%03d*.txt', numList(i)) ));  
  
%  scoreFile  = dir(fullfile(folder1, innerF, sprintf('/RESULTS_Init_%03d*.txt', numList(i)) ));  
%  scoreFile  = dir(fullfile(folder1, innerF, sprintf('/RESULTS_LPL_%03d*.txt', numList(i)) ));
  scoreFile  = dir(fullfile(folder1, innerF, sprintf('/RESULTS_K%03d_10*.txt', numList(i)) ));
%  scoreFile  = dir(fullfile(folder1, innerF, sprintf('/RESULTS_F%03d_10_*.txt', numList(i)) ));  
  fName = fullfile(folder1, innerF, scoreFile.name);
  fid = fopen ( fName );
  if fid<0
    scoreFile  = dir(fullfile(folder1, innerF, sprintf('/RESULTS_%03d*.txt', numList(i)) ));
    fName = fullfile(folder1, innerF, scoreFile.name);
    fid = fopen ( fName );
  end 
  
%  scoreFile2 = dir(fullfile(folder2, innerF, sprintf('/RESULTS%03d*.txt', numList(i)) ));    
%  scoreFile2  = dir(fullfile(folder2, innerF, sprintf('/RESULTS_Rig_%03d*.txt', numList(i)) ));
  scoreFile2 = dir(fullfile(folder2, innerF, sprintf('/RESULTS_K%03d_10*.txt', numList(i)) ));
%  scoreFile2  = dir(fullfile(folder2, innerF, sprintf('/RESULTS_F%03d_10_*.txt', numList(i)) ));  
  fName2 = fullfile(folder2, innerF, scoreFile2.name);
  fid2 = fopen ( fName2 );
  if fid2<0
    scoreFile2  = dir(fullfile(folder2, innerF, sprintf('/RESULTS_%03d*.txt', numList(i)) ));
    fName2 = fullfile(folder2, innerF, scoreFile2.name);
    fid2 = fopen ( fName2 );
  end
  
  if fid > -1 && fid2 > -1
    try
      res = reshape( fscanf(fid, '%*s 2/3/4/5 %f & %f & %f & %f'), [4,4]);
    catch
      breakHere = 1;
      fclose( fid );fclose( fid2 );
      continue;
    end

    frewind(fid);
    for iiii=1:6 fgets(fid);end
    res4 = fscanf(fid, 'DispEPE %f & %f \nFlowEPE %f & %f');
%    res4 = fscanf(fid, ' EPE %f & EPE(noc) %f');res4 = cat(1, res4(:), res4(:));
    if ~isempty(res4)
      if any(isnan(res4))
        stophere = 1;
      else
        epes1 = epes1 + res4';
      end
    end
    try
      res2 = reshape( fscanf(fid2, '%*s 2/3/4/5 %f & %f & %f & %f'), [4,4]);
    catch
      breakHere = 1;
      fclose( fid );fclose( fid2 );
    end
    res = res(:,[2,4,1,3]);
    sum = sum + res;
    elms = elms +1;
    
    res2 = res2(:,[2,4,1,3]);
    sum2 = sum2 + res2;
    elms2 = elms2 +1;

    frewind(fid2);
    for iiii=1:6 fgets(fid2);end
    res3 = fscanf(fid2, 'DispEPE %f & %f \nFlowEPE %f & %f');
%    res3 = fscanf(fid, ' EPE %f & EPE(noc) %f');res3 = cat(1, res3(:), res3(:));

    if ~isempty(res3)
      if any(isnan(res3))
        stophere = 1;
      else
        epes = epes + res3';
      end
    end
    
    fclose( fid );fclose( fid2 );
    
    if res(2,2)< res2(2,2)
      bestCombi = bestCombi + res;
    else
      bestCombi = bestCombi + res2;
    end
    %%%%%%%%% differnece of methods :
    if doDifference == 1
      res = res2 - res;
      res4 = res3-res4;      
    end
    if doDifference == 2
      res = res - res2;
      res4 = res4-res3;      
    end    
    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%
    listPic(end+1)     =  numList(i);
    allScores(1,end+1) = res(1,1);
    allScores(2,end)   = res(2,1);
    allScores(3,end)   = res(3,1);
    allScores(4,end)   = res(4,1);
    allEpe1(1:4,end+1) = res4';    
%     allScores(1,end+1) = res2(1,1)-res(1,1);
%     allScores(2,end)   = res2(2,1)-res(2,1);
%     allScores(3,end)   = res2(3,1)-res(3,1);
%     allScores(4,end)   = res2(4,1)-res(4,1);    
    pos=4;
    for ii=2:4
      for iii = 1:4
        pos = pos + 1;
        allScores(pos,end) = res(iii,ii);
%        allScores(pos,end) = res2(iii,ii)-res(iii,ii);
      end
    end
    %%%%%%%%%%%%%
    
  else
    failImage = failImage+1;
    failImages(end+1) = numList(i);
    if fid > -1
      if compareGT
        res = reshape( fscanf(fid, '%*s 2/3/4/5 %f & %f & %f & %f'), [4,4]);
        res2 = res;
      end
      fclose( fid );
    end
    if fid2 > -1
      if compareGT
        res2 = reshape( fscanf(fid2, '%*s 2/3/4/5 %f & %f & %f & %f'), [4,4]);
        res = res2;
      end
      fclose( fid2 );
    end
    if compareGT
      sum2 = sum2 + res2;
      elms2 = elms2 +1;
      sum = sum + res;
      elms = elms +1;
    end
  end
end
fullResult1 = sum ./ elms;
fullResult1.*100

fullResult2 = sum2 ./ elms2;
fullResult2.*100
fprintf('Elements%03d\n', elms2);
% fprintf('%.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f\n', 100*fullResult1(:));
% fprintf('%.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f\n', 100*fullResult2(:));
% 
% bestCombi = bestCombi ./ elms;
% fprintf('%.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f\n', 100*bestCombi(:));

fprintf('%.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f\n', 100*fullResult1(:));
fprintf('%.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f\n', 100*fullResult2(:));
% 
bestCombi = bestCombi ./ elms;
fprintf('%.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f\n', 100*bestCombi(:));


%fprintf('Occ EPE: 1:%f 2:%f   ', epeOcc./epeOHit);fprintf('Noc EPE: 1:%f 2:%f\n', epeNoc./epeNHit);

endNr = min(7, numel(numList));
[Worst5 Worst5Order] = sort(allEpe1(4,:));
Worst5Order = listPic(Worst5Order);
fprintf('Worst 5 flowEPE occ:\n %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f\n', cat( 1, Worst5Order(end:-1:end-endNr), Worst5(end:-1:end-endNr) ));
[Worst5 Worst5Order] = sort(allEpe1(3,:));
Worst5Order = listPic(Worst5Order);
fprintf('Worst 5 flowEPE noc:\n %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f\n', cat( 1, Worst5Order(end:-1:end-endNr), Worst5(end:-1:end-endNr) ));
[Worst5 Worst5Order] = sort(allEpe1(2,:));
Worst5Order = listPic(Worst5Order);
fprintf('Worst 5 dispEPE occ:\n %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f\n', cat( 1, Worst5Order(end:-1:end-endNr), Worst5(end:-1:end-endNr) ));
[Worst5 Worst5Order] = sort(allEpe1(1,:));
Worst5Order = listPic(Worst5Order);
fprintf('Worst 5 dispEPE noc:\n %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f\n', cat( 1, Worst5Order(end:-1:end-endNr), Worst5(end:-1:end-endNr) ));
%%%
endNr = min(7, numel(numList));
[Worst5 Worst5Order] = sort(allScores(2,:));
Worst5Order = listPic(Worst5Order);
fprintf('Worst 5 flow occ:\n %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f\n', cat( 1, Worst5Order(end:-1:end-endNr), 100*Worst5(end:-1:end-endNr) ));
[Worst5 Worst5Order] = sort(allScores(6,:));
Worst5Order = listPic(Worst5Order);
fprintf('Worst 5 flow noc:\n %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f\n', cat( 1, Worst5Order(end:-1:end-endNr), 100*Worst5(end:-1:end-endNr) ));
[Worst5 Worst5Order] = sort(allScores(10,:));
Worst5Order = listPic(Worst5Order);
fprintf('Worst 5 disp occ:\n %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f\n', cat( 1, Worst5Order(end:-1:end-endNr), 100*Worst5(end:-1:end-endNr) ));
[Worst5 Worst5Order] = sort(allScores(14,:));
Worst5Order = listPic(Worst5Order);
fprintf('Worst 5 disp noc:\n %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f, %03d %.2f\n', cat( 1, Worst5Order(end:-1:end-endNr), 100*Worst5(end:-1:end-endNr) ));

epes  = epes./ elms;
epes1 = epes1./ elms;
fprintf('EPED: %.3f, %.3f EPEF: %.3f %.3f\nEPED: %.3f, %.3f EPEF: %.3f %.3f\n', epes1, epes);

return
fullResult = cat( 3, fullResult, fullResult1(:,:,1) );