echo "begin: `date`"
(
set -e
cd ~/software/vogechri_sceneflow/Geiger_TGV

# start always with a different one -> use info file for parameters, flowfile by number!
test=2
imgId=3
#la=50
la=15
wa=40
ig=0
fnr=10

#la=20# if old flow method

# all three receive the same img -id

for (( fnr=10; fnr<=10; fnr+=1))
do
for (( img=0; img<=193; img+=2))
  do
  for (( ring=2; ring<=2; ring+= 1))
    do
#    for (( stt=0.5; stt<=0.9; stt+=0.2))
#      do
#      for (( pf=0.5; pf<=0.95; pf+=0.2))
#        do
#        for (( eps=1.25 ; eps<=2 ; eps+=0.25))
#          do
#          for (( la=15; la<=150; la+=20))
#            do

# flow at t+1
#            bsub -W 0:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1, $la, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/CSAD_testF2Setr2_$test/', 'par_$img.inf', 194, 1, 1, 1, $ig, 0.5, 1, $ring, 1, 1);exit"

# stereo at t training
#            bsub -W 0:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1, 10, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGFlow_la10_noConv_ii6/', 'par_$img.inf', $img, 1, 1, 6, 0, 0.5, 1, $ring, 1, 0);exit"
#            bsub -W 0:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1, 5, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGFlow_la5_noConv/', 'par_$img.inf', $img, 1, 1, 3, 0, 0.5, 1, $ring, 1, 0);exit"
# stereo at t testing
#            bsub -W 0:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1, $la, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGFlowT/', 'par_$img.inf', $img, 1, 1, 3, 0, 0.5, 1, $ring, 1, 1);exit"
# stereo: la = 30
# flow at t
#            bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGCEN125_w40_5it_noEdge/', 'par_$img.inf', $img, 1, 1, 1, 0, 0.5, 1, $ring, 0, 0);exit"

#            bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1/255, 1.5, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGMI_w40_5it_noEdge_conv/', 'par_$img.inf', [$img:$img+9], 1, 1, 5, 0, 0.5, 0, $ring,0,0);exit"

#            bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1/255, 2, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGMI_w40_5it_noEdge_conv/', 'par_$img.inf', [$img:$img+9], 1, 1, 5, 0, 0.5, 0, $ring,0,0);exit"

#            bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1/255, 2, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGNCC_w40_5it_noEdge_conv/', 'par_$img.inf', [$img:$img+4], 1, 1, 5, 0, 0.5, 1, $ring, 0, 0);exit"

#            bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1/255, 1.5, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGNCC_w40_5it_noEdge_conv/', 'par_$img.inf', [$img:$img+4], 1, 1, 5, 0, 0.5, 1, $ring, 0, 0);exit"

#            bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1/255, 0.5, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGNCC_w40_5it_noEdge_conv/', 'par_$img.inf', [$img:$img+4], 1, 1, 5, 0, 0.5, 1, $ring, 0, 0);exit"

#            bsub -W 0:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 80, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGCEN125_w40_5it_noEdge/', 'par_$img.inf', $img, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0);exit"

#            bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGCEN125_w40_5it_noEdgeNew/', 'par_$img.inf', [$img:$img+7], 1, 1, 5, 0, 0.5, 1, $ring, 0, 0);exit"

# MULTIFRAME
#            bsub -W 35:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGCSAD_w40_5it_EdgeVideoSt_la10/', 'parA_$img.inf', [$img:$img+8], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 1, 0);exit"
#            bsub -W 35:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGCSAD_w40_5it_EdgeVideoR_la10/', 'parA_$img.inf', [$img:$img+8], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 0);exit"
#             bsub -W 35:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGCSAD_w40_5it_EdgeVideoL_la10/', 'parA_$img.inf', [$img:$img+8], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"

# how much is lost if faster, ring is 1, etc. 
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCSAD_w40_i5_la10_new/', 'parA_$img.inf', [$img:$img+8], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 40, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_w40_i5_la40_good/', 'parA_$img.inf', [$img:$img+8], $fnr, 1, 1, 5, 0, 0.5, 1, 1, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 30, 20, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_w20_i10_la30_good/', 'parA_$img.inf', [$img:$img+8], $fnr, 1, 1, 10, 0, 0.5, 1, 1, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 1.5, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_w40_i5_la1x5_linMrw1st/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, 1, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 3, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_w40_i5_la3_linMrw1st/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, 1, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 25, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_w40_i5_la25_FOpatch1stFilter/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, 1, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 15, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_w40_i5_la15_FOpatch1stFilterC/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, 1, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 110, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_w40_i5_la20_FOpatch1stFilterx2/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"

# SIMPLE
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 1.7, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_w40_i5_la1x7_SimpleI3/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 2.2, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_w40_i5_la2x2_SimpleAll/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"

# Good:
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 1.66, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCCleq0_GoodIAll_05Gone/', 'parA_$img.inf', [$img:$img+3], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 1.5, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_final2nd_05Gone/', 'parA_$img.inf', [$img:$img+3], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 12, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_finalPatch_New/', 'parA_$img.inf', [$img:$img+3], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 9*25/9, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_finalPatch_NewN0/', 'parA_$img.inf', [$img:$img+3], $fnr, 1, 1, 5, 0, 0.5, 1, 2, 0, 0, 1);exit"

# weighted gradient test:
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 9, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCSAD_weightedGrad/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, 1, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 7.5, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCSAD_weightedGrad/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, 1, 0, 0, 1);exit"

# NOT WEIGHTED CSAD:
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCSAD_notweightedGrad2/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, 1, 0, 0, 1);exit"


#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 8, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_finalPatch_TruncAfterBorderHesN0_noCut75/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, 1, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 70/9, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_finalPatch_TruncAfterBorderHesN0_noCut75/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, 1, 0, 0, 1);exit"

# SCARY
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 9.5, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_redo2_noMF/', 'parA_$img.inf', [$img:$img+7], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_redo2_peakF_oldSolveCensus/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
# little difference???
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_peakF_newSolveCensus_warpParam_newBorder_eq/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 8, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_finalPatch_TruncAfterBorderHesN0_noCut75_DOITAGAIN/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, 1, 0, 0, 1);exit"

# testSet: replace last 0 (second last parameter with 1)
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 20, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCSADCheckTGV_save_m1_TrainSet/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCSADCheckTGV_save_m1_TestSet/', 'parA_$img.inf', [$img:$img+2], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 1, 1);exit"

# stereo: a 1 after ring

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 0.15/255, 15, 20, 0.9, '/cluster/scratch_xl/public/vogechri/SparseM/TGVCenTrain_sparse_Flow_q1a5/', 'parA_$img.inf', [$img:$img+1], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 0.0/255, 9, 6, 0.9, '/cluster/scratch_xl/public/vogechri/SparseM/TGVCENTrain_sparse_Stereo_Gradbcq1a5/', 'parA_$img.inf', [$img:$img+1], $fnr, 1, 1, 17, 0, 0.5, 1, $ring, 1, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/SparseM/TGVCSADTrain_sparse_la15_StereoNewY_Gradbcq1a5/', 'parA_$img.inf', [$img:$img+1], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 1, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 15, 30, 0.9, '/cluster/scratch_xl/public/vogechri/SparseM/TGVCSADTrain_sparse_la15_flowRight/', 'parA_$img.inf', [$img:$img+3], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 0);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/SparseM/TGVCSADTrain_sparse_la15_flowLeft/', 'parA_$img.inf', [$img:$img+1], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"

# best lambda = ? 12.333
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 12.333, 5, 0.9, '/cluster/scratch_xl/public/vogechri/SparseM/GeigerFlow/pf6_a1/', 'parA_$img.inf', [$img:$img+1], $fnr, 1, 1, 20, 0, 0.5, 1, $ring, 0, 0, 1);exit"
             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 12.333, 5, 0.9, '/cluster/scratch_xl/public/vogechri/SparseM/GeigerFlow/pf6_a1_rr/', 'parA_$img.inf', [$img:$img+1], $fnr, 1, 1, 20, 0, 0.5, 1, $ring, 0, 0, 0);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 13, 5, 0.9, '/cluster/scratch_xl/public/vogechri/SparseM/GeigerFlow/', 'parA_$img.inf', [$img:$img+3], $fnr, 1, 1, 20, 0, 0.5, 1, $ring, 0, 0, 1);exit"
# works:
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 0.15/255, 9.333, 5, 0.9, '/cluster/scratch_xl/public/vogechri/SparseM/GeigerStereo/', 'parA_$img.inf', [$img:$img+1], $fnr, 1, 1, 20, 0, 0.5, 1, $ring, 1, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 0.15/255, 10, 5, 0.9, '/cluster/scratch_xl/public/vogechri/SparseM/GeigerStereo/', 'parA_$img.inf', [$img:$img+3], $fnr, 1, 1, 20, 0, 0.5, 1, $ring, 0, 0, 1);exit"

# pyr_fac, tgv, ring , stereo
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 0.15/255, 8.5, 5, 0.9, '/cluster/scratch_xl/public/vogechri/SparseM/TGVCenTrain_sparse_stereo_Flip_a0_5/', 'parA_$img.inf', [$img:$img+1], $fnr, 1, 1, 20, 0, 0.5, 1, $ring, 1, 0, 0);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 0.15/255, 10.5, 5, 0.9, '/cluster/scratch_xl/public/vogechri/SparseM/TGVCenTrain_sparse_stereo_Flip_a0_5/', 'parA_$img.inf', [$img:$img+1], $fnr, 1, 1, 20, 0, 0.5, 1, $ring, 1, 0, 0);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 20, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCSADTrain_sparse_la40_w20i10/', 'parA_$img.inf', [$img:$img+1], $fnr, 1, 1, 10, 0, 0.5, 1, $ring, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCSADCheckTGV_save_unWeighted_m1/', 'parA_$img.inf', [$img:$img+1], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 8*25/9, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCCCheckTGV_save/', 'parA_$img.inf', [$img:$img+1], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 8*25/9, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCCCheckTGV_save/', 'parA_$img.inf', [$img:$img+3], $fnr, 1, 1, 5, 0, 0.5, 1, 2, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 8*49/9, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCCCheckTGV_save/', 'parA_$img.inf', [$img:$img+3], $fnr, 1, 1, 5, 0, 0.5, 1, 3, 0, 0, 1);exit"

# for all rings:
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_peakF_newSolveCensus/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGV_CENSUSVerify_New_pf/', 'parA_$img.inf', [$img:$img+2], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_Cerify/', 'parA_$img.inf', 127, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_Cerify/', 'parA_$img.inf', 131, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_Cerify/', 'parA_$img.inf', 110, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_CVerify/', 'parA_$img.inf', 111, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_CVerify/', 'parA_$img.inf', 107, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_CVerify/', 'parA_$img.inf', 114, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_CVerify/', 'parA_$img.inf', 115, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_CVerify/', 'parA_$img.inf', 121, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_CVerify/', 'parA_$img.inf', 122, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_CVerify/', 'parA_$img.inf', 127, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_CVerify/', 'parA_$img.inf', 133, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_CVerify/', 'parA_$img.inf', 134, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_CVerify/', 'parA_$img.inf', 135, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_CVerify/', 'parA_$img.inf', 142, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_CVerify/', 'parA_$img.inf', 143, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_CVerify/', 'parA_$img.inf', 150, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_CVerify/', 'parA_$img.inf', 151, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_CVerify/', 'parA_$img.inf', 154, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVCensus_CVerify/', 'parA_$img.inf', 155, $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"


#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 8/9*25, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_finalPatch_TruncAfterBorderHesN0_noCut75/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, 2, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 8/9*49, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_finalPatch_TruncAfterBorderHesN0_noCut75/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, 3, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_finalPatch_TruncAfter/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, 1, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 25*7.5/9, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_finalPatch_NewA_N0/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, 2, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 49*7.5/9, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_finalPatch_NewA_N0/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, 3, 0, 0, 1);exit"


#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 8, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_finalPatch_NewA_N0/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, 1, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 9, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_finalPatch_NewN0/', 'parA_$img.inf', [$img:$img+3], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 8*25/9, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_finalPatch_New/', 'parA_$img.inf', [$img:$img+3], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 9*25/9, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_finalPatch_New/', 'parA_$img.inf', [$img:$img+3], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 1, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_finalPatch2ndB/', 'parA_$img.inf', [$img:$img+3], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"
#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 4, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_w40_i5_la4_Simple/', 'parA_$img.inf', [$img:$img+4], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0, 1);exit"

#             bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 1.5, 40, 0.9, '/cluster/scratch_xl/public/vogechri/TGVNCC_w40_i5_la50_linOld/', 'parA_$img.inf', [$img:$img+6], $fnr, 1, 1, 5, 0, 0.5, 1, 1, 0, 0, 1);exit"
#            bsub -W 5:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGCSAD_w40_5it_EdgeVideoSt_la10/', 'parA_$img.inf', [0,1,2,3,4,5,6,7,8,9,10], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 1, 0);exit"
#            bsub -W 5:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGCSAD_w40_5it_EdgeVideoR_la10/', 'parA_$img.inf',  [0,1,2,3,4,5,6,7,8,9,10], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0);exit"

#            bsub -W 7:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1.25/255, 10, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGCSAD_w40_5it_noEdgeVideoL/', 'parB_$img.inf', [138,147,172], $fnr, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0);exit"

#            bsub -W 0:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1/255, 10, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGStereo_la10_40_5it_noEdge_pre4/', 'par_$img.inf', $img, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0);exit"

#            bsub -W 0:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1/255, 7, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGStereo_la7_40_5it_noEdge_pre2/', 'par_$img.inf', $img, 1, 1, 5, 0, 0.5, 1, $ring, 0, 0);exit"

#            bsub -W 0:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1, 10, 25, 0.9, '/cluster/scratch_xl/public/vogechri/TVGStereo_la10/', 'par_$img.inf', $img, 1, 1, 20, 0, 0.5, 1, $ring, 0, 0);exit"
# flow testing
#            bsub -W 0:59 -R "rusage[mem=4048]" matlab -nosplash -nojvm -singleCompThread -r "run_KITTI( 1, $la, $wa, 0.9, '/cluster/scratch_xl/public/vogechri/TVGStereoT/', 'par_$img.inf', $img, 1, 1, 3, 0, 0.5, 1, $ring, 0, 1);exit"




#done
#echo "the Number is $imgId"
          imgId=`expr $imgId + 1`
          done
 #       done
#      done
#    done
#  done
#done
done
done
)
echo "end: `date`"
