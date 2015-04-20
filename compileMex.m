%%% compile mex files

% no openmp: -D_NO_OPENMP
% OPTIMFLAGS = /O2 /Oy- /DNDEBUG /openmp
% for more speed compile externally with other compiler

windows = 1;
 sseenabled = 1;

% vcsf eccv14
if sseenabled
if windows
  mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -D_SEG_ -I./SceneFlowCode/CPP -I./SceneFlowCode/maxflow-v3.03.src -I. -output ./solvePWRSMulti/mex/Seg_3SM ./SceneFlowCode/VCSF_CODE/MexInterface.cpp ./SceneFlowCode/maxflow-v3.03.src/graph.cpp ./SceneFlowCode/maxflow-v3.03.src/maxflow.cpp
  mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -D_PIX_ -I./SceneFlowCode/CPP -I./SceneFlowCode/maxflow-v3.03.src -I. -output ./solvePWRSMulti/mex/Pix_3SM ./SceneFlowCode/VCSF_CODE/MexInterface.cpp ./SceneFlowCode/maxflow-v3.03.src/graph.cpp ./SceneFlowCode/maxflow-v3.03.src/maxflow.cpp
  mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -D_testLocalReplacement_ -D_SEG_ -I./SceneFlowCode/CPP -I./SceneFlowCode/maxflow-v3.03.src -I. -output ./solvePWRSMulti/mex/Seg_3SM_locreplace ./SceneFlowCode/VCSF_CODE/MexInterface.cpp ./SceneFlowCode/maxflow-v3.03.src/graph.cpp ./SceneFlowCode/maxflow-v3.03.src/maxflow.cpp
else
% -march=native
  mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -lgomp -D_SEG_ -D_sse_linux_ -I./SceneFlowCode/CPP -I./SceneFlowCode/maxflow-v3.03.src -I. -output ./solvePWRSMulti/mex/Seg_3SM ./SceneFlowCode/VCSF_CODE/MexInterface.cpp ./SceneFlowCode/maxflow-v3.03.src/graph.cpp ./SceneFlowCode/maxflow-v3.03.src/maxflow.cpp
  mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -lgomp -D_PIX_ -D_sse_linux_ -I./SceneFlowCode/CPP -I./SceneFlowCode/maxflow-v3.03.src -I. -output ./solvePWRSMulti/mex/Pix_3SM ./SceneFlowCode/VCSF_CODE/MexInterface.cpp ./SceneFlowCode/maxflow-v3.03.src/graph.cpp ./SceneFlowCode/maxflow-v3.03.src/maxflow.cpp
  mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -lgomp -D_sse_linux_ -D_testLocalReplacement_ -D_SEG_ -I./SceneFlowCode/CPP -I./SceneFlowCode/maxflow-v3.03.src -I. -output ./solvePWRSMulti/mex/Seg_3SM_locreplace ./SceneFlowCode/VCSF_CODE/MexInterface.cpp ./SceneFlowCode/maxflow-v3.03.src/graph.cpp ./SceneFlowCode/maxflow-v3.03.src/maxflow.cpp
end  
else
  % in case of no sse extension and openmp on computers:
  mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -lgomp -D_SEG_ -D_no_sse_ip_ -I./SceneFlowCode/CPP -I./SceneFlowCode/maxflow-v3.03.src -I. -output ./solvePWRSMulti/mex/Seg_3SM ./SceneFlowCode/VCSF_CODE/MexInterface.cpp ./SceneFlowCode/maxflow-v3.03.src/graph.cpp ./SceneFlowCode/maxflow-v3.03.src/maxflow.cpp
  mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -lgomp -D_PIX_ -D_no_sse_ip_ -I./SceneFlowCode/CPP -I./SceneFlowCode/maxflow-v3.03.src -I. -output ./solvePWRSMulti/mex/Pix_3SM ./SceneFlowCode/VCSF_CODE/MexInterface.cpp ./SceneFlowCode/maxflow-v3.03.src/graph.cpp ./SceneFlowCode/maxflow-v3.03.src/maxflow.cpp
  mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -lgomp -D_testLocalReplacement_ -D_SEG_ -D_sse_linux_ -I./SceneFlowCode/CPP -I./SceneFlowCode/maxflow-v3.03.src -I. -output ./solvePWRSMulti/mex/Seg_3SM_locreplace ./SceneFlowCode/VCSF_CODE/MexInterface.cpp ./SceneFlowCode/maxflow-v3.03.src/graph.cpp ./SceneFlowCode/maxflow-v3.03.src/maxflow.cpp
end

if windows==1
  mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -llibmwlapack -I./SceneFlowCode/CPP  -I./../../CSource/SceneFlowCode/MatQ/Source/ -I. -output ./solvePWRSMulti/mex/gradGrowMatlabM9 ./SceneFlowCode/MatQ/Source/MexInterface.cpp
else
  mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -llapack -I./SceneFlowCode/CPP  -I./../../CSource/SceneFlowCode/MatQ/Source/ -I. -output ./solvePWRSMulti/mex/gradGrowMatlabM9 ./SceneFlowCode/MatQ/Source/MexInterface.cpp
end

% just deliver the ids per segments
mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -D_projectMVP_labels_ -I./SceneFlowCode/CPP -I./SceneFlowCode/VCSF_CODE/ -I./SceneFlowCode/VCSF_CODE/ProjectProposals -output ./solvePWRSMulti/mex/ProjectLabels ./SceneFlowCode/VCSF_CODE/ProjectProposals/MexInterfaceProjection.cpp
% with change of normals rotations
mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -D_projectMVPs_ -D_returnProjectionOnly_ -I./SceneFlowCode/CPP -I./SceneFlowCode/VCSF_CODE/ProjectProposals -I./SceneFlowCode/VCSF_CODE/ -output ./solvePWRSMulti/mex/ProjectMvps ./SceneFlowCode/VCSF_CODE/ProjectProposals/MexInterfaceProjection.cpp
mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -D_projectMVPs_grow_ -I./SceneFlowCode/CPP -I./SceneFlowCode/VCSF_CODE/ProjectProposals -I./SceneFlowCode/VCSF_CODE/ -output ./solvePWRSMulti/mex/ProjectCoverMvps ./SceneFlowCode/VCSF_CODE/ProjectProposals/MexInterfaceProjection.cpp

%%% fitting homographies with lorentzian:
mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -I./SceneFlowCode/CPP -output ./solvePWRSMulti/mex/fitHomoAllN ./SceneFlowCode/HomFitting/MexInterface_NormAll.cpp
mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -I./SceneFlowCode/CPP -output ./solvePWRSMulti/mex/fitHomoAllR  ./SceneFlowCode/HomFitting/MexInterface_RotAll.cpp
%mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -I./SceneFlowCode/CPP -output ./solvePWRSMulti/mex/fitHomoN ./HomFitting/SceneFlowCode/MexInterface_Norm.cpp
mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -I./SceneFlowCode/CPP -output ./solvePWRSMulti/mex/fitHomoR ./SceneFlowCode/HomFitting/MexInterface_Rot.cpp

if sseenabled
if windows==1
  mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -I./SceneFlowCode/CPP -output ./solvePWRSMulti/mex/Interpol_mex ./SceneFlowCode/BiCubicInterpolation/Source/mexInterface.cpp
else
  mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -lgomp -I./SceneFlowCode/CPP -output ./solvePWRSMulti/mex/Interpol_mex ./SceneFlowCode/BiCubicInterpolation/Source/mexInterface.cpp
end
  else
  mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -D_no_sse_ip_ -D_NO_OPENMP -I./SceneFlowCode/CPP -output ./solvePWRSMulti/mex/Interpol_mex ./SceneFlowCode/BiCubicInterpolation/Source/mexInterface.cpp
end

% % tvl1 on image - tgv2 would  be better but .. who cares? NEEDS -DWIN32 in windows
if windows==1
  mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -DWIN32 -I./SceneFlowCode/CPP -I./SceneFlowCode/TvL1_Matrix -output ./solvePWRSMulti/mex/TvL1Matrix ./SceneFlowCode/TvL1_Matrix/MexInterface.cpp
else
  mex CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" -O -I./SceneFlowCode/CPP -I./SceneFlowCode/TvL1_Matrix -output ./solvePWRSMulti/mex/TvL1Matrix ./SceneFlowCode/TvL1_Matrix/MexInterface.cpp
end
