###################################################################
#                                                                 #
#                 Piecewise Rigid Scene Flow                      #
#                                                                 #
###################################################################

Copyright 2013-2015 ETH Zurich (Christoph Vogel)


ABOUT:
This software implements our approach to scene flow estimation [1,2,3].


The additional packages
 - QPBO
 or alternatively
 - maxflow
are released under the GPL License and are not be included.
Please refer to http://pub.ist.ac.at/~vnk/software.html to download that 
package and read the licensing information provided there.




DISCLAIMER:
This demo software has been rewritten for the sake of simplifying the
implementation. Therefore, the results produced by the code may differ
from those presented in the papers [1,2,3]. In fact the results should be 
a bit better on the KITTI dataset http://www.cvlibs.net/datasets/kitti/.

IMPORTANT:
If you use this software you should cite the following in any resulting publication:

    [1] Piecewise Rigid Scene Flow
        C. Vogel, K. Schindler and S. Roth
        In ICCV, Sydney, Australia, December 2013
		
    [2] View-Consistent 3D Scene Flow Estimation over Multiple Frames, 
		C. Vogel, S. Roth, K. Schindler, 
		In: 13th European Conference on Computer Vision (ECCV), Zurich, Switzerland, 2014

    [3] 3D Scene Flow with a Piecewise Rigid World Model, 
	    C. Vogel, K. Schindler, S. Roth, 
    	In: International Journal of Computer Vision (IJCV), Feb. 2015
	
For initialisation we use the method described in (with a lot less iterations):

    [4] An Evaluation of Data Costs for Optical Flow
        C. Vogel, S. Roth and K. Schindler
        In GCPR, Saarbruecken, Germany, September 2013


INSTALLING & RUNNING

1.	Download and install maxflow v3.03 from 
	http://pub.ist.ac.at/~vnk/software.html 
	and place it into ./SceneFlowCode.
	Alternatively one can also download and the recent version of QPBO 
	from http://pub.ist.ac.at/~vnk/software.html
	and place it into ./cpp
	and adjust the definition file DataDefinitionVC.h in the folder
	SceneFlowCode/VCSF_CODE.
	Please note that the code was tested with version v1.31 only. 
	Matlab code for the non-linear conjugate gradients optimization 
	procedure 'minimize' by Carl Edward Rasmussen is already included. 
	Make sure your machine has the lapack library. If you use linux
	the install should be easy, in case of Windows you may have it 
	already as part of your Matlab installation. 
	
	
2.	Start MATLAB and run compileMex.m to build the utilities binaries.
	(This step can be omitted if you are using Windows 64 bit or Unix 64 bit.)
	Reset windows=1 or windows=0 in the file, depending on your OS.
	Make sure to have lapack installed already, otherwise you need to 
	install it. 
	Adjust the compiler flags accordingly for your purposes 
	(defaults should work in most cases). 
	Additional compiler switches can be looked up in the file 
	DataDefinitionVC.h in the folder SceneFlowCode/VCSF_CODE.

	
3.	From folder solvePWRSMulti run run_pwrs_red( xx ) - 
	example given as comment in the code.

	
CHANGES
	1.0		April 19, 2014	Initial public release
	2.0		March 11, 2015	Including code from [2,3] into the release
	

Solved Issues:	

	Multi-core extension is now supported for linux also. 
	
Basic usage:

1.	Make sure that you read in the correct calibration files and images.
    As an example some images from the KITTI and NREC datasets are provided 
    along with calibration.

2. 	There are the options to run the algorithm with 3 frames. 
	So far it is assumed that the files are given as scenenumber_frame.png 
	as in KITTI. One can use three frames with the assumption of constant 
	3D velocity, 
	option: 
	p.use3Frames   = true; 
	and also to utilize the past solution as proposals, option:
	p.usePrevProps = true; 
	The procedure should fall back to the standard 2 frame procedure in 
	case there is no video data available.

3. 	So far the proposals are taken from semi-global matching and [4] for the 
	flow part.
	One can use arbitrary procedures for that, e.g. sparse matches.   
	For that the function generateProposals has to be modified accordingly.
	Proposals can simply be concatenated to the proposal set: N_prop, RT_prop.
	So far it is assumed to have one proposal per segment.
	In the file pwrsf_v4.m similar proposals are merged, such that some 
	abuse of this constraint is possible (eg. demonstrated already in 
	generateProposals).
	The cpp-code however does not have such a (unnecessary) constraint but requires
	to have proposal-id and center per proposal as input - besides a list of moving 
	planes (normal/rigid motion pair). 
	Thus new proposal algorithms could also only provide normal,rigid motion and 
	a center and refrain form the one proposal per segment idiom.
   
4. 	Standard proposal generation is performed by fitting a motion per segment 
	and to compress the information to 1k proposals, this can be adjusted.
   
5. 	Depending on the data set the grid size (initial super-pixel) has to be set.
	Rule of thumb is to use ~1850 proposals per image, so for KITTI we use
	p.gridSize= 16. 
	Other scenes might need a different parameter here, eg. 12 for 
	smaller images. This might render the refinement step pointless.
	In the file pwrsf_v4.m the parameters 
	refineLoop = 1;% run refining based on loop default: ON  if 16x16 grid
	endlevel   = 8;% refinement in 2^-1 steps, startlevel = 16, 8, .., 
	endlevel should be set accordingly.
	Standard is to subdivide the grid once and half the grid size and 
	expansion area. Here endlevel = 8 ensures this in the standard setting, 
	going from 16x16 to 8x8 superpixels. 
	This parameter might need to be changed at a different initial gridsize.
	Or turned off if the initial size is small, e.g. 10x10 super-pixels
   
6.	The expansion area can be adjusted with the parameters   
	p.gx=8; % 8 kitti - can be adjusted to image size / relative size
	p.gy=5; % 5 kitti
	Here we expand around a proposal center 8 grid cells horizontally and 5
	cells vertically in each direction. This basically trades speed (small values)
	against accuracy (larger values) of the method. 
	In the standard procedure we guessed that proposals are not reasonable for 
	segments further away than 8 and 5 cells.

7.	The behaviour of other parameters are analyzed in [3].  

8.	Why does a function loadXXX for your data be written by yourself and called
	in line 213 of run_pwrs_red.m instead of me providing it?
	The problem is that the code does assume a certain frame numbering, by scene
	and by frame within the scene. Also calibration matrices are required, which 
	must be provided by you. Furthermore the code needs to load 
	(if 3 frames are used) also the past frame of both cameras along. 
	When reasoning is done over multiple frames the code assumes that data from 
	the per-segment solution of the previous frame is available. 
	To that end the scene-number and the frame-number is used, 
	thus these must be provided by you. You can use the function loadKittiFlow.m 
	in the folder io as a guideline to create your read procedure for your data. 
	
9.	If you cannot run the provided flow library, or your Computer has more the 2 
	cores, you can download the code from https://github.com/vogechri/DataFlow 
	and compile it respectively. 
	Please change the lines 74-78 in the file 
	'solvePWRSMulti/proposals/generateProposals.m' 
	from 'TGV_flowDouble' to 'TGV_flow', twice.
	
10.	The code already delivers the disparity at the second time frame required for
	KITTI 2015 benchmark. The main computing function pwrsf_v4.m returns it along
	disparity and flow. The last line in pwrsf_v4.m:
	'flow2d =  reconstruc2dFlowHom( ref, cam(1), N_lin, Rt_lin, SegNew, 0 );'
	composes the 2D output, and the 4th component given by
	'flow(:,:,4) = cam.I(2).u(:,:,1) - ref.I(2).u(:,:,1);'
	is exactly the disparity at timestep 1.
	In the same manner one can easily compute the flow or depthmap between each frame 
	in the model, eg., in the 3-frame model one could return the disparity in the past
	or the motion between 1st and 3rd time-step. Other possibilities include to 
	reconstruct the 3D motion. To do such you can look at reconstruc3DFlowHom.m.
