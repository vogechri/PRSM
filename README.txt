###################################################################
#                                                                 #
#                 Piecewise Rigid Scene Flow                      #
#      Christoph Vogel, Konrad Schindler and Stefan Roth          #
#                          ICCV 2013                              #
#                          ECCV 2014                              #
#                          IJCV 2015                              #
#                                                                 #
#            Copyright 2013-2015 Christoph Vogel                  #
#                                                                 #
###################################################################



ABOUT:
This software implements our approach to scene flow estimation [1,2,3].


The additional packages
 - QPBO
 or alternatively
 - maxflow
are released under the GPL License and are not be included.
Please refer to http://pub.ist.ac.at/~vnk/software.html to download that 
package and read the licensing information provided there.



==========================================================================
DISCLAIMER:
This demo software has been rewritten for the sake of simplifying the
implementation. Therefore, the results produced by the code may differ
from those presented in the papers [1,2,3].
In fact the results should be a bit better on the KITTI dataset 
(http://www.cvlibs.net/datasets/kitti/).
==========================================================================


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
	Alternatively one can also
    download and the recent version of QPBO from
	http://pub.ist.ac.at/~vnk/software.html
	and place it into ./cpp
	and adjust the definition file DataDefinitionVC.h in the folder
	SceneFlowCode/VCSF_CODE.
	QPBO is needed for the ICCV version of the code only.
	Please note that the code was tested with version v1.31 only.
	
2.	Start MATLAB and run compileMex.m to build the utilities binaries.
	(This step can be omitted if you are using Windows 64 bit or Unix 64 bit.)
	Adjust the compiler flags accordingly for your purposes (defaults should work in most cases). 
	Additional compiler switsches can be looked up in the file DataDefinitionVC.h in the folder SceneFlowCode/VCSF_CODE.

	
3.	From folder solvePWRSMulti run run_pwrs_red( xx ) - example given as comment in the code.

4.	To enable plotting go to http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig
	download the export_fig package and place it in the lowest folder
	You can further download the sc-package from 
	http://www.mathworks.com/matlabcentral/fileexchange/16233-sc-powerful-image-rendering
	and place it in the same folder.
   
CHANGES
	1.0		April 19, 2014	Initial public release
	2.0		March 11, 2015	Including code from [2,3] into the release
	
	
	
	
Basic usage:

1.	Make sure that you read in the correct calibration files and images.
	As an example some KITTI images are provided along with calibration.

2. 	There are the options to run the algorithm with 3 frames. So far it is assumed that
	the files are given as scenenumber_frame.png as in Kitti.  
	One can use three frames with the assumption of constant 3D velocity,
	option: 
	p.use3Frames   = true; 
	and also to utilize the past solution as proposals, option:
	p.usePrevProps = true; 
	The procedure should fall back to the standard 2 frame procedure in 
	case there is no video data available.

3. 	So far the proposals are taken from semi-global matching and [4] for the flow part.
	One can use arbitrary procedures for that, e.g. sparse matches.   
	For that the function generateProposals has to be modified accordingly.
	Proposals can simply be concatenated to the proposal set: N_prop, RT_prop.
	So far it is assumed to have one proposal per segment.
	In the file pwrsfMulti_simpler_v3.m similar proposals are merged, such that some 
	abuse of this constraint is possible (eg. demonstrated already in generateProposals).
	The cpp-code however does not have such a (unnecessary) constraint but requires
	to have proposal-id and center per proposal as input - besides a list of moving planes
	(normal/rigid motion pair). 
	Thus new proposal algorithms could also only provide normal,rigid motion and a center and
	refrain form the one proposal per segment idiom.
   
4. 	Standard proposal generation is performed by fitting a motion per segment and to 
	compress the information to 1k proposals, this can be adjusted.
   
5. 	Depending on the data set the grid size (initial super-pixel) has to be set.
	Rule of thumb is to use ~1850 proposals per image, so for kitti we use
	p.gridSize= 16. 
	Other scenes might need a different parameter here, eg. 12 for 
	smaller images. This might render the refinement step pointless.
	In the file pwrsfMulti_simpler_v3.m the parameters
	refineLoop   = 1;% run refining based on loop default: ON  if 16x16 grid
	endlevel     = 8;% refinement in 2^-1 steps, so startlevel= 16, 8 , .., endlevel
	should be set accordingly.
	Standard is to subdivide the grid once and half the grid size and expansion area.
	Here endlevel = 8 ensures this in the standard setting, going from 16x16 to 8x8
	superpixels. This parameter might need to be changed at a different initial gridsize.
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