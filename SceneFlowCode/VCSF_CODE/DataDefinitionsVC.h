/*
Copyright (C) 2013 Christoph Vogel, PhD. Student ETH Zurich
All rights reserved.

This software is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3.0 of the License, or (at your option) any later version.

This software is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this software.*/

#ifndef _DATA_DEFINESVC_H
#define _DATA_DEFINESVC_H

/// turn openmp usage off -- also do not link the library
//#define _NO_OPENMP

/// maximal number of threads used:
#define _max_num_threads_used_ 4

/// compiler can handle sse4.1 ? - not important anyway
//#define __SSE4_1__

/// in case the computer can not handle sse(2) extensions
//#define _no_sse_ip_

/////////////// from DataPerSegment/Pixel: ///////////////
// stupid special case, if pixel moves behind camera what to do : impossible or oob
#define _impFactor_ 1.6

// special case, if pixel moves behind camera what to do : impossible or oob (yes) so off
//#define _inFrontPenalty_

// no influence 1.0 is on the image plane
#define _minDepth_ 1.0
//1.25

// only FW - no expansion possible
#define _expansion_ON_

#define _Shortcut_Speedup_

// indeed faster without but energy lower with
#define _withoutCuttingImprobable_
// commenting this leads to smaller expansion graphs constructed

// auto only for valid data terms (not oob, occ, imp) -- does nothing :) 
#define _checkValiadityAuto_

// test for data truncation -- actually gets worse? so off 
//#define _data_truncation_
// must be in conjunction with occlusion threshold !
#ifdef _data_truncation_
//#define cutOffValue 0.8
#endif

/// genscore.cpp:

/// 7x7 census - actually unfortunately hard-coded in here, although all data types are present which just need to be adjusted/checked
#define boxRadius 3 


// interacts with maxDisp and maxMot parameters and checks whether 2d motions/ disparities are in the specified range 
//
// the define turns off the check for maximal disparity and motion.
// the check prevents evaluations of the data term in places which do not make sense, so speeds up (sometimes significantly) the whole thing
// just imagine a bad proposals, covering the whole image in the other views ...
//
// ON KITTI on or off does nothing ! so really only for speedup, which it achieves.
// might be problematic for general camera configurations though
//#define __2dboundControl__
#define __depthControl__
#ifdef __depthControl__
  #define __depthControlAmplifier__  2.
#endif
//////////////////////////////////////////////

// penalize unreasonible motions (super fast mvps, geometry too close or behind camera)
// enable only for final improvement, very simple test to drop non-sense motions, better filter in advance to prevent eval of data ,etc. 
#define _MotionChecking_
//////////////////////////////////////////////

/// # of replacement steps applied -- only if _testLocalReplacement_ is defined
#ifdef _testLocalReplacement_
  #define __locRepsReplace__ 4
#endif

/// debugging:
//#define _simple_data_
//#define _FW_ONLY_

//using namespace std;
//using namespace Math;

/// defines whether the rescaling of oob pixel has to be done locally(per pixel) or globally per segment
//#define _globalRescaling_

/// debug output
//#define __QUIET__

/// use openmp or not
//#define _NO_OPENMP

/// cutoff data penalty at some point appears to be slightly worse - but faster if occlusion is on (less supermodular potentials)
//#define _data_truncation_

/// either LS-AUX or qpbo, the former is recommended so far
//#define _qpboVersion_

// if energy improves below this threshold stop -- in practice this too too low, speedup possible
//#define __Small_Improvement__ 0001
#define __Small_Improvement__ 0.501

/// debug: only depth, no motion - not sure if still working
//#define __do_depth_only__

/// 3d smoothing on/off (off: use 2d smoothing instead)
//#define __USE3D__

/// locally need to set something globally, enabling parts of vc version but placed here like ifdef ...
#define __vc_version__

#endif