/*
Copyright (C) 2013 Christoph Vogel, PhD. Student ETH Zurich

This software can be used for research purposes only.
This software or its derivatives must not be publicly distributed
without a prior consent from the author (Christoph Vogel).

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef _DATA_DEFINES_H
#define _DATA_DEFINES_H

/// compiler can handle sse4.1 ? - not important anyway
//#define __SSE4_1__

/// rotations stored wrt. center of patch == trick to receive more proposals without doing anything
#define _use_patchCenters_  

using namespace std;
using namespace Math;

/// defines whether the rescaling of oob pixel has to be done locally(per pixel) or globally per segment
//#define _globalRescaling_

/// debug output
//#define __QUIET__

/// use openmp or not
//#define _NO_OPENMP

/// cutoff data penalty at some point appears to be slightly worse - but faster if occlusion is on (less supermodular potentials)
//#define _data_truncation_

// penalize unreasonible motions (super fast, geometry too close or behind camera)
#define _motionChecking_	  	  

/// how many probing operations are done (unlabeled pixel)
#define maxProbeRuns 3
//3

/// super slow and useless if more nodes (probing if all are unlabeled does nothing)
#define limitUnsolved 1000

/// size of census data term - keep fixed (unsure if anything else but 3 is still working)
#define boxRadius 3

/// 3d smoothing on/off (off: use 2d smoothing instead)
//#define __USE3D__

/// use the approximation recommended in the paper when smoothing at superpixel level
#define _approximateSmoothing_

#endif