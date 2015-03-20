/*
Copyright (C) 2014 Christoph Vogel, PhD. Student ETH Zurich
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.
*/

#ifndef __ENERGY_MOTIONCHECKER__CPP

#include "SensibleMotionCheck.h"

//#define __doubleCheck__

//////////////////////////////////////////////////////////////////////////////
// test for mvps which do not make sense. (super fast, in front of camera)  //
//////////////////////////////////////////////////////////////////////////////

template<class Scalar>
void MotionCheck<Scalar>::
  computeGeometryConsistencyLocal(int* segImg, int _nSegments, const std::vector<int> &ids_, 
  const std::vector<int> &seg2var, P3 mc, Scalar dm, 
  std::vector<Scalar>& _scores)
{
  Scalar dm2 = dm*dm;
  P3  addDispDiff = K * mc;
  // also assume that the max ids is smaller then the amount of normals
  assert(ids_.size() == nSegments);

  _scores.clear();
  _scores.resize( nSegments, 0);
  //    std::vector<Scalar> hits(nSegments, 0);
  const int step  = 4;
  const int step2 = step*step;

  for(int i=step/2; i< width;i+=step)
    for(int j=step/2; j< height;j+=step)
    {
      int pix = j+i*height;
      int segId   = segImg[pix];
      int localId = seg2var[segId];
      if( localId <0 ) continue;
      int k = ids_[segId];

      // get rotation, translation, normal
      const P3& norm  = (*normals)[k];// new
      const M3& rotk  = (*rot)[k];
      const P3& trak  = (*tra)[k];

      P3 add = trak;
      P3 p_3d = iK * P3( i+1, j+1, 1. );// a lot faster !!??

      Scalar dp = 1. / (norm | p_3d);
      p_3d      *= dp;

      P3 p_3d_t1 = rotk * p_3d + add;
      P3 motion  = p_3d_t1 - p_3d;

#ifndef __doubleCheck__
      Scalar distance2 = motion.sqrnorm() / std::min(dp*dp, std::max( Scalar(1), p_3d_t1[2]*p_3d_t1[2]));
#else
      Scalar distance2 = motion.sqrnorm() / std::min(dp*dp, std::max( Scalar(1), std::min( p_3d_t1[2]*p_3d_t1[2], (p_3d[2]-motion[2] )*(p_3d[2]-motion[2] ) )));
#endif

      // disp max compi:
      P3 p_img  = K * p_3d_t1;
      P3 p_img2 = p_img + addDispDiff;
      p_img /= p_img[2];p_img2 /= p_img2[2];
      Scalar dist = (p_img - p_img2).sqrnorm();
      //Scalar dist = (p_img - P3( ids[j]/height+1, ids[j]%height+1, 1. )).sqrnorm();

      if (distance2 > epsilon2 || dp < Scalar(1) || dist > dm2)
      {
        _scores[segId] += step2*penalty; // per segment
        //          hits[segId] += 1.0;
      }
    }
}



template<class Scalar>
void MotionCheck<Scalar>::
  computeGeometryConsistencyPerPixelLocal( int width, int height, int trial, P4i aoi, 
  const std::vector<int>& loc2glob, int* segImg, const std::vector<P3>* normalz, 
  std::vector<P2>& scoresUn, P3 mc, Scalar dm )
{
  const P3  normT  = (*normalz)[trial];
  const M3& rotkT  = (*rot)[trial];
  const P3& trakT  = (*tra)[trial];
  const Scalar dm2 = dm*dm;
  const P3 addDispDiff = K * mc;

  int aoix = aoi[0];
  int aoiX = aoi[2];
  int aoiy = aoi[1];
  int aoiY = aoi[3];

  //#pragma omp parallel for schedule (static)
  for (int ii=aoix;ii<aoiX;ii++)
  {
    int pos = ii*height;
    for (int j=aoiy;j<aoiY;j++)
    {
      int py = j;
      //        int i = pos + j;

      int id_pix1_local = loc2glob[j-aoiy + (aoiY-aoiy)*(ii-aoix)];
      if (id_pix1_local < 0) continue; // the other pixel have an id in the mrf tored here

      int px = ii;
      int segId = segImg[ pos + j ];
      //        int id_pix_global = pos + j;

      P3 p_2d = iK * P3( px+1, py+1, 1. );

      // get rotation, translation, normal
      const P3& norm  = (*normalz)[segId];
      const M3& rotk  =  (*rot)[segId];
      const P3& trak  =  (*tra)[segId];

      Scalar dp = 1. / (norm | p_2d);
      P3 p_3d   = p_2d*dp;

      P3 p_3d_t1 = rotk * p_3d + trak;
      P3 motion  = p_3d_t1 - p_3d;

      // disp max compi 2nd timestep
      P3 p_img  = K * p_3d_t1;
      P3 p_img2 = p_img + addDispDiff;
      p_img /= p_img[2];p_img2 /= p_img2[2];
      Scalar dist = (p_img - p_img2).sqrnorm();        

#ifndef __doubleCheck__
      Scalar distance = motion.norm() / std::min(dp, std::max( Scalar(1), fabs(p_3d_t1[2])));
#else
      Scalar distance = motion.norm() / std::min(dp, std::max( Scalar(1), std::min( fabs(p_3d_t1[2]), fabs(p_3d[2]-motion[2] )  )));
#endif

      P2 pen(0,0);

      if (distance > epsilon || dp < Scalar(1) || dist > dm2)
        pen[0] = penalty;

      ///////////////////////////////////////////////////////////////////

      dp      = 1. / (normT | p_2d);
      p_3d    = p_2d*dp;

      p_3d_t1 = rotkT * p_3d + trakT;
      motion  = p_3d_t1 - p_3d;

#ifndef __doubleCheck__
      distance = motion.norm() / std::min(dp, std::max( Scalar(1), std::min( fabs(p_3d_t1[2]), fabs(p_3d[2]-motion[2] )  )));
      //        distance = motion.norm() / std::min(dp, std::max( Scalar(1), fabs(p_3d_t1[2])));
#else
      distance = motion.norm() / std::min(dp, std::max( Scalar(1), std::min( fabs(p_3d_t1[2]), fabs(p_3d[2]-motion[2] )  )));
#endif

      // disp max compi:
      p_img  = K * p_3d_t1;
      p_img2 = p_img + addDispDiff;
      p_img /= p_img[2];p_img2 /= p_img2[2];
      dist = (p_img - p_img2).sqrnorm();

      if (distance > epsilon || dp < Scalar(1) || dist > dm2)
        pen[1] = penalty;

      scoresUn[ id_pix1_local ] = pen;
    }
  }
}


template<class Scalar>
void MotionCheck<Scalar>::
computeGeometryConsistency(const mxArray*& Segments, Scalar* centers_, const std::vector<int> ids_)
{
  nSegments = (int) (mxGetM(Segments) * mxGetN(Segments));
  const int skip = 8; // eval every skip'th pixel

  // also assume that the max ids is smaller then the amount of normals
  assert(ids_.size() == nSegments);

  scores.clear();
  scores.resize( nSegments, 0);

  for (int i = 0; i < nSegments ;i++)
  {
    int k = ids_[i];

    // get rotation, translation, normal
    const P3& norm = (*normals)[k];
    const M3& rotk  = (*rot)[k];
    const P3& trak  = (*tra)[k];

    // get center:
#ifdef _use_patchCenters_
    const P3 c2d( &centers_[ 3*i ] );
    P3 c3d = c2d * (1. / (norm | c2d));
    P3 add = trak + c3d -rotk * c3d;
#else
    P3 add(0.,0.,0.);
#endif
    int* ids = (int*) mxGetPr(mxGetCell(Segments, i));
    int nIds = (int)  mxGetNumberOfElements(mxGetCell(Segments, i));

    for (int j=0;j<nIds;j+=skip)
    {
      P3 p_3d = iK * P3( ids[j]/height+1, ids[j]%height+1, 1. );// a lot faster !!??
      //        P3 p_3d = P3(&p2d[3*ids[j]]);

      Scalar dp = 1. / (norm | p_3d);
      p_3d      *= dp;

      //          P3 p_3d_t1 = rotk * (p_3d - c3d) + trak + c3d;
      P3 p_3d_t1 = rotk * p_3d + add;
      P3 motion  = p_3d_t1 - p_3d;

#ifndef __doubleCheck__
      Scalar distance2 = motion.sqrnorm() / std::min(dp*dp, std::max( Scalar(1), p_3d_t1[2]*p_3d_t1[2]));
#else
      Scalar distance2 = motion.sqrnorm() / std::min(dp*dp, std::max( Scalar(1), std::min( p_3d_t1[2]*p_3d_t1[2], (p_3d[2]-motion[2] )*(p_3d[2]-motion[2] ) )));
#endif

      if (distance2 > epsilon2 || dp < Scalar(1))
        scores[i] += Scalar(skip)*penalty;

    }
  }
}


template<class Scalar>
void MotionCheck<Scalar>::
computeGeometryConsistencyPerPixel(int width, int height, int trial, int aoix, int aoiX, int aoiy, int aoiY, 
  std::vector<int>& loc2glob, int nPixel, int* segImg, Scalar* centers )
{
  scores.clear();
  scores.resize( 2*nPixel ,0);

  const P3& normT = (*normals)[trial];
  const M3& rotkT  = (*rot)[trial];
  const P3& trakT  = (*tra)[trial];
#ifdef _use_patchCenters_
  const P3  c2dT( &centers[ 3*trial ] );
  P3 c3dT = c2dT * (1. / (normT | c2dT));
  P3 addT = trakT + c3dT -rotkT * c3dT;
#else
  P3 addT(0.,0.,0.);
#endif

#pragma omp parallel for schedule (static)
  for (int j=aoiy;j<aoiY;j++)
  {
    int py = j;
    int pos = j*width;
    int i = pos+aoix;
    for (int ii=aoix;ii<aoiX;ii++,i++)
    {
      int id_pix1_local = loc2glob[ii-aoix + (aoiX-aoix)*(j-aoiy)];
      if (id_pix1_local < 0) continue; // the other pixel have an id in the mrf tored here

      int px = ii;
      int segId = segImg[ i ];
      int id_pix_global = i;
      //        P3 centerA ( &centers[3*segId] );

      P3 p_2d = iK * P3( px, py, 1. );

      // get rotation, translation, normal
      const P3& norm = (*normals)[segId];
      const M3& rotk  = (*rot)[segId];
      const P3& trak  = (*tra)[segId];

      // get center:
#ifdef _use_patchCenters_
      const P3 c2d( &centers[ 3*segId ] );
      P3 c3d = c2d * (1. / (norm | c2d));
      P3 add = trak + c3d -rotk * c3d;
#else
      P3 add (0.,0.,0.);
#endif
      //      P3 p_3d   = P3(&p2d[3*ids[j]]);
      Scalar dp = 1. / (norm | p_2d);
      P3 p_3d   = p_2d*dp;

      //          P3 p_3d_t1 = rotk * (p_3d - c3d) + trak + c3d;
      P3 p_3d_t1 = rotk * p_3d + add;
      P3 motion  = p_3d_t1 - p_3d;

#ifndef __doubleCheck__
      Scalar distance = motion.norm() / std::min(dp, std::max( Scalar(1), fabs(p_3d_t1[2])));
#else
      Scalar distance = motion.norm() / std::min(dp, std::max( Scalar(1), std::min( fabs(p_3d_t1[2]), fabs(p_3d[2]-motion[2] )  )));
#endif

      if (distance > epsilon || dp < Scalar(1))
        scores[ 2*id_pix1_local ] = penalty;

      ///////////////////////////////////////////////////////////////////

      dp = 1. / (normT | p_2d);
      p_3d    = p_2d*dp;

      //          P3 p_3d_t1 = rotk * (p_3d - c3d) + trak + c3d;
      p_3d_t1 = rotkT * p_3d + addT;
      motion  = p_3d_t1 - p_3d;

#ifndef __doubleCheck__
      distance = motion.norm() / std::min(dp, std::max( Scalar(1), fabs(p_3d_t1[2])));
#else
      distance = motion.norm() / std::min(dp, std::max( Scalar(1), std::min( fabs(p_3d_t1[2]), fabs(p_3d[2]-motion[2] )  )));
#endif

      if (distance > epsilon || dp < Scalar(1) )
        scores[ 2*id_pix1_local+1 ] = penalty;

    }
  }
}


#endif