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

#ifndef __ENERGY_ROTTRANORPIXEL__cpp
#define __ENERGY_ROTTRANORPIXEL__cpp

#include "EvalEnergyFull2FramePixel.h"

#include "DataDefinitionsVC.h"
#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

//#define __use_m128__
//#define _Debug_out_

#ifndef _NO_OPENMP
#include <omp.h>
#endif
#include<limits.h>
#include <map>
#include <vector>

using namespace std;
using namespace Math;

//////////////////////////////////////////////////////////////////////////////////////

template<class Scalar>
void
EvalEnergyFullFramePixel<Scalar>::
set2dMotionMatrix ( Scalar* K, Scalar* Rot, Scalar* mc, Scalar pixJump_ ) 
  { 
    Pl = M3(K);
    Pr = M3(K) * M3(Rot);
    pr = M3(K) * P3(mc);
  };

template<class Scalar>
void
EvalEnergyFullFramePixel<Scalar>::
set2dMotionMatrix ( Scalar* Kl, Scalar* Kr, Scalar* Rot, Scalar* mc, Scalar pixJump_ ) 
  { 
    Pl = M3(Kl);
    Pr = M3(Kr) * M3(Rot);
    pr = M3(Kr) * P3(mc);
  };

template<class Scalar>
void
EvalEnergyFullFramePixel<Scalar>::
set2dMotionMatrix_inv ( Scalar* K, Scalar* Rot, Scalar* mc, Scalar pixJump_ ) 
  { 
    Pl = M3(K);
    Pr = M3(K) * M3(Rot).invert();
    pr = Pr * (-P3(mc));
  };

template<class Scalar>
void
EvalEnergyFullFramePixel<Scalar>::
setK (Scalar* K_)
{
    assert (homs_r.size() < 1); // this is to be called before setHoms
    // both ways are identical !!!!
    iK=M3(K_);
    //    iK=M3(K_[0], K_[3], K_[6], K_[1], K_[4], K_[7], K_[2], K_[5], K_[8]);
    iK.invert();

    // the per pixel normals have to be multiplied with this here argh
    iiKt;iiKt = iK; iiKt=iiKt.transpose();iiKt.invert();
}

template<class Scalar>
void
EvalEnergyFullFramePixel<Scalar>::
setSegImg    ( int _w, int _h, int* _segImg, int _nSegments )
{ 
  height=_h;
  width =_w;
  nSegments = _nSegments;
  segImg    = _segImg;
};


/// segimg must be flipped
template<class Scalar>
void
EvalEnergyFullFramePixel<Scalar>::
prepare( int Npixel, int nSegments_)
{
    nSegments = nSegments_;
    Idk.clear();
    Idl.clear();
    F01.clear();
    F00.clear();
    F11.clear();
    F10.clear();
    F0.clear();
    F1.clear();

    if (halfEdgeXY == NULL)
    {
      F11.reserve( 2*Npixel );
      F01.reserve( 2*Npixel );
      F10.reserve( 2*Npixel );
      F00.reserve( 2*Npixel );
      Idk.reserve( 2*Npixel );
      Idl.reserve( 2*Npixel );
    }
    else
    {
      F11.reserve( 4*Npixel );
      F01.reserve( 4*Npixel );
      F10.reserve( 4*Npixel );
      F00.reserve( 4*Npixel );
      Idk.reserve( 4*Npixel );
      Idl.reserve( 4*Npixel );
    }

    F0.reserve( Npixel );
    F1.reserve( Npixel );
  };


/*!
  compute the energy - debug info
*/
template<class Scalar>
Scalar
EvalEnergyFullFramePixel<Scalar>::
compute_score_combiDepth_consistentScore( 
const std::vector<M3>& homs_r, const std::vector<M3>& homs_l, const std::vector<M3>& homs_rt )
  {
    int aoix = 0; int aoiX = width; int aoiy = 0; int aoiY = height;

    const int    dx8[8] = {-1,0,1, 0, 1, -1, -1,  1};
    const int    dy8[8] = { 0,1,0,-1, 1, -1,  1, -1};
    const Scalar we8[8] = {1.,1.,1.,1., 1./sqrt(2.), 1./sqrt(2.), 1./sqrt(2.), 1./sqrt(2.)};//weights

    const int neighIts = (halfEdgeXY == NULL) ? 4:8;
    const int* dx= dx8;
    const int* dy= dy8;
    const Scalar* we= we8;

    assert (normals != NULL);
    assert (tra != NULL);
    assert (rot != NULL);
    assert (homs_r.size() > 1); // this is to be called before setHoms
    Scalar score (0.);

    for (int j=aoiy;j<aoiY;j++)
    {
      Scalar f00_it(0), f01_it(0), f10_it(0), f11_it(0);
      int py = j;
      for (int ii=aoix;ii<aoiX;ii++)
      {
        int px = ii;
        int id_pix1_global = py+height*px;

        int idN_a = segImg[ py+height*px ];
        P3 norm0 = (*normals)[idN_a];
        Scalar dn_i =  1. / sqrt(norm0|norm0);
        M3 roti = (*rot)[idN_a];

        P3 n0  = norm0 * dn_i;
        P3 nr0 = roti * n0;

        assert(px<width);
        assert(py<height);
        assert(px>=0);
        assert(py>=0);
        assert((*normals).size() > idN_a);
        assert((*rot).size() > idN_a);
        assert(homs_r.size() > idN_a); 
        assert(homs_l.size() > idN_a); 
        assert(0 <= idN_a); 

        // run over the 4/8 neighbours, check whether oob or inside
        for(int displace = 0;displace < neighIts; displace++)
        {
          int qx = px+dx[displace];
          int qy = py+dy[displace];

          // oob? - indeed pixel can have term with border pixel
          if (qx<0 || qx>=width || qy<0 || qy>=height) continue;

          Scalar w_it = we[displace];

          // can also be -1
          int id_pix2_global = qx*height+qy;//qx+width*qy

          // oob go further (id_pix2_local=-1) else only 1 occurence
          if(id_pix2_global > id_pix1_global ) continue; // id_pix2_local == -1 still ok

          ////////////////////////////////////////////
          Scalar edgeWeight(1);

          if (displace<4)
          {
            // one of those 2 craps
            if(abs (id_pix2_global - id_pix1_global) != 1)
              edgeWeight = halfEdgeX[ min(qx,px)*height + qy ];
            else
              edgeWeight = halfEdgeY[ px*(height-1) + min(qy,py) ];
          }
          else
          {
            if (displace<6)
              edgeWeight = halfEdgeXY [ min(qx,px)*(height-1) + min(qy,py) ];
            else
              edgeWeight = halfEdgeiXY[ min(qx,px)*(height-1) + min(qy,py) ];
          }

          int idN_b = segImg[ id_pix2_global ];

          if ( idN_a == idN_b ) continue;
          assert(0 <= idN_b);
          //---------------------------------------------------------------------
          P3 norm1 = (*normals)[idN_b];     
          Scalar dn_j =  1. / sqrt(norm1|norm1);

          M3 rotj = (*rot)[idN_b];
          P3 n1  = norm1 * dn_j;
          P3 nr1 = rotj * n1;

          P3 p2d ( (px+qx+dy[displace])/2.+1., (py+qy+dx[displace])/2.+1., 1. );
          P3 q2d ( (px+qx-dy[displace])/2.+1., (py+qy-dx[displace])/2.+1., 1. );

#ifdef _Debug_out_
          printf("edge: %f,%f, %f,%f \n", p2d[0], p2d[1],q2d[0], q2d[1] );
#endif
          // now add the normal constraint to get a linear prior
          P3 dN00 = n0 - n1;

          P3 dN00r = nr0 - nr1;

          P2 f00, f11;
          P2 g00, g11;
          P2 h00, h11;

          P2 f00_t, f00_rt, f11_t, f11_rt; // f01, f01_t, f01_rt, f10, f10_t, f10_rt,
          P2 g00_t, g00_rt, g11_t, g11_rt; // g01, g01_t, g01_rt, g10, g10_t, g10_rt,
          P2 h00_t, h00_rt, h11_t, h11_rt; // h01, h01_t, h01_rt, h10, h10_t, h10_rt,

          P2 p_i_t0_r;
          P2 q_i_t0_r;
          P2 p_i_t1_l, p_i_t1_r;
          P2 q_i_t1_l, q_i_t1_r;

          P2 p_j_t0_r;
          P2 q_j_t0_r;
          P2 p_j_t1_l, p_j_t1_r;
          P2 q_j_t1_l, q_j_t1_r;

          P2 p_k_t0_r;
          P2 q_k_t0_r;
          P2 p_k_t1_l, p_k_t1_r;
          P2 q_k_t1_l, q_k_t1_r;

          P2 p_t0_l (p2d[0], p2d[1]);
          P2 q_t0_l (q2d[0], q2d[1]);

          P3 tmp1 = homs_r[idN_a]  * p2d;p_i_t0_r = P2(tmp1[0], tmp1[1]) * (1./tmp1[2]);
          P3 tmp2 = homs_r[idN_a]  * q2d;q_i_t0_r = P2(tmp2[0], tmp2[1]) * (1./tmp2[2]);
          P3 tmp3 = homs_l[idN_a]  * p2d;p_i_t1_l = P2(tmp3[0], tmp3[1]) * (1./tmp3[2]);
          P3 tmp4 = homs_l[idN_a]  * q2d;q_i_t1_l = P2(tmp4[0], tmp4[1]) * (1./tmp4[2]);
          P3 tmp5 = homs_rt[idN_a] * p2d;p_i_t1_r = P2(tmp5[0], tmp5[1]) * (1./tmp5[2]);
          P3 tmp6 = homs_rt[idN_a] * q2d;q_i_t1_r = P2(tmp6[0], tmp6[1]) * (1./tmp6[2]);

          tmp1 = homs_r[idN_b]  * p2d;p_j_t0_r = P2(tmp1[0], tmp1[1]) * (1./tmp1[2]);
          tmp2 = homs_r[idN_b]  * q2d;q_j_t0_r = P2(tmp2[0], tmp2[1]) * (1./tmp2[2]);
          tmp3 = homs_l[idN_b]  * p2d;p_j_t1_l = P2(tmp3[0], tmp3[1]) * (1./tmp3[2]);
          tmp4 = homs_l[idN_b]  * q2d;q_j_t1_l = P2(tmp4[0], tmp4[1]) * (1./tmp4[2]);
          tmp5 = homs_rt[idN_b] * p2d;p_j_t1_r = P2(tmp5[0], tmp5[1]) * (1./tmp5[2]);
          tmp6 = homs_rt[idN_b] * q2d;q_j_t1_r = P2(tmp6[0], tmp6[1]) * (1./tmp6[2]);
          ////////////////////////////  

          get2dDiff_( p_t0_l, p_i_t1_l, p_i_t0_r, p_i_t1_r,   p_t0_l, p_j_t1_l, p_j_t0_r, p_j_t1_r,  f00, f00_t, f00_rt );
          get2dDiff_( q_t0_l, q_i_t1_l, q_i_t0_r, q_i_t1_r,   q_t0_l, q_j_t1_l, q_j_t0_r, q_j_t1_r,  f11, f11_t, f11_rt );

#ifndef __do_depth_only__
      // some Hack combine both motions 
      P4 f00_m ( f00_t[0], f00_t[1], f00_rt[0], f00_rt[1] );
      P4 f11_m ( f11_t[0], f11_t[1], f11_rt[0], f11_rt[1] );
#endif

      ////////////////////////////////////////////////////////////////////////////////
      // if it is not connected at timestep t+1 why should it be at timestep t
      // take the maximum of the depth energies at the different timesteps:
      // the depth energy at timestep t

     f00_it = sqrt (epsilon);
     if (idN_a != idN_b)
     {
       f00_it = sqrt ( ( (f00|f00) + (f11|f11) + (f00|f11) + gamma*(dN00|dN00) ) * Scalar(3.3) + epsilon);
     }


#ifndef __do_depth_only__

     Scalar f00_R = sqrt (epsilon);
     if (idN_a != idN_b)
     {
       f00_R = sqrt ( ((f00_m|f00_m) + (f11_m|f11_m) + (f00_m|f11_m) + gamma*(dN00r|dN00r)) * Scalar(3.3) + epsilon);
     }
#endif

// now weight to pixel:

#ifndef __do_depth_only__
          f00_it = (min(f00_it, depthJump) + rotWeight*min(f00_R, rotJump)) * w_it * edgeWeight;
#else
          f00_it = min(f00_it, depthJump)  * w_it * edgeWeight;//maybe get rid of that (edge weights)
#endif

          // test just add a little bit here, scale it from the outside
          f00_it += (idN_a != idN_b) * edgeWeight * w_it * edgeScale;// * depthJump *.1;


              Idk.push_back(id_pix1_global);
              Idl.push_back(id_pix2_global);
              F00.push_back(f00_it);
              if ( f00_it != f00_it || std::numeric_limits<Scalar>::infinity()==f00_it )
              {
                printf("Bug in smooth f00 .., id:%d, idB:%d, pixIds: %d, %d; d0:%.1f\n", idN_a, idN_b, id_pix1_global, id_pix2_global, sqrt(norm0|norm0));
              }
            score += max(max(max(f00_it, f01_it), f10_it), f11_it);
            }// for displace
          }// for y
        }// for x
    if ( score != score || std::numeric_limits<Scalar>::infinity()==score )
    {
      printf ("WHAT WHY %f \n", score);
      mexEvalString("drawnow");
      int breakhere = 0;
    }
    return score;
  }

#define _use_safe_version_to_test_against_
#ifdef _use_safe_version_to_test_against_
///////////////////// new version without local rotation stuff ////////////////////
//compute_score_2D_homs
/// PURE 2d: returns the worst score possible, input is the normal (id) of the binary one and the P3 id configuration of the current solution

template<class Scalar>
Scalar
EvalEnergyFullFramePixel<Scalar>::
compute_score_combiDepth_consistent(const int id, const P4i& box, const std::vector<int>& loc2glob,// const int localSize, 
    const std::vector<M3>& homs_r, const std::vector<M3>& homs_l, const std::vector<M3>& homs_rt, //)//, std::vector<int>& seg2Prop )
    std::vector<Scalar>& F00, std::vector<Scalar>& F01, std::vector<Scalar>& F10, std::vector<Scalar>& F11,
    std::vector<Scalar>& F0, std::vector<Scalar>& F1, std::vector<int>& Idk, std::vector<int>& Idl )

  {
    int aoix = box[0]; int aoiX = box[2]; int aoiy = box[1]; int aoiY = box[3];

    const int    dx8[8] = {-1,0,1, 0, 1, -1, -1,  1};
    const int    dy8[8] = { 0,1,0,-1, 1, -1,  1, -1};
    const Scalar we8[8] = {1.,1.,1.,1., 1./sqrt(2.), 1./sqrt(2.), 1./sqrt(2.), 1./sqrt(2.)};//weights

    const int neighIts = (halfEdgeXY == NULL) ? 4:8;
    const int* dx= dx8;
    const int* dy= dy8;
    const Scalar* we= we8;

    assert (normals != NULL);
    assert (tra != NULL);
    assert (rot != NULL);
    assert (homs_r.size() > 1); // this is to be called before setHoms
    assert((*normals).size() > id);
    assert((*rot).size() > id);
    assert(homs_r.size() > id); 
    assert(homs_l.size() > id); 

    Scalar score (0.);
    
    P3 normFix = (*normals)[id];
    Scalar dn_k =  1. / sqrt(normFix|normFix);
    M3 rotk = (*rot)[id];
    P3 nk   = normFix * dn_k;
    P3 nrk  = rotk * nk;

    for (int j=aoiy;j<aoiY;j++)
    {
      Scalar f00_it(0), f01_it(0), f10_it(0), f11_it(0);
      int py = j;
      for (int ii=aoix;ii<aoiX;ii++)
      {
        int id_pix1_local = loc2glob[(j-aoiy)+ (aoiY-aoiy)*(ii-aoix)];//local id

        if (id_pix1_local < 0) continue; // the other pixel have an id in the mrf tored here

        int px = ii;
        int id_pix1_global = py+height*px;

        int idN_a = segImg[ py+height*px ];
        P3 norm0 = (*normals)[idN_a];
        Scalar dn_i =  1. / sqrt(norm0|norm0);
        M3 roti = (*rot)[idN_a];

        P3 n0  = norm0 * dn_i;
        P3 nr0 = roti * n0;

        assert(px<width);
        assert(py<height);
        assert(px>=0);
        assert(py>=0);
        assert((*normals).size() > idN_a);
        assert((*rot).size() > idN_a);
        assert(homs_r.size() > idN_a); 
        assert(homs_l.size() > idN_a); 
        assert(0 <= idN_a); 

        // run over the 4/8 neighbours, check whether oob or inside
        for(int displace = 0;displace < neighIts; displace++)
        {
          int qx = px+dx[displace];
          int qy = py+dy[displace];

          // oob? - indeed pixel can have term with border pixel
          if (qx<0 || qx>=width || qy<0 || qy>=height) continue;

          Scalar w_it = we[displace];

          // can also be -1
          int id_pix2_global = qx*height+qy;

          int id_pix2_local;

          if ( qx >= aoix && qx < aoiX && qy >= aoiy && qy < aoiY )
            id_pix2_local = loc2glob[(qy-aoiy)+ (aoiY-aoiy)*(qx-aoix)];
          else
            id_pix2_local = -1;

          // oob go further (id_pix2_local=-1) else only 1 occurence
          if(id_pix2_local > id_pix1_local ) continue; // id_pix2_local == -1 still ok

          ////////////////////////////////////////////
          Scalar edgeWeight(1);

          if (displace<4)
          {
            // one of those 2 craps
            if(abs (id_pix2_global - id_pix1_global) != 1)
              edgeWeight = halfEdgeX[ min(qx,px)*height + qy ];
            else
              edgeWeight = halfEdgeY[ px*(height-1) + min(qy,py) ];
          }
          else
          {
            if (displace<6)
              edgeWeight = halfEdgeXY [ min(qx,px)*(height-1) + min(qy,py) ];
            else
              edgeWeight = halfEdgeiXY[ min(qx,px)*(height-1) + min(qy,py) ];
          }

          int idN_b = segImg[ id_pix2_global ];

          if ( id == idN_a && id == idN_b ) continue;
          assert(0 <= idN_b);
          //---------------------------------------------------------------------
          {

          P3 norm1 = (*normals)[idN_b];     
          Scalar dn_j =  1. / sqrt(norm1|norm1);

          M3 rotj = (*rot)[idN_b];
          P3 n1  = norm1 * dn_j;
          P3 nr1 = rotj * n1;

          P3 p2d ( (px+qx+dy[displace])/2.+1., (py+qy+dx[displace])/2.+1., 1. );
          P3 q2d ( (px+qx-dy[displace])/2.+1., (py+qy-dx[displace])/2.+1., 1. );

#ifdef _Debug_out_
          printf("edge: %f,%f, %f,%f \n", p2d[0], p2d[1],q2d[0], q2d[1] );
#endif
          // now add the normal constraint to get a linear prior
          P3 dN00 = n0 - n1;
          P3 dN01 = n0 - nk;
          P3 dN10 = nk - n1;
          P3 dN11(0,0,0);// = normFix*dn_k - normFix *dn_k;// that's dull

          P3 dN00r = nr0 - nr1;
          P3 dN01r = nr0 - nrk;
          P3 dN10r = nrk - nr1;
          P3 dN11r(0,0,0);

          P2 f00, f11;
          P2 g00, g11;
          P2 h00, h11;

          P2 f00_t, f00_rt, f11_t, f11_rt; 
          P2 g00_t, g00_rt, g11_t, g11_rt; 
          P2 h00_t, h00_rt, h11_t, h11_rt; 

          P2 p_i_t0_r;
          P2 q_i_t0_r;
          P2 p_i_t1_l, p_i_t1_r;
          P2 q_i_t1_l, q_i_t1_r;

          P2 p_j_t0_r;
          P2 q_j_t0_r;
          P2 p_j_t1_l, p_j_t1_r;
          P2 q_j_t1_l, q_j_t1_r;

          P2 p_k_t0_r;
          P2 q_k_t0_r;
          P2 p_k_t1_l, p_k_t1_r;
          P2 q_k_t1_l, q_k_t1_r;

          P2 p_t0_l (p2d[0], p2d[1]);
          P2 q_t0_l (q2d[0], q2d[1]);

          P3 tmp1 = homs_r[idN_a]  * p2d;p_i_t0_r = P2(tmp1[0], tmp1[1]) * (1./tmp1[2]);
          P3 tmp2 = homs_r[idN_a]  * q2d;q_i_t0_r = P2(tmp2[0], tmp2[1]) * (1./tmp2[2]);
          P3 tmp3 = homs_l[idN_a]  * p2d;p_i_t1_l = P2(tmp3[0], tmp3[1]) * (1./tmp3[2]);
          P3 tmp4 = homs_l[idN_a]  * q2d;q_i_t1_l = P2(tmp4[0], tmp4[1]) * (1./tmp4[2]);
          P3 tmp5 = homs_rt[idN_a] * p2d;p_i_t1_r = P2(tmp5[0], tmp5[1]) * (1./tmp5[2]);
          P3 tmp6 = homs_rt[idN_a] * q2d;q_i_t1_r = P2(tmp6[0], tmp6[1]) * (1./tmp6[2]);

          tmp1 = homs_r[idN_b]  * p2d;p_j_t0_r = P2(tmp1[0], tmp1[1]) * (1./tmp1[2]);
          tmp2 = homs_r[idN_b]  * q2d;q_j_t0_r = P2(tmp2[0], tmp2[1]) * (1./tmp2[2]);
          tmp3 = homs_l[idN_b]  * p2d;p_j_t1_l = P2(tmp3[0], tmp3[1]) * (1./tmp3[2]);
          tmp4 = homs_l[idN_b]  * q2d;q_j_t1_l = P2(tmp4[0], tmp4[1]) * (1./tmp4[2]);
          tmp5 = homs_rt[idN_b] * p2d;p_j_t1_r = P2(tmp5[0], tmp5[1]) * (1./tmp5[2]);
          tmp6 = homs_rt[idN_b] * q2d;q_j_t1_r = P2(tmp6[0], tmp6[1]) * (1./tmp6[2]);

          tmp1 = homs_r[id]  * p2d;p_k_t0_r = P2(tmp1[0], tmp1[1]) * (1./tmp1[2]);
          tmp2 = homs_r[id]  * q2d;q_k_t0_r = P2(tmp2[0], tmp2[1]) * (1./tmp2[2]);
          tmp3 = homs_l[id]  * p2d;p_k_t1_l = P2(tmp3[0], tmp3[1]) * (1./tmp3[2]);
          tmp4 = homs_l[id]  * q2d;q_k_t1_l = P2(tmp4[0], tmp4[1]) * (1./tmp4[2]);
          tmp5 = homs_rt[id] * p2d;p_k_t1_r = P2(tmp5[0], tmp5[1]) * (1./tmp5[2]);
          tmp6 = homs_rt[id] * q2d;q_k_t1_r = P2(tmp6[0], tmp6[1]) * (1./tmp6[2]);
          ////////////////////////////  

          get2dDiff_( p_t0_l, p_i_t1_l, p_i_t0_r, p_i_t1_r,   p_t0_l, p_j_t1_l, p_j_t0_r, p_j_t1_r,  f00, f00_t, f00_rt );
          get2dDiff_( p_t0_l, p_i_t1_l, p_i_t0_r, p_i_t1_r,   p_t0_l, p_k_t1_l, p_k_t0_r, p_k_t1_r,  g00, g00_t, g00_rt );

          get2dDiff_( p_t0_l, p_k_t1_l, p_k_t0_r, p_k_t1_r,   p_t0_l, p_j_t1_l, p_j_t0_r, p_j_t1_r,  h00, h00_t, h00_rt );

          get2dDiff_( q_t0_l, q_i_t1_l, q_i_t0_r, q_i_t1_r,   q_t0_l, q_j_t1_l, q_j_t0_r, q_j_t1_r,  f11, f11_t, f11_rt );
          get2dDiff_( q_t0_l, q_i_t1_l, q_i_t0_r, q_i_t1_r,   q_t0_l, q_k_t1_l, q_k_t0_r, q_k_t1_r,  g11, g11_t, g11_rt );

          get2dDiff_( q_t0_l, q_k_t1_l, q_k_t0_r, q_k_t1_r,   q_t0_l, q_j_t1_l, q_j_t0_r, q_j_t1_r,  h11, h11_t, h11_rt );

#ifndef __do_depth_only__
      P4 f00_m ( f00_t[0], f00_t[1], f00_rt[0], f00_rt[1] );
      P4 f11_m ( f11_t[0], f11_t[1], f11_rt[0], f11_rt[1] );
      P4 g00_m ( g00_t[0], g00_t[1], g00_rt[0], g00_rt[1] );
      P4 g11_m ( g11_t[0], g11_t[1], g11_rt[0], g11_rt[1] );
      P4 h00_m ( h00_t[0], h00_t[1], h00_rt[0], h00_rt[1] );
      P4 h11_m ( h11_t[0], h11_t[1], h11_rt[0], h11_rt[1] );
#endif
      ////////////////////////////////////////////////////////////////////////////////
      // if it is not connected at timestep t+1 why should it be at timestep t
      // take the maximum of the depth energies at the different timesteps:
      // the depth energy at timestep t

     f01_it = sqrt ( ( (g00|g00) + (g11|g11) + (g00|g11) + gamma*(dN01|dN01) ) * Scalar(3.3) + epsilon);
     f11_it = sqrt (epsilon);
     if (idN_a != idN_b)
     {
       f00_it = sqrt ( ( (f00|f00) + (f11|f11) + (f00|f11) + gamma*(dN00|dN00) ) * Scalar(3.3) + epsilon);
       f10_it = sqrt ( ( (h00|h00) + (h11|h11) + (h00|h11) + gamma*(dN10|dN10) ) * Scalar(3.3) + epsilon);
     }
     else
     {
       f00_it = f11_it; 
       f10_it = f01_it;
     }

#ifndef __do_depth_only__

     Scalar f01_R = sqrt ( ((g00_m|g00_m) + (g11_m|g11_m) + (g00_m|g11_m) + gamma*(dN01r|dN01r)) * Scalar(3.3) + epsilon);
     Scalar f11_R = sqrt (epsilon);
     Scalar f00_R;
     Scalar f10_R;
     if (idN_a != idN_b)
     {
       f00_R = sqrt ( ((f00_m|f00_m) + (f11_m|f11_m) + (f00_m|f11_m) + gamma*(dN00r|dN00r)) * Scalar(3.3) + epsilon);
       f10_R = sqrt ( ((h00_m|h00_m) + (h11_m|h11_m) + (h00_m|h11_m) + gamma*(dN10r|dN10r)) * Scalar(3.3) + epsilon);
     }
     else
     {
       f00_R = f11_R; 
       f10_R = f01_R;
     }
#endif

// now weight to pixel:

#ifndef __do_depth_only__
          f00_it = (min(f00_it, depthJump) + rotWeight*min(f00_R, rotJump)) * w_it * edgeWeight;
          f10_it = (min(f10_it, depthJump) + rotWeight*min(f10_R, rotJump)) * w_it * edgeWeight;
          f01_it = (min(f01_it, depthJump) + rotWeight*min(f01_R, rotJump)) * w_it * edgeWeight;
          f11_it = (min(f11_it, depthJump) + rotWeight*min(f11_R, rotJump)) * w_it * edgeWeight;
#else
          f00_it = min(f00_it, depthJump)  * w_it * edgeWeight;//maybe get rid of that (edge weights)
          f10_it = min(f10_it, depthJump)  * w_it * edgeWeight;
          f01_it = min(f01_it, depthJump)  * w_it * edgeWeight;
          f11_it = min(f11_it, depthJump)  * w_it * edgeWeight;
#endif
          }

          // test just add a little bit here, scale it from the outside
          f00_it += (idN_a != idN_b) * edgeWeight * w_it * edgeScale;
          f10_it += (id    != idN_b) * edgeWeight * w_it * edgeScale;
          f01_it += (idN_a != id   ) * edgeWeight * w_it * edgeScale;

          {
            if (id_pix2_local < 0)
            {
              if ( f00_it != f00_it || std::numeric_limits<Scalar>::infinity()==f00_it || f10_it != f10_it || std::numeric_limits<Scalar>::infinity()==f10_it  )
              {
                printf("Bug in smooth f0 and f1, id:%d, id0 %d, id1 %d\n", id, idN_a, id);
              }

              if (id_pix1_local > F0.size() || id_pix1_local > F1.size())
                printf("Bug in smooth f0 and f1, idpixLocal:%d, sizes: %d,%d\n", id_pix1_local, F0.size(), F1.size() );

              F0[id_pix1_local] += f00_it;
              F1[id_pix1_local] += f10_it;
            }
            else
            {
              Idk.push_back(id_pix1_local);
              Idl.push_back(id_pix2_local);
              F00.push_back(f00_it);
              F10.push_back(f10_it);
              F01.push_back(f01_it);
              F11.push_back(f11_it);
              if ( f00_it != f00_it || std::numeric_limits<Scalar>::infinity()==f00_it || f10_it != f10_it || std::numeric_limits<Scalar>::infinity()==f10_it || 
                   f01_it != f01_it || std::numeric_limits<Scalar>::infinity()==f01_it || f11_it != f11_it || std::numeric_limits<Scalar>::infinity()==f11_it )
              {
                printf("Bug in smooth f00 .., id:%d, idB:%d box %d,%d,%d,%d, pixIds: %d, %d; df:%.1f, d0:%.1f\n", id, idN_b, box[0], box[1], box[2], box[3], id_pix1_local, id_pix2_local, sqrt(normFix|normFix), sqrt(norm0|norm0));//, sqrt(norm1|norm1) d1:%.1f, 
              }

            }
            score += max(max(max(f00_it, f01_it), f10_it), f11_it);
          }
        }
      }
    }
    if ( score != score || std::numeric_limits<Scalar>::infinity()==score )
    {
      printf ("WHAT WHY %f \n", score);
      mexEvalString("drawnow");
      int breakhere = 0;
    }
    return score;
  }
#endif


/// for a general camera setup: p and q are the 3d points at timstep t0 (of the same 2d point)
template<class Scalar>
inline void
EvalEnergyFullFramePixel<Scalar>::
get2dDiff( const P3& p_t0, const P3& q_t0, const P3& p_t1, const P3& q_t1, P2& dispDifft0, P2& dispDifft1, P2& flowDiff ) const 
  {
    P3 pixP_t0_l = Pl * p_t0;pixP_t0_l /= pixP_t0_l[2];
    P3 pixQ_t0_l = Pl * q_t0;pixQ_t0_l /= pixQ_t0_l[2];

    P3 pixP_t1_l = Pl * p_t1;pixP_t1_l /= pixP_t1_l[2];
    P3 pixQ_t1_l = Pl * q_t1;pixQ_t1_l /= pixQ_t1_l[2];

    P3 pixP_t0_r = Pr * p_t0 +pr;pixP_t0_r /= pixP_t0_r[2];
    P3 pixQ_t0_r = Pr * q_t0 +pr;pixQ_t0_r /= pixQ_t0_r[2];

    P3 pixP_t1_r = Pr * p_t1 +pr; pixP_t1_r /= pixP_t1_r[2];
    P3 pixQ_t1_r = Pr * q_t1 +pr; pixQ_t1_r /= pixQ_t1_r[2];

    P3 difft0 = (pixP_t0_r-pixP_t0_l) - (pixQ_t0_r-pixQ_t0_l);
    P3 difft1 = (pixP_t1_r-pixP_t1_l) - (pixQ_t1_r-pixQ_t1_l);
    P3 diffF  = (pixP_t1_l-pixP_t0_l) - (pixQ_t1_l-pixQ_t0_l);
    ///////////////////

    flowDiff   = P2(diffF[0],  diffF[1]);
    dispDifft0 = P2(difft0[0], difft0[1]);
    dispDifft1 = P2(difft1[0]-difft0[0], difft1[1]-difft0[1]); // correct but non sense different difference to cam so also different disp difference or not ?
  }

template<class Scalar>
inline void
EvalEnergyFullFramePixel<Scalar>::
get2dDiff_x2disp( const P2& pixP_t0_l, const P2& pixP_t1_l, const P2& pixP_t0_r, const P2& pixP_t1_r, 
                  const P2& pixQ_t0_l, const P2& pixQ_t1_l, const P2& pixQ_t0_r, const P2& pixQ_t1_r,
                  P4& dispDiff, P2& flowDiff ) const
  {
    P2 difft0   = (pixP_t0_r-pixP_t0_l) - (pixQ_t0_r-pixQ_t0_l);
    P2 difft1   = (pixP_t1_r-pixP_t1_l) - (pixQ_t1_r-pixQ_t1_l);
       flowDiff = (pixP_t1_l-pixP_t0_l) - (pixQ_t1_l-pixQ_t0_l);
    ///////////////////

    dispDiff[0] = difft0[0];
    dispDiff[1] = difft0[1];
    dispDiff[2] = difft1[0];
    dispDiff[3] = difft1[1];
  }

template<class Scalar>
inline void
EvalEnergyFullFramePixel<Scalar>::
get2dDiff_x1disp( const P2& pixP_t0_l, const P2& pixP_t1_l, const P2& pixP_t0_r, const P2& pixP_t1_r, 
                  const P2& pixQ_t0_l, const P2& pixQ_t1_l, const P2& pixQ_t0_r, const P2& pixQ_t1_r,
                  P4& dispDiff, P2& flowDiff ) const
  {
    P2 difft0   = (pixP_t0_r-pixP_t0_l) - (pixQ_t0_r-pixQ_t0_l);
    flowDiff = (pixP_t1_l-pixP_t0_l) - (pixQ_t1_l-pixQ_t0_l);
    ///////////////////

    dispDiff[0] = difft0[0];
    dispDiff[1] = difft0[1];
#ifdef    _dispDiffSmoothing_
    P2 difft1   = (pixP_t1_r-pixP_t1_l) - (pixQ_t1_r-pixQ_t1_l);
    dispDiff[2] = difft1[0];
    dispDiff[3] = difft1[1];
#endif
  };

  /// for a general camera setup: p and q are the 3d points at timstep t0 (of the same 2d point)
template<class Scalar>
inline void
EvalEnergyFullFramePixel<Scalar>::
get2dProjections( const P3& p_t0, P2& pixP_l, P2& pixP_r) const 
{
    P3 pixP_l_3 = Pl * p_t0;     pixP_l = P2(pixP_l_3[0], pixP_l_3[1]); pixP_l  /= pixP_l_3[2];
    P3 pixP_r_3 = Pr * p_t0 +pr; pixP_r = P2(pixP_r_3[0], pixP_r_3[1]); pixP_r  /= pixP_r_3[2];
};

template<class Scalar>
inline void
EvalEnergyFullFramePixel<Scalar>::
get2dDiff_( const P2& pixP_t0_l, const P2& pixP_t1_l, const P2& pixP_t0_r, const P2& pixP_t1_r, 
            const P2& pixQ_t0_l, const P2& pixQ_t1_l, const P2& pixQ_t0_r, const P2& pixQ_t1_r,
            P2& dispDifft0, P2& dispDifft1, P2& flowDiff ) const 
{
    P2 difft0   = (pixP_t0_r-pixP_t0_l) - (pixQ_t0_r-pixQ_t0_l);
    P2 difft1   = (pixP_t1_r-pixP_t1_l) - (pixQ_t1_r-pixQ_t1_l);
       flowDiff = (pixP_t1_l-pixP_t0_l) - (pixQ_t1_l-pixQ_t0_l);
    ///////////////////

//    flowDiff   = P2(diffF[0],  diffF[1]);
    dispDifft0 = P2(difft0[0], difft0[1]);
    dispDifft1 = P2(difft1[0]-difft0[0], difft1[1]-difft0[1]); // correct but non sense different difference to cam so also different disp difference or not ?
};

#undef __do_depth_only__

#endif