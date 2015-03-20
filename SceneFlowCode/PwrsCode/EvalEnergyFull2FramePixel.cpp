/*
Copyright (C) 2013 Christoph Vogel, PhD. Student ETH Zurich
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

///////////////////////////////////////////////////////
////////// Smoothness term per pixel        ///////////
///////////////////////////////////////////////////////
#ifndef __ENERGY_ROTTRANORPIXEL__CPP
#define __ENERGY_ROTTRANORPIXEL__CPP

#include "EvalEnergyFull2FramePixel.h"

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"
#include "DataDefinitions.h"
#ifndef _NO_OPENMP
#include <omp.h>
#endif

#include "mex.h" // include after VectorT and Mat3x3T

#include<limits.h>

#include <map>
#include <vector>

using namespace std;
using namespace Math;

//////////////////////////////////////////////////////////////////////////////////////

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
  }

/// PURE 3d: returns the worst score possible, input is the normal (id) of the binary one and the P3 id configuration of the current solution
#ifdef __USE3D__
template<class Scalar>
Scalar
EvalEnergyFullFramePixel<Scalar>::
compute_score_combiDepth( int id, int aoix, int aoiX, int aoiy, int aoiY, std::vector<int>& loc2glob, int localSize )
#else
template<class Scalar>
Scalar
EvalEnergyFullFramePixel<Scalar>::
compute_score_3D( int id, int aoix, int aoiX, int aoiy, int aoiY, std::vector<int>& loc2glob, int localSize )
#endif
  {
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

    Scalar score (0.);

    P3 normFix = (*normals)[id];
    Scalar dn_k =  1. / sqrt(normFix|normFix);

    M3 rotk = (*rot)[id];
    P3 trak = (*tra)[id];

    F0.clear();
    F1.clear();
    F0.resize(localSize, 0);
    F1.resize(localSize, 0);
    F00.clear();
    F01.clear();
    F10.clear();
    F11.clear();
    Idk.clear();
    Idl.clear();


#pragma omp parallel for schedule (static)
    for (int j=aoiy;j<aoiY;j++)
    {
      Scalar f00_it, f01_it, f10_it, f11_it;
      int py = j;
      int pos = j*width;int i = pos+aoix;
      for (int ii=aoix;ii<aoiX;ii++,i++)
      {
        int id_pix1_local = loc2glob[ii-aoix + (aoiX-aoix)*(j-aoiy)];

        int id_pix1_global = i;

        if (id_pix1_local < 0) continue; // the other pixel have an id in the mrf tored here

        int px = ii;
        int idN_a = segImg[ i ];
        P3 norm0 = (*normals)[idN_a];
        Scalar dn_i =  1. / sqrt(norm0|norm0);
        M3 roti = (*rot)[idN_a];
        P3 trai = (*tra)[idN_a];
#ifdef _use_patchCenters_
        P3 centerA ( &centers[3*idN_a] );
        P3 a_i = centerA * (1. / (norm0 | centerA));
        P3 a_k = centerA * (1. / (normFix | centerA));
#else
        P3 a_i(0.,0.,0.);
        P3 a_k(0.,0.,0.);
#endif
        // run over the 4/8 neighbours, check whether oob or inside
        for(int displace = 0;displace < neighIts; displace++)
        {
          int qx = px+dx[displace];
          int qy = py+dy[displace];

          // oob? - indeed pixel can have term with border pixel
          if (qx<0 || qx>=width || qy<0 || qy>=height) continue;

          Scalar w_it = we[displace];

          // can also be -1
          int id_pix2_global = qx+width*qy;

          int id_pix2_local;

          if ( qx >= aoix && qx < aoiX && qy >= aoiy && qy < aoiY )
            id_pix2_local = loc2glob[(qx-aoix)+ (aoiX-aoix)*(qy-aoiy)];
          else
            id_pix2_local = -1;

          // oob go further (id_pix2_local=-1) else only 1 occurence
          if(id_pix2_local > id_pix1_local ) continue; // id_pix2_local == -1 still ok

          Scalar edgeWeight(0);
          if (displace<4)
          {
            if(abs (id_pix2_global - id_pix1_global) == 1)
              edgeWeight = halfEdgeX[ min(qx,px) + (width-1)*qy ];
            else
              edgeWeight = halfEdgeY[ px + width*min(qy,py) ];
          }
          else
          {
            if (displace<6)
              edgeWeight = halfEdgeXY [ min(qx,px) + (width-1)*min(qy,py) ];
            else
              edgeWeight = halfEdgeiXY[ min(qx,px) + (width-1)*min(qy,py) ];
          }

          int idN_b = segImg[ qx+width*qy ];
//          if (displace < 4)
          {

          P3 norm1 = (*normals)[idN_b];
          Scalar dn_j =  1. / sqrt(norm1|norm1);

          M3 rotj = (*rot)[idN_b];
          P3 traj = (*tra)[idN_b];

          P3 p2d = iK * P3( (px+qx+dy[displace])/2., (py+qy+dx[displace])/2., 1. );
          P3 q2d = iK * P3( (px+qx-dy[displace])/2., (py+qy-dx[displace])/2., 1. );

#ifndef __do_depth_only__
#ifdef _use_patchCenters_
          P3 centerB ( &centers[3*idN_b] );
          P3 b_j = centerB * (1. / (norm1 | centerB));
          P3 b_k = centerB * (1. / (normFix | centerB));
#else
          P3 b_j = (0.,0.,0.);
          P3 b_k = (0.,0.,0.);
#endif
#endif
          Scalar dp_k = 1. / std::max(normFix | p2d, 0.00000001 );
          Scalar dq_k = 1. / std::max(normFix | q2d, 0.00000001 );
          Scalar dp_i = 1. / std::max(norm0   | p2d, 0.00000001 );
          Scalar dq_i = 1. / std::max(norm0   | q2d, 0.00000001 );
          Scalar dp_j = 1. / std::max(norm1   | p2d, 0.00000001 );
          Scalar dq_j = 1. / std::max(norm1   | q2d, 0.00000001 );

          if ( dp_k != dp_k ||  dq_k != dq_k ||  dp_i != dp_i  ||  dp_j != dp_j  ||  dq_i != dq_i  ||  dq_j != dq_j  )
          {
            int breakhere = 0;
            printf ("Normal all 0 %f \n", dp_i);
            mexEvalString("drawnow");
          }

          P3 p_i = p2d * dp_i;
          P3 q_j = q2d * dq_j;
          P3 p_j = p2d * dp_j;
          P3 q_i = q2d * dq_i;
          P3 p_ki = p2d * dp_k;
          P3 q_ki = q2d * dq_k;
          P3 p_kj = p2d * dp_k;
          P3 q_kj = q2d * dq_k;
          P3 p_k = p2d * dp_k;
          P3 q_k = q2d * dq_k;

          // now add the normal constraint to get a linear prior
          P3 dN00 = norm0*dn_i   - norm1   *dn_j;
          P3 dN01 = norm0*dn_i   - normFix *dn_k;
          P3 dN10 = normFix*dn_k - norm1   *dn_j;
          P3 dN11(0,0,0);// = normFix*dn_k - normFix *dn_k;// that's dull

          P3 f00 = p_i - p_j;
          P3 f11 = q_i - q_j;

          P3 g00 = p_i - p_kj;
          P3 g11 = q_i - q_kj;

          P3 h00 = p_ki - p_j;
          P3 h11 = q_ki - q_j;

          P3 e00 = p_ki - p_kj;
          P3 e11 = q_ki - q_kj;

          f00_it = sqrt ( ( (f00|f00) + (f11|f11) + (f00|f11) + gamma*(dN00|dN00) ) * Scalar(3.3) + epsilon);
          f10_it = sqrt ( ( (h00|h00) + (h11|h11) + (h00|h11) + gamma*(dN10|dN10) ) * Scalar(3.3) + epsilon);
          f01_it = sqrt ( ( (g00|g00) + (g11|g11) + (g00|g11) + gamma*(dN01|dN01) ) * Scalar(3.3) + epsilon);
          f11_it = sqrt (epsilon);

#ifndef __do_depth_only__
          P3 pi_t1 = roti*(p_i-a_i)+trai+a_i;
          P3 pj_t1 = rotj*(p_j-b_j)+traj+b_j;

          P3 qi_t1 = roti*(q_i-a_i)+trai+a_i;
          P3 qj_t1 = rotj*(q_j-b_j)+traj+b_j;

          P3 pkA_t1 = rotk*(p_k-a_k)+trak+a_k;
          P3 qkA_t1 = rotk*(q_k-a_k)+trak+a_k;

          P3 pkB_t1 = rotk*(p_k-b_k)+trak+b_k;
          P3 qkB_t1 = rotk*(q_k-b_k)+trak+b_k;

          P3 dN00r = roti * norm0*dn_i   - rotj * norm1*dn_j;
          P3 dN01r = roti * norm0*dn_i   - rotk * normFix*dn_k;
          P3 dN10r = rotk * normFix*dn_k - rotj * norm1*dn_j;
          P3 dN11r(0,0,0);// = rotk * normFix*dn_k - rotk * normFix*dn_k;

          // in the paper is the 3d motion difference, so:  no penalty at global motion
          P3 f00_rt = pi_t1 - pj_t1 - f00;
          P3 f11_rt = qi_t1 - qj_t1 - f11;

          P3 g00_rt = pi_t1 - pkB_t1 - g00;
          P3 g11_rt = qi_t1 - qkB_t1 - g11;

          P3 h00_rt = pkA_t1 - pj_t1 - h00;
          P3 h11_rt = qkA_t1 - qj_t1 - h11;

          P3 e00_rt = pkA_t1 - pkB_t1 - e00;
          P3 e11_rt = qkA_t1 - qkB_t1 - e11;

      ////////////////////////////////////////////////////////////////////////////////
      // if it is not connected at timestep t+1 why should it be at timestep t
      // take the maximum of the depth energies at the different timesteps:
      // the depth energy at timestep t

      Scalar f00_R = sqrt ( ((f00_rt|f00_rt) + (f11_rt|f11_rt) + (f00_rt|f11_rt) + gamma*(dN00r|dN00r)) * Scalar(3.3) + epsilon);
      Scalar f10_R = sqrt ( ((h00_rt|h00_rt) + (h11_rt|h11_rt) + (h00_rt|h11_rt) + gamma*(dN10r|dN10r)) * Scalar(3.3) + epsilon);
      Scalar f01_R = sqrt ( ((g00_rt|g00_rt) + (g11_rt|g11_rt) + (g00_rt|g11_rt) + gamma*(dN01r|dN01r)) * Scalar(3.3) + epsilon);
      Scalar f11_R = sqrt ( ((e00_rt|e00_rt) + (e11_rt|e11_rt) + (e00_rt|e11_rt) + gamma*(dN11r|dN11r)) * Scalar(3.3) + epsilon);
#endif

// now weight to pixel:

#ifndef __do_depth_only__

          f00_it = (min(f00_it, depthJump) + rotWeight*min(f00_R, rotJump)) * w_it;// * edgeWeight;
          f10_it = (min(f10_it, depthJump) + rotWeight*min(f10_R, rotJump)) * w_it;// * edgeWeight;
          f01_it = (min(f01_it, depthJump) + rotWeight*min(f01_R, rotJump)) * w_it;// * edgeWeight;
          f11_it = (min(f11_it, depthJump) + rotWeight*min(f11_R, rotJump)) * w_it;// * edgeWeight;
#else
          f00_it = min(f00_it, depthJump)  * w_it * edgeWeight;
          f10_it = min(f10_it, depthJump)  * w_it * edgeWeight;
          f01_it = min(f01_it, depthJump)  * w_it * edgeWeight;
          f11_it = min(f11_it, depthJump)  * w_it * edgeWeight;
#endif
          }

          // test just add a little bit here, scale it from the outside
          f00_it += (idN_a != idN_b) * edgeWeight * w_it;
          f10_it += (id    != idN_b) * edgeWeight * w_it;
          f01_it += (idN_a != id   ) * edgeWeight * w_it;

#pragma omp critical
          {
            if (id_pix2_local < 0)
            {
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
            }

            score += max(max(max(f00_it, f01_it), f10_it), f11_it);

            if ( score != score || std::numeric_limits<Scalar>::infinity()==score )
            {
              printf ("Normal 0 vector %f \n", score);
              mexEvalString("drawnow");
              int breakhere = 0;
            }
          }

        }
      }
    }
    return score;
  }

  //compute_score_2D_homs
  /// PURE 2d: returns the worst score possible, input is the normal (id) of the binary one and the P3 id configuration of the current solution
#ifndef __USE3D__
template<class Scalar>
Scalar
EvalEnergyFullFramePixel<Scalar>::
compute_score_combiDepth(const int id, int aoix, int aoiX, int aoiy, int aoiY, std::vector<int>& loc2glob, int localSize )
#else
template<class Scalar>
Scalar
EvalEnergyFullFramePixel<Scalar>::
compute_score_2D(const int id, int aoix, int aoiX, int aoiy, int aoiY, std::vector<int>& loc2glob, int localSize )
#endif
  {
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
    
    P3 normFix = (*normals)[id];
    Scalar dn_k =  1. / sqrt(normFix|normFix);

    M3 rotk = (*rot)[id];
    P3 nk   = normFix * (1./dn_k);
    P3 nrk  = rotk * nk;
    
    F0.clear();
    F1.clear();
    F0.resize(localSize, 0);
    F1.resize(localSize, 0);
    F00.clear();
    F01.clear();
    F10.clear();
    F11.clear();
    Idk.clear();
    Idl.clear();

#pragma omp parallel for schedule (static)
    for (int j=aoiy;j<aoiY;j++)
    {
      Scalar f00_it, f01_it, f10_it, f11_it;
      int py = j;
      int pos = j*width;int i = pos+aoix;
      for (int ii=aoix;ii<aoiX;ii++,i++)
      {
        int id_pix1_local = loc2glob[ii-aoix + (aoiX-aoix)*(j-aoiy)];

        int id_pix1_global = i;

        if (id_pix1_local < 0) continue; // the other pixel have an id in the mrf tored here

        int px = ii;
        int idN_a = segImg[ i ];
        P3 norm0 = (*normals)[idN_a];
        Scalar dn_i =  1. / sqrt(norm0|norm0);
        M3 roti = (*rot)[idN_a];

        P3 n0  = norm0 * (1./dn_i);
        P3 nr0 = roti * n0;

        // run over the 4/8 neighbours, check whether oob or inside
        for(int displace = 0;displace < neighIts; displace++)
        {
          int qx = px+dx[displace];
          int qy = py+dy[displace];

          // oob? - indeed pixel can have term with border pixel
          if (qx<0 || qx>=width || qy<0 || qy>=height) continue;

          Scalar w_it = we[displace];

          // can also be -1
          int id_pix2_global = qx+width*qy;

          int id_pix2_local;

          if ( qx >= aoix && qx < aoiX && qy >= aoiy && qy < aoiY )
            id_pix2_local = loc2glob[(qx-aoix)+ (aoiX-aoix)*(qy-aoiy)];
          else
            id_pix2_local = -1;

          // oob go further (id_pix2_local=-1) else only 1 occurence
          if(id_pix2_local > id_pix1_local ) continue; // id_pix2_local == -1 still ok

          Scalar edgeWeight(0);
          if (displace<4)
          {
            if(abs (id_pix2_global - id_pix1_global) == 1)
              edgeWeight = halfEdgeX[ min(qx,px) + (width-1)*qy ];
            else
              edgeWeight = halfEdgeY[ px + width*min(qy,py) ];
          }
          else
          {
            if (displace<6)
              edgeWeight = halfEdgeXY [ min(qx,px) + (width-1)*min(qy,py) ];
            else
              edgeWeight = halfEdgeiXY[ min(qx,px) + (width-1)*min(qy,py) ];
          }

          int idN_b = segImg[ qx+width*qy ];
//          if (displace < 4)
          {

          P3 norm1 = (*normals)[idN_b];     
          Scalar dn_j =  1. / sqrt(norm1|norm1);

          M3 rotj = (*rot)[idN_b];
          P3 n1  = norm1 * (1./dn_j);
          P3 nr1 = rotj * n1;

          P3 p2d = P3( (px+qx+dy[displace])/2., (py+qy+dx[displace])/2., 1. );
          P3 q2d = P3( (px+qx-dy[displace])/2., (py+qy-dx[displace])/2., 1. );

          // now add the normal constraint to get a linear prior
          P3 dN00 = n0 - n1;
          P3 dN01 = n0 - nk;
          P3 dN10 = nk - n1;
          P3 dN11(0,0,0);// = normFix*dn_k - normFix *dn_k;// that's dull

          P3 dN00r = nr0 - nr1;
          P3 dN01r = nr0 - nrk;
          P3 dN10r = nrk - nr1;
          P3 dN11r(0,0,0);

         P2 f00(0.,0.), f11(0.,0.);
         P2 g00(0.,0.), g11(0.,0.);
         P2 h00(0.,0.), h11(0.,0.);

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

          P2 p_t0_l, q_t0_l;

          P3 tmp1 = p2d;                p_t0_l = P2(tmp1[0], tmp1[1]) * (1./tmp1[2]);
          P3 tmp2 = q2d;                q_t0_l = P2(tmp2[0], tmp2[1]) * (1./tmp2[2]);

          tmp1 = homs_r[idN_a]  * p2d;p_i_t0_r = P2(tmp1[0], tmp1[1]) * (1./tmp1[2]);
          tmp2 = homs_r[idN_a]  * q2d;q_i_t0_r = P2(tmp2[0], tmp2[1]) * (1./tmp2[2]);
          tmp1 = homs_l[idN_a]  * p2d;p_i_t1_l = P2(tmp1[0], tmp1[1]) * (1./tmp1[2]);
          tmp2 = homs_l[idN_a]  * q2d;q_i_t1_l = P2(tmp2[0], tmp2[1]) * (1./tmp2[2]);
          tmp1 = homs_rt[idN_a] * p2d;p_i_t1_r = P2(tmp1[0], tmp1[1]) * (1./tmp1[2]);
          tmp2 = homs_rt[idN_a] * q2d;q_i_t1_r = P2(tmp2[0], tmp2[1]) * (1./tmp2[2]);

          tmp1 = homs_r[idN_b]  * p2d;p_j_t0_r = P2(tmp1[0], tmp1[1]) * (1./tmp1[2]);
          tmp2 = homs_r[idN_b]  * q2d;q_j_t0_r = P2(tmp2[0], tmp2[1]) * (1./tmp2[2]);
          tmp1 = homs_l[idN_b]  * p2d;p_j_t1_l = P2(tmp1[0], tmp1[1]) * (1./tmp1[2]);
          tmp2 = homs_l[idN_b]  * q2d;q_j_t1_l = P2(tmp2[0], tmp2[1]) * (1./tmp2[2]);
          tmp1 = homs_rt[idN_b] * p2d;p_j_t1_r = P2(tmp1[0], tmp1[1]) * (1./tmp1[2]);
          tmp2 = homs_rt[idN_b] * q2d;q_j_t1_r = P2(tmp2[0], tmp2[1]) * (1./tmp2[2]);

          tmp1 = homs_r[id]  * p2d;p_k_t0_r = P2(tmp1[0], tmp1[1]) * (1./tmp1[2]);
          tmp2 = homs_r[id]  * q2d;q_k_t0_r = P2(tmp2[0], tmp2[1]) * (1./tmp2[2]);
          tmp1 = homs_l[id]  * p2d;p_k_t1_l = P2(tmp1[0], tmp1[1]) * (1./tmp1[2]);
          tmp2 = homs_l[id]  * q2d;q_k_t1_l = P2(tmp2[0], tmp2[1]) * (1./tmp2[2]);
          tmp1 = homs_rt[id] * p2d;p_k_t1_r = P2(tmp1[0], tmp1[1]) * (1./tmp1[2]);
          tmp2 = homs_rt[id] * q2d;q_k_t1_r = P2(tmp2[0], tmp2[1]) * (1./tmp2[2]);
          ////////////////////////////  

          get2dDiff_Fast( p_t0_l, p_i_t1_l, p_i_t0_r, p_i_t1_r,   p_t0_l, p_j_t1_l, p_j_t0_r, p_j_t1_r,  f00, f00_t, f00_rt );
          get2dDiff_Fast( p_t0_l, p_i_t1_l, p_i_t0_r, p_i_t1_r,   p_t0_l, p_k_t1_l, p_k_t0_r, p_k_t1_r,  g00, g00_t, g00_rt );

          get2dDiff_Fast( p_t0_l, p_k_t1_l, p_k_t0_r, p_k_t1_r,   p_t0_l, p_j_t1_l, p_j_t0_r, p_j_t1_r,  h00, h00_t, h00_rt );

          get2dDiff_Fast( q_t0_l, q_i_t1_l, q_i_t0_r, q_i_t1_r,   q_t0_l, q_j_t1_l, q_j_t0_r, q_j_t1_r,  f11, f11_t, f11_rt );
          get2dDiff_Fast( q_t0_l, q_i_t1_l, q_i_t0_r, q_i_t1_r,   q_t0_l, q_k_t1_l, q_k_t0_r, q_k_t1_r,  g11, g11_t, g11_rt );

          get2dDiff_Fast( q_t0_l, q_k_t1_l, q_k_t0_r, q_k_t1_r,   q_t0_l, q_j_t1_l, q_j_t0_r, q_j_t1_r,  h11, h11_t, h11_rt );

#ifndef __do_depth_only__
      // combine both motions 
      P4 f00_m ( f00_t[0], f00_t[1], f00_rt[0], f00_rt[1] );
      P4 f11_m ( f11_t[0], f11_t[1], f11_rt[0], f11_rt[1] );
      P4 g00_m ( g00_t[0], g00_t[1], g00_rt[0], g00_rt[1] );
      P4 g11_m ( g11_t[0], g11_t[1], g11_rt[0], g11_rt[1] );
      P4 h00_m ( h00_t[0], h00_t[1], h00_rt[0], h00_rt[1] );
      P4 h11_m ( h11_t[0], h11_t[1], h11_rt[0], h11_rt[1] );
#endif

      ////////////////////////////////////////////////////////////////////////////////
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
      // the same:
      //      f00_it = sqrt ( ( (f00|f00) + (f11|f11) + (f00|f11) + gamma*(dN00|dN00) ) * Scalar(3.3) + epsilon);
      //      f10_it = sqrt ( ( (h00|h00) + (h11|h11) + (h00|h11) + gamma*(dN10|dN10) ) * Scalar(3.3) + epsilon);
      //      f01_it = sqrt ( ( (g00|g00) + (g11|g11) + (g00|g11) + gamma*(dN01|dN01) ) * Scalar(3.3) + epsilon);
      //      f11_it = sqrt (epsilon);

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
          f00_it = (min(f00_it, depthJump) + rotWeight*min(f00_R, rotJump)) * w_it;
          f10_it = (min(f10_it, depthJump) + rotWeight*min(f10_R, rotJump)) * w_it;
          f01_it = (min(f01_it, depthJump) + rotWeight*min(f01_R, rotJump)) * w_it;
          f11_it = (min(f11_it, depthJump) + rotWeight*min(f11_R, rotJump)) * w_it;
#else
          f00_it = min(f00_it, depthJump)  * w_it ;
          f10_it = min(f10_it, depthJump)  * w_it ;
          f01_it = min(f01_it, depthJump)  * w_it ;
          f11_it = min(f11_it, depthJump)  * w_it ;
#endif
          }
          // test just add a little bit here, scale it from the outside
          f00_it += (idN_a != idN_b) * edgeWeight * w_it ;
          f10_it += (id    != idN_b) * edgeWeight * w_it ;
          f01_it += (idN_a != id   ) * edgeWeight * w_it ;

#pragma omp critical
          {
            if (id_pix2_local < 0)
            {
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
            }

            score += max(max(max(f00_it, f01_it), f10_it), f11_it);

            if ( score != score || std::numeric_limits<Scalar>::infinity()==score )
            {
              printf ("Normal 0 vector %f \n", score);
              mexEvalString("drawnow");
              int breakhere = 0;
            }
          }

        }
      }
    }
    return score;
  }

#endif
