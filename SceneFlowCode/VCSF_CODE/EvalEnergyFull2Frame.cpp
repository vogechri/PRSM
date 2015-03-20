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

#ifndef __ENERGY_ROTTRANOR__cpp
#define __ENERGY_ROTTRANOR__cpp

///////////////////////////////////////////////////////
////////// First version of faster iterations///////////
////////////////////////////////////////////////////////

#include "DataDefinitionsVC.h"
#include "EvalEnergyFull2Frame.h"
#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#ifndef _NO_OPENMP
#include <omp.h>
#endif

#include "mex.h" // include after VectorT and Mat3x3T

#include <map>
#include <vector>

using namespace std;
using namespace Math;

/// segimg must be flipped
template<class Scalar>
void
EvalEnergyFullFrame<Scalar>::
buildWeights(int* segImg, int width, int height, int nSegments)
{
//    const int    dx4[4] = {-1,0,1,0};
//    const int    dy4[4] = {0,1,0,-1};
//    const Scalar we4[4] = {1.,1.,1.,1.};//weights

    edgeWeightMap.clear();
    edgeWeightMap.resize(nSegments);

    const int    dx8[8] = {-1,0,1, 0, 1, -1, -1,  1};
    const int    dy8[8] = { 0,1,0,-1, 1, -1,  1, -1};
    const Scalar we8[8] = {1.,1.,1.,1., 1./sqrt(2.), 1./sqrt(2.), 1./sqrt(2.), 1./sqrt(2.)};//weights

    const int neighIts = (halfEdgeXY == NULL) ? 4:8;
    const int* dx= dx8;
    const int* dy= dy8;
    const Scalar* we= we8;

    // paralellism shoud be outside anyway
//#pragma omp parallel for schedule (static)
    for (int j=0;j<height;j++)
    {
      int py=j;
      for (int i=0;i<width;i++)
      {
        int posA = j*width+i;

        int px=i;
        int idN_a = segImg[ posA ];
        // run over the 4/8 neighbours, check whether oob or inside
        for(int displace = 0;displace < neighIts; displace++)
        {
          int qx = px+dx[displace];
          int qy = py+dy[displace];

          // oob?
          if (qx<0 || qx>=width || qy<0 || qy>=height) continue;
          int posB = qy*width+qx;

          if( posB > posA ) continue; // since stored in both arrays - right ?
		  
          int idN_b = segImg[ posB ];

          if (idN_a == idN_b) continue;

          Scalar w_it = we[displace];

          // can also be -1
          int id_pix2_global = qx+width*qy;

          Scalar edgeWeight(0);
          if (displace<4)
          {
            if(abs (px - qx) == 1)
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

          /// now put into the map:
          Scalar theWeight = edgeWeight * w_it;
          
          typename std::map< int, Scalar>& segiMap = edgeWeightMap[idN_a];
          typename std::map< int, Scalar>::iterator s_it, s_end (segiMap.end());
          s_it= segiMap.find( idN_b );
         
          if (s_it != segiMap.end())
            s_it->second += theWeight;
          else
            segiMap.insert(typename std::pair<int,Scalar> (idN_b, theWeight) );

          std::map< int, Scalar>& segjMap = edgeWeightMap[idN_b];
          s_end = segjMap.end();
          s_it  = segjMap.find( idN_a );
         
          if (s_it != segjMap.end())
            s_it->second += theWeight;
          else
            segjMap.insert(typename std::pair<int,Scalar> (idN_a, theWeight) );
        }
      }
    }
  }

  /// segimg must be flipped, builds edgeweightMap, summing over shared edge taking image edges into account
template<class Scalar>
void EvalEnergyFullFrame<Scalar>::
buildWeightsVC(int* segImg, int width, int height, int nSegments)
  {
    edgeWeightMap.clear();
    edgeWeightMap.resize(nSegments);

    const int    dx8[8] = {-1,0,1, 0, 1, -1, -1,  1};
    const int    dy8[8] = { 0,1,0,-1, 1, -1,  1, -1};
    const Scalar we8[8] = {1.,1.,1.,1., 1./sqrt(2.), 1./sqrt(2.), 1./sqrt(2.), 1./sqrt(2.)};//weights

    const int neighIts = (halfEdgeXY == NULL) ? 4:8;
    const int* dx= dx8;
    const int* dy= dy8;
    const Scalar* we= we8;

    for (int j=0;j<height;j++)
    {
      int py=j;
      for (int i=0;i<width;i++)
      {
        int px=i;
        ////////// start new ///////////////
        int id_pix1_global = py+height*px;
        int idN_a          = segImg[ id_pix1_global ];

        for(int displace = 0;displace < neighIts; displace++)
        {
          int qx = px+dx[displace];
          int qy = py+dy[displace];

          // oob? - indeed pixel can have term with border pixel
          if (qx<0 || qx>=width || qy<0 || qy>=height) continue;

          Scalar w_it = we[displace];

          // can also be -1
          int id_pix2_global = qx*height+qy;//qx+width*qy;

          int idN_b = segImg[ id_pix2_global ];
          if (idN_a == idN_b) continue;

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
          /////////////// end new /////////////////////////


          /// now put into the map: summing over shared edge
          Scalar theWeight = edgeWeight * w_it;
          
          typename std::map< int, Scalar>& segiMap = edgeWeightMap[idN_a];
          typename std::map< int, Scalar>::iterator s_it, s_end (segiMap.end());
          s_it= segiMap.find( idN_b );
         
          if (s_it != segiMap.end())
            s_it->second += theWeight;
          else
            segiMap.insert(typename std::pair<int,Scalar> (idN_b, theWeight) );

          std::map< int, Scalar>& segjMap = edgeWeightMap[idN_b];
          s_end = segjMap.end();
          s_it  = segjMap.find( idN_a );
         
          if (s_it != segjMap.end())
            s_it->second += theWeight;
          else
            segjMap.insert(typename std::pair<int,Scalar> (idN_a, theWeight) );
        }
      }
    }
  }

template<class Scalar>
void
EvalEnergyFullFrame<Scalar>::
prepare( int nSegments_, const mxArray* edges_, const mxArray* weights_, Scalar* centers_ )
  {
    nSegments = nSegments_;
    weights.clear();
    edgesP.clear();
    edgesQ.clear();
    Idk.clear();
    Idl.clear();
    centersA.clear();
    centersB.clear();

    // reserves not enough but i push_back, so ...
    weights.reserve( nSegments );
    edgesP.reserve(  nSegments );
    edgesQ.reserve(  nSegments );
    Idk.reserve( nSegments );
    Idl.reserve( nSegments );
    centersA.reserve( nSegments );
    centersB.reserve( nSegments );

    for (int i = 0; i < nSegments ;i++)
    {
      Scalar* wts   = (Scalar*) mxGetPr( mxGetCell(weights_, i) );
      int nIds      = (int)     mxGetNumberOfElements( mxGetCell(weights_, i) );
      Scalar* edge  = (Scalar*) mxGetPr( mxGetCell(edges_, i) );

//      P3 centerA = P3( centers_[ 2*i ], centers_[ 2*i+1], 1.0 );
      P3 centerA = P3( &centers_[ 3*i ] );

      for (int j =0; j < nIds ;j++)
      {
        int id = (int) (edge [5*j ])-1;

        if (id < i) continue;// can t we even break ?

        Scalar w = wts[j];

        P3 centerB = P3( &centers_[ 3*id ] );

        P3 p2d ( edge [5*j+1 ], edge [5*j+2 ], 1. );
        P3 q2d ( edge [5*j+3 ], edge [5*j+4 ], 1. );

        edgesP.push_back(p2d);
        edgesQ.push_back(q2d);

        centersA.push_back(centerA);
        centersB.push_back(centerB);

        Idk.push_back(i);
        Idl.push_back(id);

        F11.push_back(sqrt(epsilon)*w);

        weights.push_back(w);
      }}

    F01.resize(F11.size());
    F00.resize(F11.size());
    F10.resize(F11.size());
  }

/// segimage flipped so also width and height
template<class Scalar>
void
EvalEnergyFullFrame<Scalar>::
prepareOwnWeights( int nSegments_, const mxArray* edges_, Scalar* centers_, int* segImg, int width, int height)
  {
    nSegments = nSegments_;
    weights.clear();
    edgesP.clear();
    edgesQ.clear();
    Idk.clear();
    Idl.clear();
    centersA.clear();
    centersB.clear();

#ifndef __vc_version__
    buildWeights( segImg, width, height,  nSegments);
#else
    buildWeightsVC( segImg, width, height,  nSegments);  
#endif

    // reserves not enough but i push_back, so ...
    weights.reserve( nSegments );
    edgesP.reserve(  nSegments );
    edgesQ.reserve(  nSegments );
    Idk.reserve( nSegments );
    Idl.reserve( nSegments );
    centersA.reserve( nSegments );
    centersB.reserve( nSegments );

    for (int i = 0; i < nSegments ;i++)
    {
      int nIds      = (int)     mxGetNumberOfElements( mxGetCell(edges_, i) )/5;
      Scalar* edge  = (Scalar*) mxGetPr( mxGetCell(edges_, i) );

      P3 centerA = P3( &centers_[ 3*i ] );

      for (int j =0; j < nIds ;j++)
      {
        int id = (int) (edge [5*j ])-1;

        if (id < i) continue;// can t we even break ?

        Scalar w = (edgeWeightMap[i])[id];

        P3 centerB = P3( &centers_[ 3*id ] );

        P3 p2d ( edge [5*j+1 ], edge [5*j+2 ], 1. );
        P3 q2d ( edge [5*j+3 ], edge [5*j+4 ], 1. );

        edgesP.push_back(p2d);
        edgesQ.push_back(q2d);

        centersA.push_back(centerA);
        centersB.push_back(centerB);

        Idk.push_back(i);
        Idl.push_back(id);

        F11.push_back(sqrt(epsilon)*w);

        weights.push_back(w);
      }}

    F01.resize(F11.size());
    F00.resize(F11.size());
    F10.resize(F11.size());
  }

#ifdef __USE3D__
template<class Scalar>
Scalar
EvalEnergyFullFrame<Scalar>::
compute_score_Fuse00( std::vector<int>& curr_ids )
{
compute_score_Fuse( curr_ids, curr_ids );
}
#endif
  
#ifdef __USE3D__
template<class Scalar>
Scalar
EvalEnergyFullFrame<Scalar>::
compute_score_Fuse( std::vector<int>& curr_ids, std::vector<int>& curr_ids2 )
#else
template<class Scalar>
Scalar
EvalEnergyFullFrame<Scalar>::
compute_score_Fuse3D( std::vector<int>& curr_ids, std::vector<int>& curr_ids2 )
#endif
  {
    if (normals == NULL) return 0;

    Scalar score (0.);

#pragma omp parallel for schedule (static)
    for (int www = 0; www < weights.size(); www++)
    {
      Scalar* w_it   = &(weights[www]);
      Scalar* f00_it = &(F00[www]);
      Scalar* f01_it = &(F01[www]);
      Scalar* f10_it = &(F10[www]);
      Scalar* f11_it = &(F11[www]);

      int* idl_it = &(Idl[www]);
      int* idk_it = &(Idk[www]);

      P3* eP_it = &(edgesP[www]);
      P3* eQ_it = &(edgesQ[www]);
      P3* cA_it = &(centersA[www]);
      P3* cB_it = &(centersB[www]);

      int idN_a = curr_ids[ *idk_it ];
      int idN_b = curr_ids[ *idl_it ];

      P3 norm0 = (*normals)[idN_a];
      P3 norm1 = (*normals)[idN_b];

      M3 roti = (*rot)[idN_a];
      M3 rotj = (*rot)[idN_b];
      P3 trai = (*tra)[idN_a];
      P3 traj = (*tra)[idN_b];

      int idK_a = curr_ids2[ *idk_it ];
      int idK_b = curr_ids2[ *idl_it ];
      P3 normFix_i = (*normals)[idK_a];
      P3 normFix_j = (*normals)[idK_b];

      Scalar dn_ki =  1. / sqrt(normFix_i|normFix_i);
      Scalar dn_kj =  1. / sqrt(normFix_j|normFix_j);

      M3 rotki = (*rot)[idK_a];
      M3 rotkj = (*rot)[idK_b];
      P3 traki = (*tra)[idK_a];
      P3 trakj = (*tra)[idK_b];

      P3 p2d ( *eP_it );
      P3 q2d ( *eQ_it );

      P3 centerA ( *cA_it );
      P3 centerB ( *cB_it );

      P3 a_i  = centerA * (1. / (norm0 | centerA));
      P3 b_j  = centerB * (1. / (norm1 | centerB));
      P3 a_k  = centerA * (1. / (normFix_i | centerA));
      P3 b_k  = centerB * (1. / (normFix_j | centerB));

      Scalar dp_ki = 1. / max(normFix_i | p2d, 0.00000001 );
      Scalar dq_ki = 1. / max(normFix_i | q2d, 0.00000001 );
      Scalar dp_kj = 1. / max(normFix_j | p2d, 0.00000001 );
      Scalar dq_kj = 1. / max(normFix_j | q2d, 0.00000001 );
      Scalar dp_i  = 1. / max(norm0     | p2d, 0.00000001 );
      Scalar dq_i  = 1. / max(norm0     | q2d, 0.00000001 );
      Scalar dp_j  = 1. / max(norm1     | p2d, 0.00000001 );
      Scalar dq_j  = 1. / max(norm1     | q2d, 0.00000001 );

      //  patch 'depth'
      Scalar dn_i =  1. / sqrt(norm0|norm0);
      Scalar dn_j =  1. / sqrt(norm1|norm1);

      P3 p_i  = p2d * dp_i;
      P3 q_j  = q2d * dq_j;
      P3 p_j  = p2d * dp_j;
      P3 q_i  = q2d * dq_i;
      P3 p_ki = p2d * dp_ki;
      P3 q_ki = q2d * dq_ki;
      P3 p_kj = p2d * dp_kj;
      P3 q_kj = q2d * dq_kj;


      ///////// 2nd part rotations:
      P3 p_i_0 = p_i - a_i;
      P3 q_i_0 = q_i - a_i;

      P3 p_j_0 = p_j - b_j;
      P3 q_j_0 = q_j - b_j;

      P3 p_k_0A = p_ki - a_k;
      P3 q_k_0A = q_ki - a_k;

      P3 p_k_0B = p_kj - b_k;
      P3 q_k_0B = q_kj - b_k;


      // now add the normal constraint to get a linear prior
      P3 dN00 = norm0*dn_i      - norm1     *dn_j;
      P3 dN01 = norm0*dn_i      - normFix_j *dn_kj;
      P3 dN10 = normFix_i*dn_ki - norm1     *dn_j;
      P3 dN11 = normFix_i*dn_ki - normFix_j *dn_kj;

      P3 dN00r = roti  * norm0*dn_i      - rotj  * norm1    *dn_j;
      P3 dN01r = roti  * norm0*dn_i      - rotkj * normFix_j*dn_kj;
      P3 dN10r = rotki * normFix_i*dn_ki - rotj  * norm1    *dn_j;
      P3 dN11r = rotki * normFix_i*dn_ki - rotkj * normFix_j*dn_kj;

      P3 f00 = p_i - p_j;
      P3 f11 = q_i - q_j;

      P3 g00 = p_i - p_kj;
      P3 g11 = q_i - q_kj;
      
      P3 h00 = p_ki - p_j;
      P3 h11 = q_ki - q_j;

      P3 e00 = p_ki - p_kj;
      P3 e11 = q_ki - q_kj;

#ifndef __do_depth_only__

      P3 pi_t1 = roti*p_i_0+trai+a_i; 
      P3 pj_t1 = rotj*p_j_0+traj+b_j;

      P3 qi_t1 = roti*q_i_0+trai+a_i; 
      P3 qj_t1 = rotj*q_j_0+traj+b_j;

      P3 pkA_t1 = rotki*p_k_0A+traki+a_k;
      P3 pkB_t1 = rotkj*p_k_0B+trakj+b_k;

      P3 qkA_t1 = rotki*q_k_0A+traki+a_k;
      P3 qkB_t1 = rotkj*q_k_0B+trakj+b_k;

	  // in the paper is the 3d motion difference, so:  no penalty at global motion
  	  P3 f00_rt = pi_t1 - pj_t1 - f00;
      P3 f11_rt = qi_t1 - qj_t1 - f11;

      P3 g00_rt = pi_t1 - pkB_t1 - g00;
      P3 g11_rt = qi_t1 - qkB_t1 - g11;
      
      P3 h00_rt = pkA_t1 - pj_t1 - h00;
      P3 h11_rt = qkA_t1 - qj_t1 - h11;

      P3 e00_rt = pkA_t1 - pkB_t1 - e00;
      P3 e11_rt = qkA_t1 - qkB_t1 - e11;

	  
#endif
      ////////////////////////////////////////////////////////////////////////////////
      // if it is not connected at timestep t+1 why should it be at timestep t
      // take the maximum of the depth energies at the different timesteps:
      // the depth energy at timestep t
      *f00_it = sqrt ( ( (f00|f00) + (f11|f11) + (f00|f11) + gamma*(dN00|dN00) ) * Scalar(3.3) + epsilon);
      *f10_it = sqrt ( ( (h00|h00) + (h11|h11) + (h00|h11) + gamma*(dN10|dN10) ) * Scalar(3.3) + epsilon);
      *f01_it = sqrt ( ( (g00|g00) + (g11|g11) + (g00|g11) + gamma*(dN01|dN01) ) * Scalar(3.3) + epsilon);
      *f11_it = sqrt ( ( (e00|e00) + (e11|e11) + (e00|e11) + gamma*(dN11|dN11) ) * Scalar(3.3) + epsilon);

      // the depth energy at timestep t+1:
#ifndef __do_depth_only__
      Scalar f00_R = sqrt ( ( (f00_rt|f00_rt) + (f11_rt|f11_rt) + (f00_rt|f11_rt) + gamma*(dN00r|dN00r))* Scalar(3.3) + epsilon);
      Scalar f10_R = sqrt ( ( (h00_rt|h00_rt) + (h11_rt|h11_rt) + (h00_rt|h11_rt) + gamma*(dN10r|dN10r))* Scalar(3.3) + epsilon);
      Scalar f01_R = sqrt ( ( (g00_rt|g00_rt) + (g11_rt|g11_rt) + (g00_rt|g11_rt) + gamma*(dN01r|dN01r))* Scalar(3.3) + epsilon);
      Scalar f11_R = sqrt ( ( (e00_rt|e00_rt) + (e11_rt|e11_rt) + (e00_rt|e11_rt) + gamma*(dN11r|dN11r))* Scalar(3.3) + epsilon);
#endif


#pragma omp critical
      {
#ifndef __do_depth_only__

      *f00_it = (min(*f00_it, depthJump) + rotWeight * min(f00_R, rotJump) ) * *w_it;
      *f10_it = (min(*f10_it, depthJump) + rotWeight * min(f10_R, rotJump) ) * *w_it;
      *f01_it = (min(*f01_it, depthJump) + rotWeight * min(f01_R, rotJump) ) * *w_it;
      *f11_it = (min(*f11_it, depthJump) + rotWeight * min(f11_R, rotJump) ) * *w_it;

#else
// depth only
        *f00_it = (min(*f00_it, depthJump) ) * *w_it;
        *f11_it = (min(*f11_it, depthJump) ) * *w_it;
        *f01_it = (min(*f01_it, depthJump) ) * *w_it;
        *f10_it = (min(*f10_it, depthJump) ) * *w_it;
#endif

      score += max(max(max(*f00_it, *f01_it), *f10_it), *f11_it);

      if ( score != score || std::numeric_limits<Scalar>::infinity()==score )
      {
        printf ("0 vec as normal vector %f \n", score);
        mexEvalString("drawnow");
        int breakhere = 0;
      }

      }
    }
    return score;
  }

#ifdef __USE3D__
template<class Scalar>
Scalar
EvalEnergyFullFrame<Scalar>::
compute_score_Fuse2D( std::vector<int>& curr_ids, std::vector<int>& curr_ids2 )
#else
template<class Scalar>
Scalar
EvalEnergyFullFrame<Scalar>::
compute_score_Fuse( std::vector<int>& curr_ids, std::vector<int>& curr_ids2 )
#endif
  {
    if (normals == NULL) return 0;

    // HACK we should assume that here else 2d stuff is too complicated

    Scalar score (0.);

#pragma omp parallel for schedule (static)
    for (int www = 0; www < weights.size(); www++)
    {
      Scalar* w_it   = &(weights[www]);
      Scalar* f00_it = &(F00[www]);
      Scalar* f01_it = &(F01[www]);
      Scalar* f10_it = &(F10[www]);
      Scalar* f11_it = &(F11[www]);

      int* idl_it = &(Idl[www]);
      int* idk_it = &(Idk[www]);

      P3* eP_it = &(edgesP[www]);
      P3* eQ_it = &(edgesQ[www]);
      P3* cA_it = &(centersA[www]);
      P3* cB_it = &(centersB[www]);

      int idN_a = curr_ids[ *idk_it ];
      int idN_b = curr_ids[ *idl_it ];

      P3 norm0 = (*normals)[idN_a];
      P3 norm1 = (*normals)[idN_b];

      M3 roti = (*rot)[idN_a];
      M3 rotj = (*rot)[idN_b];
      P3 trai = (*tra)[idN_a];
      P3 traj = (*tra)[idN_b];

      int idK_a = curr_ids2[ *idk_it ];
      int idK_b = curr_ids2[ *idl_it ];
      P3 normFix_i = (*normals)[idK_a];
      P3 normFix_j = (*normals)[idK_b];

      Scalar dn_ki =  1. / sqrt(normFix_i|normFix_i);
      Scalar dn_kj =  1. / sqrt(normFix_j|normFix_j);

      M3 rotki = (*rot)[idK_a];
      M3 rotkj = (*rot)[idK_b];
      P3 traki = (*tra)[idK_a];
      P3 trakj = (*tra)[idK_b];

      P3 p2d ( *eP_it );
      P3 q2d ( *eQ_it );

#ifdef _use_patchCenters_
      P3 centerA ( *cA_it );
      P3 centerB ( *cB_it );

      P3 a_i  = centerA * (1. / (norm0 | centerA));
      P3 b_j  = centerB * (1. / (norm1 | centerB));
      P3 a_k  = centerA * (1. / (normFix_i | centerA));
      P3 b_k  = centerB * (1. / (normFix_j | centerB));
#endif

      Scalar dp_ki = 1. / max(normFix_i | p2d, 0.00000001 );
      Scalar dq_ki = 1. / max(normFix_i | q2d, 0.00000001 );
      Scalar dp_kj = 1. / max(normFix_j | p2d, 0.00000001 );
      Scalar dq_kj = 1. / max(normFix_j | q2d, 0.00000001 );
      Scalar dp_i  = 1. / max(norm0     | p2d, 0.00000001 );
      Scalar dq_i  = 1. / max(norm0     | q2d, 0.00000001 );
      Scalar dp_j  = 1. / max(norm1     | p2d, 0.00000001 );
      Scalar dq_j  = 1. / max(norm1     | q2d, 0.00000001 );


      //  patch 'depth'
      Scalar dn_i =  1. / max(sqrt(norm0|norm0), 0.00000001);
      Scalar dn_j =  1. / max(sqrt(norm1|norm1), 0.00000001);

      P3 p_i  = p2d * dp_i;
      P3 q_j  = q2d * dq_j;
      P3 p_j  = p2d * dp_j;
      P3 q_i  = q2d * dq_i;
      P3 p_ki = p2d * dp_ki;
      P3 q_ki = q2d * dq_ki;
      P3 p_kj = p2d * dp_kj;
      P3 q_kj = q2d * dq_kj;

#ifdef _use_patchCenters_
      ///////// 2nd part rotations:
      P3 p_i_0 = p_i - a_i;
      P3 q_i_0 = q_i - a_i;

      P3 p_j_0 = p_j - b_j;
      P3 q_j_0 = q_j - b_j;

      P3 p_k_0A = p_ki - a_k;
      P3 q_k_0A = q_ki - a_k;

      P3 p_k_0B = p_kj - b_k;
      P3 q_k_0B = q_kj - b_k;
#else
      P3 p_i_0 = p_i;
      P3 q_i_0 = q_i;

      P3 p_j_0 = p_j;
      P3 q_j_0 = q_j;

      P3 p_k_0A = p_ki;
      P3 q_k_0A = q_ki;

      P3 p_k_0B = p_kj;
      P3 q_k_0B = q_kj;
#endif

      // now add the normal constraint to get a linear prior
      P3 dN00 = norm0*dn_i      - norm1     *dn_j;
      P3 dN01 = norm0*dn_i      - normFix_j *dn_kj;
      P3 dN10 = normFix_i*dn_ki - norm1     *dn_j;
      P3 dN11 = normFix_i*dn_ki - normFix_j *dn_kj;

      P3 dN00r = roti  * norm0*dn_i      - rotj  * norm1    *dn_j;
      P3 dN01r = roti  * norm0*dn_i      - rotkj * normFix_j*dn_kj;
      P3 dN10r = rotki * normFix_i*dn_ki - rotj  * norm1    *dn_j;
      P3 dN11r = rotki * normFix_i*dn_ki - rotkj * normFix_j*dn_kj;

#ifdef __do_depth_only__

      P1 f00 ( getDisparity( dp_i, dp_j ) );
      P1 f10 ( getDisparity( dq_i, dq_j ) );
      P1 f01 = f00;
      P1 f11 = f10;

      P1 g00 ( getDisparity( dp_i, dp_k ));
      P1 g10 ( getDisparity( dq_i, dq_k ));
      P1 g01 = g00;
      P1 g11 = g10;
      
      P1 h00 ( getDisparity( dp_j, dp_k ));
      P1 h10 ( getDisparity( dq_j, dq_k ));
      P1 h01 = h00;
      P1 h11 = h10;

#else

#ifdef _use_patchCenters_
      P3 pi_t1 = roti*p_i_0+trai+a_i; 
      P3 pj_t1 = rotj*p_j_0+traj+b_j;

      P3 qi_t1 = roti*q_i_0+trai+a_i; 
      P3 qj_t1 = rotj*q_j_0+traj+b_j;

      P3 pkA_t1 = rotki*p_k_0A+traki+a_k;
      P3 pkB_t1 = rotkj*p_k_0B+trakj+b_k;

      P3 qkA_t1 = rotki*q_k_0A+traki+a_k;
      P3 qkB_t1 = rotkj*q_k_0B+trakj+b_k;
#else
      P3 pi_t1 = roti*p_i_0+trai; 
      P3 pj_t1 = rotj*p_j_0+traj;

      P3 qi_t1 = roti*q_i_0+trai; 
      P3 qj_t1 = rotj*q_j_0+traj;

      P3 pkA_t1 = rotki*p_k_0A+traki;
      P3 pkB_t1 = rotkj*p_k_0B+trakj;

      P3 qkA_t1 = rotki*q_k_0A+traki;
      P3 qkB_t1 = rotkj*q_k_0B+trakj;
#endif

      P2 f00(0.,0.), f11(0.,0.);
      P2 g00(0.,0.), g11(0.,0.);
      P2 h00(0.,0.), h11(0.,0.);
      P2 e00(0.,0.), e11(0.,0.);

      P2 f00_t, f00_rt, f11_t, f11_rt;
      P2 g00_t, g00_rt, g11_t, g11_rt;
      P2 h00_t, h00_rt, h11_t, h11_rt;
      P2 e00_t, e00_rt, e11_t, e11_rt;

      P2 p_i_t0_l, p_i_t0_r;
      P2 q_i_t0_l, q_i_t0_r;
      P2 p_i_t1_l, p_i_t1_r;
      P2 q_i_t1_l, q_i_t1_r;

      P2 p_j_t0_l, p_j_t0_r;
      P2 q_j_t0_l, q_j_t0_r;
      P2 p_j_t1_l, p_j_t1_r;
      P2 q_j_t1_l, q_j_t1_r;

      P2 p_ki_t0_l, p_ki_t0_r;
      P2 q_ki_t0_l, q_ki_t0_r;
      P2 p_kj_t0_l, p_kj_t0_r;
      P2 q_kj_t0_l, q_kj_t0_r;
      P2 p_kA_t1_l, p_kA_t1_r;
      P2 q_kA_t1_l, q_kA_t1_r;
      P2 p_kB_t1_l, p_kB_t1_r;
      P2 q_kB_t1_l, q_kB_t1_r;

      get2dProjections( p_i, p_i_t0_l, p_i_t0_r );
      get2dProjections( q_i, q_i_t0_l, q_i_t0_r );
      get2dProjections( pi_t1, p_i_t1_l, p_i_t1_r );
      get2dProjections( qi_t1, q_i_t1_l, q_i_t1_r );

      get2dProjections( p_j, p_j_t0_l, p_j_t0_r );
      get2dProjections( q_j, q_j_t0_l, q_j_t0_r );
      get2dProjections( pj_t1, p_j_t1_l, p_j_t1_r );
      get2dProjections( qj_t1, q_j_t1_l, q_j_t1_r );

      get2dProjections( p_ki, p_ki_t0_l, p_ki_t0_r );
      get2dProjections( q_ki, q_ki_t0_l, q_ki_t0_r );
      get2dProjections( p_kj, p_kj_t0_l, p_kj_t0_r );
      get2dProjections( q_kj, q_kj_t0_l, q_kj_t0_r );
      get2dProjections( pkA_t1, p_kA_t1_l, p_kA_t1_r );
      get2dProjections( qkA_t1, q_kA_t1_l, q_kA_t1_r );
      get2dProjections( pkB_t1, p_kB_t1_l, p_kB_t1_r );
      get2dProjections( qkB_t1, q_kB_t1_l, q_kB_t1_r );

      get2dDiff_Fast( p_i_t0_l, p_i_t1_l, p_i_t0_r, p_i_t1_r,   p_j_t0_l, p_j_t1_l, p_j_t0_r, p_j_t1_r,     f00, f00_t, f00_rt );
      get2dDiff_Fast( p_i_t0_l, p_i_t1_l, p_i_t0_r, p_i_t1_r,   p_kj_t0_l, p_kB_t1_l, p_kj_t0_r, p_kB_t1_r,   g00, g00_t, g00_rt );
      
      get2dDiff_Fast( p_ki_t0_l, p_kA_t1_l, p_ki_t0_r, p_kA_t1_r,  p_j_t0_l, p_j_t1_l, p_j_t0_r, p_j_t1_r,    h00, h00_t, h00_rt );
      get2dDiff_Fast( p_ki_t0_l, p_kA_t1_l, p_ki_t0_r, p_kA_t1_r,  p_kj_t0_l, p_kB_t1_l, p_kj_t0_r, p_kB_t1_r,  e00, e00_t, e00_rt );

      get2dDiff_Fast( q_i_t0_l, q_i_t1_l, q_i_t0_r, q_i_t1_r,   q_j_t0_l, q_j_t1_l, q_j_t0_r, q_j_t1_r,     f11, f11_t, f11_rt );
      get2dDiff_Fast( q_i_t0_l, q_i_t1_l, q_i_t0_r, q_i_t1_r,   q_kj_t0_l, q_kB_t1_l, q_kj_t0_r, q_kB_t1_r,   g11, g11_t, g11_rt );

      get2dDiff_Fast( q_ki_t0_l, q_kA_t1_l, q_ki_t0_r, q_kA_t1_r,  q_j_t0_l, q_j_t1_l, q_j_t0_r, q_j_t1_r,    h11, h11_t, h11_rt );
      get2dDiff_Fast( q_ki_t0_l, q_kA_t1_l, q_ki_t0_r, q_kA_t1_r,  q_kj_t0_l, q_kB_t1_l, q_kj_t0_r, q_kB_t1_r,  e11, e11_t, e11_rt );

      // combine both motions 
      P4 f00_m ( f00_t[0], f00_t[1], f00_rt[0], f00_rt[1] );
      P4 f11_m ( f11_t[0], f11_t[1], f11_rt[0], f11_rt[1] );
      P4 g00_m ( g00_t[0], g00_t[1], g00_rt[0], g00_rt[1] );
      P4 g11_m ( g11_t[0], g11_t[1], g11_rt[0], g11_rt[1] );
      P4 h00_m ( h00_t[0], h00_t[1], h00_rt[0], h00_rt[1] );
      P4 h11_m ( h11_t[0], h11_t[1], h11_rt[0], h11_rt[1] );
      P4 e00_m ( e00_t[0], e00_t[1], e00_rt[0], e00_rt[1] );
      P4 e11_m ( e11_t[0], e11_t[1], e11_rt[0], e11_rt[1] );

#endif
      ////////////////////////////////////////////////////////////////////////////////
      *f00_it = sqrt ( ( (f00|f00) + (f11|f11) + (f00|f11) + gamma*(dN00|dN00) ) * Scalar(3.3) + epsilon);
      *f10_it = sqrt ( ( (h00|h00) + (h11|h11) + (h00|h11) + gamma*(dN10|dN10) ) * Scalar(3.3) + epsilon);
      *f01_it = sqrt ( ( (g00|g00) + (g11|g11) + (g00|g11) + gamma*(dN01|dN01) ) * Scalar(3.3) + epsilon);
      *f11_it = sqrt ( ( (e00|e00) + (e11|e11) + (e00|e11) + gamma*(dN11|dN11) ) * Scalar(3.3) + epsilon);

      // the depth energy at timestep t+1:
#ifndef __do_depth_only__
       Scalar f00_M = sqrt ( ((f00_m|f00_m) + (f11_m|f11_m) + (f00_m|f11_m) + gamma*(dN00r|dN00r)) * Scalar(3.3) + epsilon);
       Scalar f10_M = sqrt ( ((h00_m|h00_m) + (h11_m|h11_m) + (h00_m|h11_m) + gamma*(dN10r|dN10r)) * Scalar(3.3) + epsilon);
       Scalar f01_M = sqrt ( ((g00_m|g00_m) + (g11_m|g11_m) + (g00_m|g11_m) + gamma*(dN01r|dN01r)) * Scalar(3.3) + epsilon);
       Scalar f11_M = sqrt ( ((e00_m|e00_m) + (e11_m|e11_m) + (e00_m|e11_m) + gamma*(dN11r|dN11r)) * Scalar(3.3) + epsilon);
#endif

#pragma omp critical
      {
#ifndef __do_depth_only__
      *f00_it = (min(*f00_it, depthJump) + rotWeight * min(f00_M, rotJump) ) * *w_it;
      *f10_it = (min(*f10_it, depthJump) + rotWeight * min(f10_M, rotJump) ) * *w_it;
      *f01_it = (min(*f01_it, depthJump) + rotWeight * min(f01_M, rotJump) ) * *w_it;
      *f11_it = (min(*f11_it, depthJump) + rotWeight * min(f11_M, rotJump) ) * *w_it;
#else

        *f00_it = (min(*f00_it, depthJump) ) * *w_it;
        *f11_it = (min(*f11_it, depthJump) ) * *w_it;
        *f01_it = (min(*f01_it, depthJump) ) * *w_it;
        *f10_it = (min(*f10_it, depthJump) ) * *w_it;

#endif

      score += max(max(max(*f00_it, *f01_it), *f10_it), *f11_it);

      if ( score != score || std::numeric_limits<Scalar>::infinity()==score )
      {
        printf ("0Vec as Normal Vector %f \n", score);
        mexEvalString("drawnow");
        int breakhere = 0;
      }

      }
    }
    return score;
  }
  
  /// faster version if only the current solution is to be evaluated - called from getEnergy , etc.
#ifdef __USE3D__
template<class Scalar>
Scalar
EvalEnergyFullFrame<Scalar>::
compute_score_Fuse2D00( std::vector<int>& curr_ids, std::vector<int>& curr_ids2 )
#else
template<class Scalar>
Scalar
EvalEnergyFullFrame<Scalar>::
compute_score_Fuse00( std::vector<int>& curr_ids )
#endif
  {
    if (normals == NULL) return 0;
    Scalar score (0.);

//#pragma omp parallel for schedule (static)
    for (int www = 0; www < weights.size(); www++)
    {
      Scalar* w_it   = &(weights[www]);
      Scalar* f00_it = &(F00[www]);

      int* idl_it = &(Idl[www]);
      int* idk_it = &(Idk[www]);

      P3* eP_it = &(edgesP[www]);
      P3* eQ_it = &(edgesQ[www]);
      P3* cA_it = &(centersA[www]);
      P3* cB_it = &(centersB[www]);

      int idN_a = curr_ids[ *idk_it ];
      int idN_b = curr_ids[ *idl_it ];

      P3 norm0 = (*normals)[idN_a];
      P3 norm1 = (*normals)[idN_b];

      M3 roti = (*rot)[idN_a];
      M3 rotj = (*rot)[idN_b];
      P3 trai = (*tra)[idN_a];
      P3 traj = (*tra)[idN_b];

      P3 p2d ( *eP_it );
      P3 q2d ( *eQ_it );

#ifdef _use_patchCenters_
      P3 centerA ( *cA_it );
      P3 centerB ( *cB_it );

      P3 a_i  = centerA * (1. / (norm0 | centerA));
      P3 b_j  = centerB * (1. / (norm1 | centerB));
#else
      P3 a_i  (0.,0.,0.);
      P3 b_j  (0.,0.,0.);
#endif

      Scalar dp_i  = 1. / max(norm0     | p2d, 0.00000001 );
      Scalar dq_i  = 1. / max(norm0     | q2d, 0.00000001 );
      Scalar dp_j  = 1. / max(norm1     | p2d, 0.00000001 );
      Scalar dq_j  = 1. / max(norm1     | q2d, 0.00000001 );

      //  patch 'depth'
      Scalar dn_i =  1. / max(sqrt(norm0|norm0), 0.00000001);
      Scalar dn_j =  1. / max(sqrt(norm1|norm1), 0.00000001);

      P3 p_i  = p2d * dp_i;
      P3 q_j  = q2d * dq_j;
      P3 p_j  = p2d * dp_j;
      P3 q_i  = q2d * dq_i;

      ///////// 2nd part rotations:
      P3 p_i_0 = p_i - a_i;
      P3 q_i_0 = q_i - a_i;

      P3 p_j_0 = p_j - b_j;
      P3 q_j_0 = q_j - b_j;

      // now add the normal constraint to get a linear prior
      P3 dN00  = norm0*dn_i         - norm1     *dn_j;
      P3 dN00r = roti  * norm0*dn_i - rotj  * norm1    *dn_j;

#ifdef __do_depth_only__

      P1 f00 ( getDisparity( dp_i, dp_j ) );
      P1 f10 ( getDisparity( dq_i, dq_j ) );
      P1 f01 = f00;
      P1 f11 = f10;

#else

      P3 pi_t1 = roti*p_i_0+trai+a_i; 
      P3 pj_t1 = rotj*p_j_0+traj+b_j;

      P3 qi_t1 = roti*q_i_0+trai+a_i; 
      P3 qj_t1 = rotj*q_j_0+traj+b_j;

      P2 f00, f00_t, f00_rt, f11, f11_t, f11_rt;

      P2 p_i_t0_l, p_i_t0_r;
      P2 q_i_t0_l, q_i_t0_r;
      P2 p_i_t1_l, p_i_t1_r;
      P2 q_i_t1_l, q_i_t1_r;

      P2 p_j_t0_l, p_j_t0_r;
      P2 q_j_t0_l, q_j_t0_r;
      P2 p_j_t1_l, p_j_t1_r;
      P2 q_j_t1_l, q_j_t1_r;

      get2dProjections( p_i, p_i_t0_l, p_i_t0_r );
      get2dProjections( q_i, q_i_t0_l, q_i_t0_r );
      get2dProjections( pi_t1, p_i_t1_l, p_i_t1_r );
      get2dProjections( qi_t1, q_i_t1_l, q_i_t1_r );

      get2dProjections( p_j, p_j_t0_l, p_j_t0_r );
      get2dProjections( q_j, q_j_t0_l, q_j_t0_r );
      get2dProjections( pj_t1, p_j_t1_l, p_j_t1_r );
      get2dProjections( qj_t1, q_j_t1_l, q_j_t1_r );

      get2dDiff_Fast( p_i_t0_l, p_i_t1_l, p_i_t0_r, p_i_t1_r,   p_j_t0_l, p_j_t1_l, p_j_t0_r, p_j_t1_r,     f00, f00_t, f00_rt );
      get2dDiff_Fast( q_i_t0_l, q_i_t1_l, q_i_t0_r, q_i_t1_r,   q_j_t0_l, q_j_t1_l, q_j_t0_r, q_j_t1_r,     f11, f11_t, f11_rt );

      // combine both motions 
      P4 f00_m ( f00_t[0], f00_t[1], f00_rt[0], f00_rt[1] );
      P4 f11_m ( f11_t[0], f11_t[1], f11_rt[0], f11_rt[1] );
#endif
      ////////////////////////////////////////////////////////////////////////////////
      *f00_it = sqrt ( ( (f00|f00) + (f11|f11) + (f00|f11) + gamma*(dN00|dN00) ) * Scalar(3.3) + epsilon);

      // the depth energy at timestep t+1:
#ifndef __do_depth_only__
       Scalar f00_M = sqrt ( ((f00_m|f00_m) + (f11_m|f11_m) + (f00_m|f11_m) + gamma*(dN00r|dN00r)) * Scalar(3.3) + epsilon);
#endif


//#pragma omp critical
      {
#ifndef __do_depth_only__
      *f00_it = (min(*f00_it, depthJump) + rotWeight * min(f00_M, rotJump) ) * *w_it;
#else
        *f00_it = (min(*f00_it, depthJump) ) * *w_it;
#endif
      score += *f00_it;
      }
    }
    return score;
  }
#endif
