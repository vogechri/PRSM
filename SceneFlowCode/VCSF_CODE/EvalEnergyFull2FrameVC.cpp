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

#ifndef __ENERGY_ROTTRANORVC__cpp
#define __ENERGY_ROTTRANORVC__cpp

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

/// 3d not fully implemented so on 
#define _usePure_2d
/// just above 0 
#define tinyValue 0.00000001

using namespace std;
using namespace Math;

/// segimg must be flipped
template<class Scalar>
Scalar
EvalEnergyFullFrame<Scalar>::
  compute_score_Fuse_local( const std::vector<int>& curr_ids, const std::vector<int>& curr_ids2, const std::vector<int>& localIds, 
    std::vector<Scalar>& _F00, std::vector<Scalar>& _F01, std::vector<Scalar>& _F10, std::vector<Scalar>& _F11,
    std::vector<int>& _idk, std::vector<int>& _idl, bool fillVecs )
  {
    if (normals == NULL) return 0;
    Scalar score (0.);

    for (int www = 0; www < weights.size(); www++)
    {
      const int* const idl_it = &(Idl[www]);
      const int* const idk_it = &(Idk[www]);

      const Scalar* const w_it   = &(weights[www]);

      // here local to avoid overriding in parallel from outside :)
      Scalar temp_f00(0),temp_f10(0),temp_f01(0),temp_f11(0);
      Scalar* const f00_it = (fillVecs) ? &(temp_f00) : &(F00[www]);
      Scalar* const f01_it = (fillVecs) ? &(temp_f01) : &(F01[www]);
      Scalar* const f10_it = (fillVecs) ? &(temp_f10) : &(F10[www]);
      Scalar* const f11_it = (fillVecs) ? &(temp_f11) : &(F11[www]);

      if(localIds[*idk_it] <0 && localIds[*idl_it] <0) //skip if both oobs
      {
        *f00_it =0;*f10_it =0;*f01_it =0;*f11_it =0;
        continue;
      }

      P3* eP_it = &(edgesP[www]);
      P3* eQ_it = &(edgesQ[www]);
//      P3* cA_it = &(centersA[www]);
//      P3* cB_it = &(centersB[www]);

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
//      P3 addki = normFix_a * (dn_ki * gamma);
//      P3 addkj = normFix_b * (dn_kj * gamma);

      M3 rotki = (*rot)[idK_a];
      M3 rotkj = (*rot)[idK_b];
      P3 traki = (*tra)[idK_a];
      P3 trakj = (*tra)[idK_b];

      P3 p2d ( *eP_it );
      P3 q2d ( *eQ_it );

//      P3 centerA ( *cA_it );
//      P3 centerB ( *cB_it );

//      P3 a_i  = centerA * (1. / (norm0 | centerA));
//      P3 b_j  = centerB * (1. / (norm1 | centerB));
//      P3 a_k  = centerA * (1. / (normFix_i | centerA));
//      P3 b_k  = centerB * (1. / (normFix_j | centerB));

      Scalar dp_ki = 1. / max(normFix_i | p2d, tinyValue );
      Scalar dq_ki = 1. / max(normFix_i | q2d, tinyValue );
      Scalar dp_kj = 1. / max(normFix_j | p2d, tinyValue );
      Scalar dq_kj = 1. / max(normFix_j | q2d, tinyValue );
      Scalar dp_i  = 1. / max(norm0     | p2d, tinyValue );
      Scalar dq_i  = 1. / max(norm0     | q2d, tinyValue );
      Scalar dp_j  = 1. / max(norm1     | p2d, tinyValue );
      Scalar dq_j  = 1. / max(norm1     | q2d, tinyValue );


      //  patch 'depth'
      Scalar dn_i =  1. / max(sqrt(norm0|norm0), tinyValue);
      Scalar dn_j =  1. / max(sqrt(norm1|norm1), tinyValue);

//      P3 addI = norm0 * (dn_i * gamma);
//      P3 addJ = norm1 * (dn_j * gamma);

      P3 p_i  = p2d * dp_i;
      P3 q_j  = q2d * dq_j;
      P3 p_j  = p2d * dp_j;
      P3 q_i  = q2d * dq_i;
      P3 p_ki = p2d * dp_ki;
      P3 q_ki = q2d * dq_ki;
      P3 p_kj = p2d * dp_kj;
      P3 q_kj = q2d * dq_kj;


      ///////// 2nd part rotations:
      P3 p_i_0 = p_i;// - a_i;
      P3 q_i_0 = q_i;// - a_i;

      P3 p_j_0 = p_j;// - b_j;
      P3 q_j_0 = q_j;// - b_j;

      P3 p_k_0A = p_ki;// - a_k;
      P3 q_k_0A = q_ki;// - a_k;

      P3 p_k_0B = p_kj;// - b_k;
      P3 q_k_0B = q_kj;// - b_k;

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

      P3 pi_t1 = roti*p_i_0+trai;//+a_i; 
      P3 pj_t1 = rotj*p_j_0+traj;//+b_j;

      P3 qi_t1 = roti*q_i_0+trai;//+a_i; 
      P3 qj_t1 = rotj*q_j_0+traj;//+b_j;

      P3 pkA_t1 = rotki*p_k_0A+traki;//+a_k;
      P3 pkB_t1 = rotkj*p_k_0B+trakj;//+b_k;

      P3 qkA_t1 = rotki*q_k_0A+traki;//+a_k;
      P3 qkB_t1 = rotkj*q_k_0B+trakj;//+b_k;


#ifdef _usePure_2d

#ifdef _x2DispSmoothing_
      P4 f00(0,0,0,0), f11(0,0,0,0);
      P4 g00(0,0,0,0), g11(0,0,0,0);
      P4 h00(0,0,0,0), h11(0,0,0,0);
      P4 e00(0,0,0,0), e11(0,0,0,0);
#else
      P2 f00, f11;
      P2 g00, g11;
      P2 h00, h11;
      P2 e00, e11;
#endif
      P2 f00_t(0,0), f00_rt(0,0), f11_t(0,0), f11_rt(0,0);
      P2 g00_t(0,0), g00_rt(0,0), g11_t(0,0), g11_rt(0,0);
      P2 h00_t(0,0), h00_rt(0,0), h11_t(0,0), h11_rt(0,0);
      P2 e00_t(0,0), e00_rt(0,0), e11_t(0,0), e11_rt(0,0);

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

#ifdef _x2DispSmoothing_
      f00_rt = P2(0,0); f11_rt = P2(0,0);
      g00_rt = P2(0,0); g11_rt = P2(0,0);
      h00_rt = P2(0,0); h11_rt = P2(0,0);
      e00_rt = P2(0,0); e11_rt = P2(0,0);

      get2dDiff_x2disp( p_i_t0_l, p_i_t1_l, p_i_t0_r, p_i_t1_r,   p_j_t0_l, p_j_t1_l, p_j_t0_r, p_j_t1_r,     f00, f00_t);
      get2dDiff_x2disp( p_i_t0_l, p_i_t1_l, p_i_t0_r, p_i_t1_r,   p_kj_t0_l, p_kB_t1_l, p_kj_t0_r, p_kB_t1_r,   g00, g00_t );
      
      get2dDiff_x2disp( p_ki_t0_l, p_kA_t1_l, p_ki_t0_r, p_kA_t1_r,  p_j_t0_l, p_j_t1_l, p_j_t0_r, p_j_t1_r,    h00, h00_t );
      get2dDiff_x2disp( p_ki_t0_l, p_kA_t1_l, p_ki_t0_r, p_kA_t1_r,  p_kj_t0_l, p_kB_t1_l, p_kj_t0_r, p_kB_t1_r,  e00, e00_t );

      get2dDiff_x2disp( q_i_t0_l, q_i_t1_l, q_i_t0_r, q_i_t1_r,   q_j_t0_l, q_j_t1_l, q_j_t0_r, q_j_t1_r,     f11, f11_t );
      get2dDiff_x2disp( q_i_t0_l, q_i_t1_l, q_i_t0_r, q_i_t1_r,   q_kj_t0_l, q_kB_t1_l, q_kj_t0_r, q_kB_t1_r,   g11, g11_t );

      get2dDiff_x2disp( q_ki_t0_l, q_kA_t1_l, q_ki_t0_r, q_kA_t1_r,  q_j_t0_l, q_j_t1_l, q_j_t0_r, q_j_t1_r,    h11, h11_t );
      get2dDiff_x2disp( q_ki_t0_l, q_kA_t1_l, q_ki_t0_r, q_kA_t1_r,  q_kj_t0_l, q_kB_t1_l, q_kj_t0_r, q_kB_t1_r,  e11, e11_t );
#else
      get2dDiff_vc( p_i_t0_l, p_i_t1_l, p_i_t0_r, p_i_t1_r,   p_j_t0_l, p_j_t1_l, p_j_t0_r, p_j_t1_r,     f00, f00_t, f00_rt );
      get2dDiff_vc( p_i_t0_l, p_i_t1_l, p_i_t0_r, p_i_t1_r,   p_kj_t0_l, p_kB_t1_l, p_kj_t0_r, p_kB_t1_r,   g00, g00_t, g00_rt );
      
      get2dDiff_vc( p_ki_t0_l, p_kA_t1_l, p_ki_t0_r, p_kA_t1_r,  p_j_t0_l, p_j_t1_l, p_j_t0_r, p_j_t1_r,    h00, h00_t, h00_rt );
      get2dDiff_vc( p_ki_t0_l, p_kA_t1_l, p_ki_t0_r, p_kA_t1_r,  p_kj_t0_l, p_kB_t1_l, p_kj_t0_r, p_kB_t1_r,  e00, e00_t, e00_rt );

      get2dDiff_vc( q_i_t0_l, q_i_t1_l, q_i_t0_r, q_i_t1_r,   q_j_t0_l, q_j_t1_l, q_j_t0_r, q_j_t1_r,     f11, f11_t, f11_rt );
      get2dDiff_vc( q_i_t0_l, q_i_t1_l, q_i_t0_r, q_i_t1_r,   q_kj_t0_l, q_kB_t1_l, q_kj_t0_r, q_kB_t1_r,   g11, g11_t, g11_rt );

      get2dDiff_vc( q_ki_t0_l, q_kA_t1_l, q_ki_t0_r, q_kA_t1_r,  q_j_t0_l, q_j_t1_l, q_j_t0_r, q_j_t1_r,    h11, h11_t, h11_rt );
      get2dDiff_vc( q_ki_t0_l, q_kA_t1_l, q_ki_t0_r, q_kA_t1_r,  q_kj_t0_l, q_kB_t1_l, q_kj_t0_r, q_kB_t1_r,  e11, e11_t, e11_rt );
#endif

      // combine both motions 
      P4 f00_m ( f00_t[0], f00_t[1], f00_rt[0], f00_rt[1] );
      P4 f11_m ( f11_t[0], f11_t[1], f11_rt[0], f11_rt[1] );
      P4 g00_m ( g00_t[0], g00_t[1], g00_rt[0], g00_rt[1] );
      P4 g11_m ( g11_t[0], g11_t[1], g11_rt[0], g11_rt[1] );
      P4 h00_m ( h00_t[0], h00_t[1], h00_rt[0], h00_rt[1] );
      P4 h11_m ( h11_t[0], h11_t[1], h11_rt[0], h11_rt[1] );
      P4 e00_m ( e00_t[0], e00_t[1], e00_rt[0], e00_rt[1] );
      P4 e11_m ( e11_t[0], e11_t[1], e11_rt[0], e11_rt[1] );
#else

      P2 f00, f11;
      P2 g00, g11;
      P2 h00, h11;
      P2 e00, e11;
      P3 f00_rt, f11_rt;
      P3 g00_rt, g11_rt;
      P3 h00_rt, h11_rt;
      P3 e00_rt, e11_rt;

      get2d_3dDiff( p_i, p_j, pi_t1, pj_t1, f00, f00_rt );
      get2d_3dDiff( q_i, q_j, qi_t1, qj_t1, f11, f11_rt );

      // pk of right segment
      get2d_3dDiff( p_i, p_k, pi_t1, pkB_t1, g00, g00_rt );
      get2d_3dDiff( q_i, q_k, qi_t1, qkB_t1, g11, g11_rt );

      // pk of left segment
      get2d_3dDiff( p_k, p_j, pkA_t1, pj_t1, h00, h00_rt );
      get2d_3dDiff( q_k, q_j, qkA_t1, qj_t1, h11, h11_rt );

      //
      get2d_3dDiff( p_k, p_k, pkA_t1, pkB_t1, e00, e00_rt );
      get2d_3dDiff( q_k, q_k, qkA_t1, qkB_t1, e11, e11_rt );

#endif

#endif
      ////////////////////////////////////////////////////////////////////////////////
      // if it is not connected at timestep t+1 why should it be at timestep t
      // take the maximum of the depth energies at the different timesteps:
      // the depth energy at timestep t
#ifdef __AllInOne__
      *f00_it = sqrt ( ( (f00|f00) + (f11|f11) + (f00|f11) + gamma*(dN00|dN00) + rotWeight*((f00_m|f00_m) + (f11_m|f11_m) + (f00_m|f11_m) + gamma*(dN00r|dN00r)) ) * Scalar(3.3) + epsilon);
      *f10_it = sqrt ( ( (h00|h00) + (h11|h11) + (h00|h11) + gamma*(dN10|dN10) + rotWeight*((h00_m|h00_m) + (h11_m|h11_m) + (h00_m|h11_m) + gamma*(dN10r|dN10r)) ) * Scalar(3.3) + epsilon);
      *f01_it = sqrt ( ( (g00|g00) + (g11|g11) + (g00|g11) + gamma*(dN01|dN01) + rotWeight*((g00_m|g00_m) + (g11_m|g11_m) + (g00_m|g11_m) + gamma*(dN01r|dN01r)) ) * Scalar(3.3) + epsilon);
      *f11_it = sqrt ( ( (e00|e00) + (e11|e11) + (e00|e11) + gamma*(dN11|dN11) + rotWeight*((e00_m|e00_m) + (e11_m|e11_m) + (e00_m|e11_m) + gamma*(dN11r|dN11r)) ) * Scalar(3.3) + epsilon);
#else

      *f00_it = sqrt ( ( (f00|f00) + (f11|f11) + (f00|f11) + gamma*(dN00|dN00) ) * Scalar(3.3) + epsilon);
      *f10_it = sqrt ( ( (h00|h00) + (h11|h11) + (h00|h11) + gamma*(dN10|dN10) ) * Scalar(3.3) + epsilon);
      *f01_it = sqrt ( ( (g00|g00) + (g11|g11) + (g00|g11) + gamma*(dN01|dN01) ) * Scalar(3.3) + epsilon);
      *f11_it = sqrt ( ( (e00|e00) + (e11|e11) + (e00|e11) + gamma*(dN11|dN11) ) * Scalar(3.3) + epsilon);

      // the depth energy at timestep t+1:
#ifndef __do_depth_only__

#ifdef _usePure_2d

       Scalar f00_M = sqrt ( ((f00_m|f00_m) + (f11_m|f11_m) + (f00_m|f11_m) + gamma*(dN00r|dN00r)) * Scalar(3.3) + epsilon);
       Scalar f10_M = sqrt ( ((h00_m|h00_m) + (h11_m|h11_m) + (h00_m|h11_m) + gamma*(dN10r|dN10r)) * Scalar(3.3) + epsilon);
       Scalar f01_M = sqrt ( ((g00_m|g00_m) + (g11_m|g11_m) + (g00_m|g11_m) + gamma*(dN01r|dN01r)) * Scalar(3.3) + epsilon);
       Scalar f11_M = sqrt ( ((e00_m|e00_m) + (e11_m|e11_m) + (e00_m|e11_m) + gamma*(dN11r|dN11r)) * Scalar(3.3) + epsilon);

#else
       Scalar f00_T(0), f10_T(0), f01_T(0), f11_T(0);
#endif
#endif
#endif

      {
#ifndef __AllInOne__

#ifndef __do_depth_only__
#ifdef __do_max_of_both__

      *f00_it = (max(min(*f00_it, depthJump), min(f00_T, depthJump)) + rotWeight*min(f00_R, rotJump)) * *w_it;
      *f10_it = (max(min(*f10_it, depthJump), min(f10_T, depthJump)) + rotWeight*min(f10_R, rotJump)) * *w_it;
      *f01_it = (max(min(*f01_it, depthJump), min(f01_T, depthJump)) + rotWeight*min(f01_R, rotJump)) * *w_it;
      *f11_it = (max(min(*f11_it, depthJump), min(f11_T, depthJump)) + rotWeight*min(f11_R, rotJump)) * *w_it;
#else

      *f00_it = (min(*f00_it, depthJump) + rotWeight * min(f00_M, rotJump) ) * *w_it;
      *f10_it = (min(*f10_it, depthJump) + rotWeight * min(f10_M, rotJump) ) * *w_it;
      *f01_it = (min(*f01_it, depthJump) + rotWeight * min(f01_M, rotJump) ) * *w_it;
      *f11_it = (min(*f11_it, depthJump) + rotWeight * min(f11_M, rotJump) ) * *w_it;

#endif
#else

        *f00_it = (min(*f00_it, depthJump) ) * *w_it;
        *f11_it = (min(*f11_it, depthJump) ) * *w_it;
        *f01_it = (min(*f01_it, depthJump) ) * *w_it;
        *f10_it = (min(*f10_it, depthJump) ) * *w_it;

#endif

#else // all in one:
        *f00_it = (min(*f00_it, depthJump) ) * *w_it;
        *f11_it = (min(*f11_it, depthJump) ) * *w_it;
        *f01_it = (min(*f01_it, depthJump) ) * *w_it;
        *f10_it = (min(*f10_it, depthJump) ) * *w_it;
#endif

if (fillVecs)
{
    _F00.push_back(*f00_it);
    _F01.push_back(*f01_it); 
    _F10.push_back(*f10_it); 
    _F11.push_back(*f11_it);
    _idk.push_back( localIds[*idk_it] );
    _idl.push_back( localIds[*idl_it] );
}
      score += max(max(max(*f00_it, *f01_it), *f10_it), *f11_it);

      if ( score != score || std::numeric_limits<Scalar>::infinity()==score )
      {
        printf ("Infinity in local smoothness %f \n", score);
//        mexEvalString("drawnow");
        int breakhere = 0;
        P3 crap = (*normals)[-1];
        score = 10000.;
      }

      }
    }
    return score;
  };

#undef tinyValue
#endif