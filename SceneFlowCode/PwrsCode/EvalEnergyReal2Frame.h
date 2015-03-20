////////////////////////////////////////////////////////
////////// First version of faster iterations///////////
////////////////////////////////////////////////////////
#ifndef __ENERGY_ROTTRANORRTF__
#define __ENERGY_ROTTRANORRTF__

#include "Math/VectorT.h"
#include "Math/Mat3x3T.h"

#ifndef _NO_OPENMP
#include <omp.h>
#endif
#include "mex.h" // include after VectorT and Mat3x3T

//#define _Debug_out_

#include <map>
#include <vector>

//#define __USE3D__

//#define __do_depth_only__
//#undef __do_depth_only__

#define _usePure_2d
//#undef _usePure_2d

#define __combi2dSmoothMotion__


using namespace std;
using namespace Math;

///////////////////////////////////////////////////////////////////////////////

/*
allows : 
setNormals( std::vector<P3>* normals_ ) 
prepare( int nSegments, mxArray* edges_,  mxArray* weights_ )
compute_score( int id, std::vector<int>& curr_ids )
std::vector<Scalar>& getF11()
etc.
*/
/*! provides the possiblity to evaluate the edge energy bz iterating over all edges, stores per edge energy in array 
* needs to store the current solution: seg i-j : score and p3d/ (p2d ?) 
* 
* store pixel pairs : i,j: seg i,j
*/
//template<typename Scalar> class EvalEnergyRealFrame
template<typename Scalar> class EvalEnergyFullFrame
{
public:

  typedef unsigned int iType;
  typedef Math::VectorT<Scalar, 4> P4;
  typedef Math::VectorT<Scalar, 3> P3;
  typedef Math::VectorT<Scalar, 2> P2;
  typedef Math::VectorT<Scalar, 1> P1;
  typedef Math::Mat3x3T<Scalar>    M3; 

  EvalEnergyFullFrame(): gamma(1.0), epsilon(0.00), normals(NULL), tra(), rot(), rotJump(20.), depthJump(20.), nSegments(0), weightMap(),
    epix(0), Pl(1.,0.,0., 0.,1.,0., 0.,0.,1.), Pr(1.,0.,0., 0.,1.,0., 0.,0.,1.), pr(0.,0.,0.), simpleProjection(false), iK(1., 0., 0., 0., 1., 0., 0.,0.,1.)
  {};

  ~EvalEnergyFullFrame(){};

  void set2dMotionMatrix ( Scalar* K, Scalar* Rot, Scalar* mc, Scalar pixJump_ = Scalar(1.) ) 
  { 
    Pl = M3(K);
    Pr = M3(K) * M3(Rot);
    pr = M3(K) * P3(mc);

// assumes rot = eye(3);
//    P3 motion2d = M3(K) *  M3(Rot) * P3(mc);
    P3 motion2d = M3(K) * P3(mc);
    epix = motion2d[0] / pixJump_;

    simpleProjection = false;
    if (Rot[0] == Scalar(1.) && Rot[4] == Scalar(1.) && Rot[8] == Scalar(1.) )
      simpleProjection = true;

    iK=M3(K);
    //    iK=M3(K_[0], K_[3], K_[6], K_[1], K_[4], K_[7], K_[2], K_[5], K_[8]);
    iK.invert();
  };

  void setEpsilon (Scalar epsilon_)           {epsilon = epsilon_;}

  void setGamma (Scalar gamma_)               {gamma = gamma_;}

  void setRotWeight (Scalar rotWeight_)       { rotWeight = rotWeight_;};

  void setRotJump   (Scalar rotJump_)         {rotJump = rotJump_;}
  void setDepthJump (Scalar depthJump_)       {depthJump = depthJump_;}

  void setNormals( const std::vector<P3>* normals_ ) {normals = normals_;}
  void setRotTra(  const std::vector<M3>* rot_, const std::vector<P3>* tra_ ) {tra = tra_;rot = rot_;}

  const std::vector<int>& getIdl()  { return Idl; }
  const std::vector<int>& getIdk()  { return Idk; }

  const std::vector<Scalar>& getF11()  { return F11; }
  const std::vector<Scalar>& getF00()  { return F00; }
  const std::vector<Scalar>& getF10()  { return F10; }
  const std::vector<Scalar>& getF01()  { return F01; }

  const std::vector<P3>& getTra()      { return tra; }
  const std::vector<M3>& getRot()      { return rot; }
//  const std::vector<P3>& getCenters()  { return centers; }

  void setHalfEdges( Scalar* _halfEdgeX, Scalar* _halfEdgeY) {halfEdgeX = _halfEdgeX; halfEdgeY = _halfEdgeY;}
  /// weights on the cross egges in a 8 neighbourhood
  void setCrossHalfEdges( Scalar* _halfEdgeXY, Scalar* _halfEdgeiXY) {halfEdgeXY = _halfEdgeXY; halfEdgeiXY = _halfEdgeiXY;}

  /// segimg must be flipped
  void buildWeights(int* segImg, int width, int height, int nSegments)
  {
//    const int    dx4[4] = {-1,0,1,0};
//    const int    dy4[4] = {0,1,0,-1};
//    const Scalar we4[4] = {1.,1.,1.,1.};//weights

    edgeWeightMap.clear();
    edgeWeightMap.resize(nSegments);
    weightMap.clear();
    weightMap.reserve( nSegments );
    weightMap.resize( nSegments );

    const int    dx8[8] = {-1,0,1, 0, 1, -1, -1,  1};
    const int    dy8[8] = { 0,1,0,-1, 1, -1,  1, -1};
    const Scalar we8[8] = {1.,1.,1.,1., 1./sqrt(2.), 1./sqrt(2.), 1./sqrt(2.), 1./sqrt(2.)};//weights

    const int neighIts = (halfEdgeXY == NULL) ? 4:8;
    const int* dx= dx8;
    const int* dy= dy8;
    const Scalar* we= we8;

    int maxSmallWidth  = width *2 + 2 + 2;
    int maxSmallHeight = height*2 + 2 + 2;

    //stores per segments all th ep2d points already copied to prevent pixel to be used twice
    std::vector< std::map<int,int> > vec_pixel_id_store;
    vec_pixel_id_store.resize( nSegments );

    segIJ_W.clear();
    segIJ_P_ids.clear();
    segIJ_Q_ids.clear();
    segIJ_P_ids.resize( Idk.size() );
    segIJ_Q_ids.resize( Idk.size() );
    segIJ_W.resize( Idk.size() );
    segP2d.clear();
    segP2d.resize(nSegments);

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

          int idN_b = segImg[ posB ];
          if (idN_a >= idN_b) continue; // same segment

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

          // store the pixel 2d pair per segment, the weight, the ids of the pixel involved in the segment-pixel lists
          // how to get rid of the double pixel

          int storePos    = segIJMap[idN_a*nSegments + idN_b];
          
          if (storePos < 0) 
          {
            int weird = 0;
            if(displace <4)
              weird = 1;
            continue;
          }
          //int storeSegAPos = segIMap[idN_a];
          //int storeSegBPos = segIMap[idN_b];
          int storeSegAPos = idN_a;
          int storeSegBPos = idN_b;

//          std::map<int,int> &pixel_id_store = vec_pixel_id_store[storePos];

          //vec_pixel_id_store: stores the ids of the pixel already saved for the segment
          // a vector of maps pixid -> store id

          // does the pixel already exist?
          std::map<int,int> &pixel_id_storeA = vec_pixel_id_store[storeSegAPos];
          std::map<int,int> &pixel_id_storeB = vec_pixel_id_store[storeSegBPos];

          // py is from 0 to height
          // px from 0 to width
          // qx is from -1 to width+1
          // dy and dx are from -1..1
          // so px+qx+dy[displace]) is from -2 to width*2 + 2

          int smallPixId_p = (2+ px+qx+dy[displace]) * maxSmallHeight + py+qy+dx[displace] +2;
          int smallPixId_q = (2+ px+qx-dy[displace]) * maxSmallHeight + py+qy-dx[displace] +2;

          std::map<int,int>::iterator A_end(pixel_id_storeA.end()), B_end(pixel_id_storeB.end());
          // and if so where is it stored ?
          std::map<int,int>::iterator posP_A = pixel_id_storeA.find( smallPixId_p );
          std::map<int,int>::iterator posQ_A = pixel_id_storeA.find( smallPixId_q );
          std::map<int,int>::iterator posP_B = pixel_id_storeB.find( smallPixId_p );
          std::map<int,int>::iterator posQ_B = pixel_id_storeB.find( smallPixId_q );

          // the two points involved:
          // use a map per segment with the keys:
          // id is (px+qx+dy[displace])/2*2 = (px+qx+dy[displace])+ (py+qy-dx[displace])*width, width = 2*widthX+widthY
          P3 p2d = iK * P3( (px+qx+dy[displace])/2., (py+qy+dx[displace])/2., 1. );
          P3 q2d = iK * P3( (px+qx-dy[displace])/2., (py+qy-dx[displace])/2., 1. );
          // store these ,s.t. points are unique 
          // store ids in pair (i,j) vector

          ///////////////////////////////////////////////////////
          // segIJ_Q_ids[storePos] stores the ids of q_seg(a) and q_seg(b)
          // segIJ_P_ids[storePos] stores the ids of p_seg(a) and p_seg(b)
          ///////////////////////////////////////////
          int id_p2d_A(-1);
          int id_p2d_B(-1);

          if (posP_A != A_end)
            id_p2d_A = posP_A->second;
          else // append
          {
            segP2d[idN_a].push_back( p2d ); // p2d per seg
            // 4 ids! pa,pb,qa,qb
            id_p2d_A = segP2d[idN_a].size()-1;
            pixel_id_storeA[smallPixId_p] = id_p2d_A;
          }
          if (posP_B != B_end)
            id_p2d_B = posP_B->second;
          else // append
          {
            segP2d[idN_b].push_back( p2d ); // p2d per seg
            // 4 ids! pa,pb,qa,qb
            id_p2d_B = segP2d[idN_b].size()-1;
            pixel_id_storeB[smallPixId_p] = id_p2d_B;
          }
          segIJ_P_ids[storePos].push_back( std::pair<int,int> (id_p2d_A, id_p2d_B) ); // ids of the segments

          if (posQ_A != A_end)
            id_p2d_A = posQ_A->second;
          else // append
          {
            segP2d[idN_a].push_back( q2d ); // p2d per seg
            // 4 ids! pa,pb,qa,qb
            id_p2d_A = segP2d[idN_a].size()-1;
            pixel_id_storeA[smallPixId_q] = id_p2d_A;
          }
          if (posQ_B != B_end)
            id_p2d_B = posQ_B->second;
          else // append
          {
            segP2d[idN_b].push_back( q2d ); // p2d per seg
            // 4 ids! pa,pb,qa,qb
            id_p2d_B = segP2d[idN_b].size()-1;
            pixel_id_storeB[smallPixId_q] = id_p2d_B;
          }
          segIJ_Q_ids[storePos].push_back( std::pair<int,int> (id_p2d_A, id_p2d_B )); // ids of the segments
          segIJ_W[storePos].push_back( theWeight ); 
          ///////////////////////////////////////////////////////
          ///////////////////////////////////////////////////////
        }
      }
    }
  }

  /// segimage flipped so also width and height
  void initMapping( int nSegments_, const mxArray* edges_, Scalar* centers_, int* segImg, int width, int height )
  {
    nSegments = nSegments_;
    Idk.clear();
    Idl.clear();
    centersA.clear();

    // reserves not enough but i push_back, so ...
    Idk.reserve( nSegments );   // for each interaction a pair
    Idl.reserve( nSegments );    
    centersA.reserve( nSegments ); // 1 per segment

    segIJMap.clear();
    segIJMap.resize(nSegments*nSegments, -1);

    int ijNumber = 0;

    for (int i = 0; i < nSegments ;i++)
    {
//      Scalar* wts   = (Scalar*) mxGetPr( mxGetCell(weights_, i) );
      int nIds      = (int)     mxGetNumberOfElements( mxGetCell(edges_, i) )/5;
      Scalar* edge  = (Scalar*) mxGetPr( mxGetCell(edges_, i) );

      P3 centerA = P3( &centers_[ 3*i ] );
      centersA.push_back(centerA);

      for (int j =0; j < nIds ;j++)
      {
        int id = (int) (edge [5*j ])-1;

        if (id < i) continue;// can t we even break ?
        
        segIJMap[i*nSegments+id] = ijNumber++;

//        P3 centerB = P3( &centers_[ 3*id ] );
//        centersB.push_back(centerB);

        Idk.push_back(i);
        Idl.push_back(id);
      }}

    F01.resize(Idk.size());
    F00.resize(Idk.size());
    F10.resize(Idk.size());
    F11.resize(Idk.size());

    // now build per pixel pairs
    buildWeights( segImg, width, height,  nSegments);
  }

  /*
  /// segimage flipped so also width and height
  void prepareOwnWeights( int nSegments_, const mxArray* edges_, Scalar* centers_, int* segImg, int width, int height)
  {
    nSegments = nSegments_;
    weights.clear();
    edgesP.clear();
    edgesQ.clear();
    Idk.clear();
    Idl.clear();
    centersA.clear();
    centersB.clear();
    weightMap.clear();

    buildWeights( segImg, width, height,  nSegments);


    // reserves not enough but i push_back, so ...
    weights.reserve( nSegments );
    edgesP.reserve(  nSegments );
    edgesQ.reserve(  nSegments );
    Idk.reserve( nSegments );
    Idl.reserve( nSegments );
    centersA.reserve( nSegments );
    centersB.reserve( nSegments );
    weightMap.reserve( nSegments );
    weightMap.resize( nSegments );

    for (int i = 0; i < nSegments ;i++)
    {
//      Scalar* wts   = (Scalar*) mxGetPr( mxGetCell(weights_, i) );
      int nIds      = (int)     mxGetNumberOfElements( mxGetCell(edges_, i) )/5;
      Scalar* edge  = (Scalar*) mxGetPr( mxGetCell(edges_, i) );

//      P3 centerA = P3( centers_[ 2*i ], centers_[ 2*i+1], 1.0 );
      P3 centerA = P3( &centers_[ 3*i ] );

      for (int j =0; j < nIds ;j++)
      {
        int id = (int) (edge [5*j ])-1;

        if (id < i) continue;// can t we even break ?

        Scalar w = (edgeWeightMap[i])[id];
//        Scalar w = wts[j];

//        P3 centerB = P3( centers_[ 2*id ], centers_[ 2*id+1], 1.0 );
        P3 centerB = P3( &centers_[ 3*id ] );

#ifdef _Debug_out_
        printf("weight: %f", w);
#endif
#ifdef _Debug_out_
        printf("- id: %d", id);
#endif

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
*/
  size_t edges_num() {return Idk.size();};


  /// set up the 2d projections segP2d_0l, .. in order to 
  void compute_InitMappings( std::vector<int>& curr_ids )
  {
    n0.resize(nSegments);
    n1.resize(nSegments);

    segP2d_0l.resize(nSegments);
    segP2d_1l.resize(nSegments);
    segP2d_0r.resize(nSegments);
    segP2d_1r.resize(nSegments);


#pragma omp parallel for schedule (static)
    for (int segi = 0; segi < curr_ids.size(); segi++)
    {
      int idN_a = curr_ids[ segi ];
      std::vector< P3 >& p2ds = segP2d[segi];

      P3 norm0 = (*normals)[idN_a];
      M3 roti  = (*rot)[idN_a];
      P3 trai  = (*tra)[idN_a];

      P3 centerA ( centersA[segi] );
      P3 a_i = centerA * (1. / (norm0 | centerA));
      a_i -= roti*a_i - trai;

      Scalar dn_i =  1. / sqrt(norm0|norm0);
      P3 N00  = norm0*dn_i;
      P3 dN00 = roti * N00;

      n0 [segi] =  N00;
      n1 [segi] = dN00;

      segP2d_0l[segi].resize( p2ds.size() );
      segP2d_1l[segi].resize( p2ds.size() );
      segP2d_0r[segi].resize( p2ds.size() );
      segP2d_1r[segi].resize( p2ds.size() );

      std::vector<P2>& segP2d_0l_I = segP2d_0l[segi];
      std::vector<P2>& segP2d_1l_I = segP2d_1l[segi];
      std::vector<P2>& segP2d_0r_I = segP2d_0r[segi];
      std::vector<P2>& segP2d_1r_I = segP2d_1r[segi];

      // here 
      for (int i=0;i< p2ds.size();i++)
      {
        P3 p2d = p2ds[i];

        Scalar dp_i = 1. / max(norm0   | p2d, 0.00000001 );
        P3 p_i   = p2d * dp_i;
        //P3 pi_t1 = roti*(p_i-a_i)+trai+a_i; 
        P3 pi_t1 = roti*p_i+a_i; 

        P2 pi_t0_l; P2 pi_t0_r;
        P2 pi_t1_l; P2 pi_t1_r;

//        if (simpleProjection)
//        {
        get2dprojections( p_i,    pi_t0_l,  pi_t0_r );
        get2dprojections( pi_t1,  pi_t1_l,  pi_t1_r );
//      }
        segP2d_0l_I[i] = pi_t0_l;
        segP2d_1l_I[i] = pi_t1_l;
        segP2d_0r_I[i] = pi_t0_r;
        segP2d_1r_I[i] = pi_t1_r;
      }
    }
  }

  // not used ?
  void update( std::vector<int>& updateSolution, int trial )
  {
#pragma omp parallel for schedule (static)
    for (int segi = 0; segi < nSegments; segi++)
    {
      if (updateSolution[segi] == trial)
      {
        std::copy( t_segP2d_0l[segi].begin(), t_segP2d_0l[segi].end(), segP2d_0l[segi].begin() );
        std::copy( t_segP2d_1l[segi].begin(), t_segP2d_1l[segi].end(), segP2d_1l[segi].begin() );
        std::copy( t_segP2d_0r[segi].begin(), t_segP2d_0r[segi].end(), segP2d_0r[segi].begin() );
        std::copy( t_segP2d_1r[segi].begin(), t_segP2d_1r[segi].end(), segP2d_1r[segi].begin() );

        n0 [segi] = t_n0;
        n1 [segi] = t_n1;
      }
    }  
  }

  /// set up the 2d projections segP2d_0l, .. in order to 
  void compute_InitTrial( int trial )
  {
    n0.resize(nSegments);
    n1.resize(nSegments);

    t_segP2d_0l.resize(nSegments);
    t_segP2d_1l.resize(nSegments);
    t_segP2d_0r.resize(nSegments);
    t_segP2d_1r.resize(nSegments);

    int idN_a = trial;
    P3 norm0 = (*normals)[idN_a];
    M3 roti  = (*rot)[idN_a];
    P3 trai  = (*tra)[idN_a];
    Scalar dn_i =  1. / sqrt(norm0|norm0);

    P3 N00  = norm0*dn_i;
    P3 dN00 = roti * N00;
    t_n0 =  N00;
    t_n1 = dN00;

#pragma omp parallel for schedule (static)
    for (int segi = 0; segi < nSegments; segi++)
    {
      std::vector< P3 >& p2ds = segP2d[segi];
      P3 centerA ( centersA[segi] );

      P3 a_i = centerA * (1. / (norm0 | centerA));
      a_i -= roti*a_i - trai;

      t_segP2d_0l[segi].resize( p2ds.size() );
      t_segP2d_1l[segi].resize( p2ds.size() );
      t_segP2d_0r[segi].resize( p2ds.size() );
      t_segP2d_1r[segi].resize( p2ds.size() );

      std::vector<P2>& t_segP2d_0l_I = t_segP2d_0l[segi];
      std::vector<P2>& t_segP2d_1l_I = t_segP2d_1l[segi];
      std::vector<P2>& t_segP2d_0r_I = t_segP2d_0r[segi];
      std::vector<P2>& t_segP2d_1r_I = t_segP2d_1r[segi];

      // here 
      for (int i=0;i< p2ds.size();i++)
      {
        P3 p2d = p2ds[i];
        Scalar dp_i = 1. / max(norm0   | p2d, 0.00000001 );

        P3 p_i   = p2d * dp_i;
        P3 pi_t1 = roti*p_i+a_i; 

        P2 pi_t0_l; P2 pi_t0_r;
        P2 pi_t1_l; P2 pi_t1_r;

//        if (simpleProjection)
//        {
        get2dprojections( p_i,    pi_t0_l,  pi_t0_r );
        get2dprojections( pi_t1,  pi_t1_l,  pi_t1_r );
//        }
        t_segP2d_0l_I[i] = pi_t0_l;
        t_segP2d_1l_I[i] = pi_t1_l;
        t_segP2d_0r_I[i] = pi_t0_r;
        t_segP2d_1r_I[i] = pi_t1_r;
      }
    }
  }


  void compute_InitTrial( const std::vector<int> &trialId )
  {
    n0.resize(nSegments);
    n1.resize(nSegments);
    tv_n0.resize(nSegments);
    tv_n1.resize(nSegments);

    t_segP2d_0l.resize(nSegments);
    t_segP2d_1l.resize(nSegments);
    t_segP2d_0r.resize(nSegments);
    t_segP2d_1r.resize(nSegments);

#pragma omp parallel for schedule (static)
    for (int segi = 0; segi < nSegments; segi++)
    {
      int idN_a = trialId[segi];
      P3 norm0 = (*normals)[idN_a];
      M3 roti  = (*rot)[idN_a];
      P3 trai  = (*tra)[idN_a];
      Scalar dn_i =  1. / sqrt(norm0|norm0);

      P3 N00  = norm0*dn_i;
      P3 dN00 = roti * N00;
      t_n0 =  N00;
      t_n1 = dN00;

      tv_n0[segi] =  N00;
      tv_n1[segi] = dN00;

      std::vector< P3 >& p2ds = segP2d[segi];
      P3 centerA ( centersA[segi] );

      P3 a_i = centerA * (1. / (norm0 | centerA));
      a_i -= roti*a_i - trai;

      t_segP2d_0l[segi].resize( p2ds.size() );
      t_segP2d_1l[segi].resize( p2ds.size() );
      t_segP2d_0r[segi].resize( p2ds.size() );
      t_segP2d_1r[segi].resize( p2ds.size() );

      std::vector<P2>& t_segP2d_0l_I = t_segP2d_0l[segi];
      std::vector<P2>& t_segP2d_1l_I = t_segP2d_1l[segi];
      std::vector<P2>& t_segP2d_0r_I = t_segP2d_0r[segi];
      std::vector<P2>& t_segP2d_1r_I = t_segP2d_1r[segi];

      // here 
      for (int i=0;i< p2ds.size();i++)
      {
        P3 p2d = p2ds[i];
        Scalar dp_i = 1. / max(norm0   | p2d, 0.00000001 );

        P3 p_i   = p2d * dp_i;
        P3 pi_t1 = roti*p_i+a_i; 

        P2 pi_t0_l; P2 pi_t0_r;
        P2 pi_t1_l; P2 pi_t1_r;

//        if (simpleProjection)
//        {
        get2dprojections( p_i,    pi_t0_l,  pi_t0_r );
        get2dprojections( pi_t1,  pi_t1_l,  pi_t1_r );
//        }
        t_segP2d_0l_I[i] = pi_t0_l;
        t_segP2d_1l_I[i] = pi_t1_l;
        t_segP2d_0r_I[i] = pi_t0_r;
        t_segP2d_1r_I[i] = pi_t1_r;
      }
    }
  }
/*
  Scalar compute_score_combiDepth_orig( int id, std::vector<int>& dummy )
  {
    compute_InitTrial( id );
    Scalar score =0;

#pragma omp parallel for schedule (static)
    for( int Sij = 0; Sij < Idk.size(); Sij++)
    {
      Scalar f00_it(0.);
      Scalar f01_it(0.);
      Scalar f10_it(0.);
      Scalar f11_it(0.);

      std::vector<std::pair<int,int> >& pids = segIJ_P_ids[Sij];
      std::vector<std::pair<int,int> >& qids = segIJ_Q_ids[Sij];
      std::vector< Scalar >& sij_weights     = segIJ_W[Sij];

      const int segI = Idk[Sij];
      const int segJ = Idl[Sij];

      const std::vector< P2 >& segP2d_0l_I = segP2d_0l[segI];
      const std::vector< P2 >& segP2d_0l_J = segP2d_0l[segJ];
      const std::vector< P2 >& segP2d_0r_I = segP2d_0r[segI];
      const std::vector< P2 >& segP2d_0r_J = segP2d_0r[segJ];
      const std::vector< P2 >& segP2d_1l_I = segP2d_1l[segI];
      const std::vector< P2 >& segP2d_1l_J = segP2d_1l[segJ];
      const std::vector< P2 >& segP2d_1r_I = segP2d_1r[segI];
      const std::vector< P2 >& segP2d_1r_J = segP2d_1r[segJ];

      const std::vector< P2 >& t_segP2d_0l_I = t_segP2d_0l[segI];
      const std::vector< P2 >& t_segP2d_0l_J = t_segP2d_0l[segJ];
      const std::vector< P2 >& t_segP2d_0r_I = t_segP2d_0r[segI];
      const std::vector< P2 >& t_segP2d_0r_J = t_segP2d_0r[segJ];
      const std::vector< P2 >& t_segP2d_1l_I = t_segP2d_1l[segI];
      const std::vector< P2 >& t_segP2d_1l_J = t_segP2d_1l[segJ];
      const std::vector< P2 >& t_segP2d_1r_I = t_segP2d_1r[segI];
      const std::vector< P2 >& t_segP2d_1r_J = t_segP2d_1r[segJ];

      //P3 dN00r = roti * norm0*dn_i   - rotj * norm1*dn_j;
      //P3 dN01r = roti * norm0*dn_i   - rotk * normFix*dn_k;
      //P3 dN10r = rotk * normFix*dn_k - rotj * norm1*dn_j;
      //P3 dN11r(0.,0.,0.) ;//= rotk * normFix*dn_k - rotk * normFix*dn_k;

            tv_n0[segi] =  N00;
      tv_n1[segi] = dN00;


      P3 dN00r = gamma*(n1[ segI ] - n1[ segJ ]);
      P3 dN00  = gamma*(n0[ segI ] - n0[ segJ ]);
      P3 dN01r = gamma*(n1[ segI ] - tv_n1[segJ]);
      P3 dN01  = gamma*(n0[ segI ] - tv_n0[segJ]);
      P3 dN10r = gamma*(tv_n1[segI] - n1[ segJ ]);
      P3 dN10  = gamma*(tv_n0[segJ] - n0[ segJ ]);
      P3 dN11r = gamma*(tv_n1[segI] - tv_n1[segJ]);
      P3 dN11  = gamma*(tv_n0[segJ] - tv_n0[segJ]);

      for (int pq_ids = 0; pq_ids < pids.size(); pq_ids++ )
      {
        P2 f00, f00_t, f00_rt, f11, f11_t, f11_rt;
        P2 g00, g00_t, g00_rt, g11, g11_t, g11_rt;
        P2 h00, h00_t, h00_rt, h11, h11_t, h11_rt;
        P2 e00, e00_t, e00_rt, e11, e11_t, e11_rt;

          P2 pi_t0_l = segP2d_0l_I[ pids[pq_ids].first ];
          P2 pi_t1_l = segP2d_1l_I[ pids[pq_ids].first ];
          P2 pi_t0_r = segP2d_0r_I[ pids[pq_ids].first ];
          P2 pi_t1_r = segP2d_1r_I[ pids[pq_ids].first ];

          P2 qi_t0_l = segP2d_0l_I[ qids[pq_ids].first ];
          P2 qi_t1_l = segP2d_1l_I[ qids[pq_ids].first ];
          P2 qi_t0_r = segP2d_0r_I[ qids[pq_ids].first ];
          P2 qi_t1_r = segP2d_1r_I[ qids[pq_ids].first ];

          P2 pj_t0_l = segP2d_0l_J[ pids[pq_ids].second];
          P2 pj_t1_l = segP2d_1l_J[ pids[pq_ids].second];
          P2 pj_t0_r = segP2d_0r_J[ pids[pq_ids].second];
          P2 pj_t1_r = segP2d_1r_J[ pids[pq_ids].second];

          P2 qj_t0_l = segP2d_0l_J[ qids[pq_ids].second];
          P2 qj_t1_l = segP2d_1l_J[ qids[pq_ids].second];
          P2 qj_t0_r = segP2d_0r_J[ qids[pq_ids].second];
          P2 qj_t1_r = segP2d_1r_J[ qids[pq_ids].second];


          P2 pk_t0_l  = t_segP2d_0l_I[ pids[pq_ids].first ];// same as for second
          P2 pkA_t1_l = t_segP2d_1l_I[ pids[pq_ids].first ];
          P2 pk_t0_r  = t_segP2d_0r_I[ pids[pq_ids].first ];
          P2 pkA_t1_r = t_segP2d_1r_I[ pids[pq_ids].first ];
          P2 pkB_t1_l = t_segP2d_1l_J[ pids[pq_ids].second];
          P2 pkB_t1_r = t_segP2d_1r_J[ pids[pq_ids].second];

          P2 qk_t0_l  = t_segP2d_0l_I[ qids[pq_ids].first ];
          P2 qkA_t1_l = t_segP2d_1l_I[ qids[pq_ids].first ];
          P2 qk_t0_r  = t_segP2d_0r_I[ qids[pq_ids].first ];
          P2 qkA_t1_r = t_segP2d_1r_I[ qids[pq_ids].first ];
          P2 qkB_t1_l = t_segP2d_1l_J[ qids[pq_ids].second];
          P2 qkB_t1_r = t_segP2d_1r_J[ qids[pq_ids].second];


          // four times pi qj pk qk combinations
          get2dDiff_fast( pi_t0_l, pi_t0_r, pi_t1_l, pi_t1_r,
                          pj_t0_l, pj_t0_r, pj_t1_l, pj_t1_r,
                          f00, f00_t, f00_rt );
          get2dDiff_fast( qi_t0_l, qi_t0_r, qi_t1_l, qi_t1_r,
                          qj_t0_l, qj_t0_r, qj_t1_l, qj_t1_r,
                          f11, f11_t, f11_rt );

          get2dDiff_fast( pi_t0_l, pi_t0_r, pi_t1_l, pi_t1_r,
                          pk_t0_l, pk_t0_r, pkB_t1_l, pkB_t1_r,
                          g00, g00_t, g00_rt );
          get2dDiff_fast( qi_t0_l, qi_t0_r, qi_t1_l, qi_t1_r,
                          qk_t0_l, qk_t0_r, qkB_t1_l, qkB_t1_r,
                          g11, g11_t, g11_rt );

          get2dDiff_fast( pk_t0_l, pk_t0_r, pkA_t1_l, pkA_t1_r,
                          pj_t0_l, pj_t0_r, pj_t1_l, pj_t1_r,
                          h00, h00_t, h00_rt );
          get2dDiff_fast( qk_t0_l, qk_t0_r, qkA_t1_l, qkA_t1_r,
                          qj_t0_l, qj_t0_r, qj_t1_l, qj_t1_r,
                          h11, h11_t, h11_rt );

          get2dDiff_fast( pk_t0_l, pk_t0_r, pkA_t1_l, pkA_t1_r,
                          pk_t0_l, pk_t0_r, pkB_t1_l, pkB_t1_r,
                          e00, e00_t, e00_rt );
          get2dDiff_fast( qk_t0_l, qk_t0_r, qkA_t1_l, qkA_t1_r,
                          qk_t0_l, qk_t0_r, qkB_t1_l, qkB_t1_r,
                          e11, e11_t, e11_rt );

#ifdef __combi2dSmoothMotion__
          // some Hack combine both motions 
          P4 f00_m ( f00_t[0], f00_t[1], f00_rt[0], f00_rt[1] );
          P4 f11_m ( f11_t[0], f11_t[1], f11_rt[0], f11_rt[1] );
          P4 g00_m ( g00_t[0], g00_t[1], g00_rt[0], g00_rt[1] );
          P4 g11_m ( g11_t[0], g11_t[1], g11_rt[0], g11_rt[1] );
          P4 h00_m ( h00_t[0], h00_t[1], h00_rt[0], h00_rt[1] );
          P4 h11_m ( h11_t[0], h11_t[1], h11_rt[0], h11_rt[1] );
          P4 e00_m ( e00_t[0], e00_t[1], e00_rt[0], e00_rt[1] );
          P4 e11_m ( e11_t[0], e11_t[1], e11_rt[0], e11_rt[1] );

          Scalar f00_M = sqrt ( ((f00_m|f00_m) + (f11_m|f11_m) + (f00_m|f11_m) + (dN00r|dN00r)) * Scalar(3.3) + epsilon);
          Scalar f10_M = sqrt ( ((h00_m|h00_m) + (h11_m|h11_m) + (h00_m|h11_m) + (dN10r|dN10r)) * Scalar(3.3) + epsilon);
          Scalar f01_M = sqrt ( ((g00_m|g00_m) + (g11_m|g11_m) + (g00_m|g11_m) + (dN01r|dN01r)) * Scalar(3.3) + epsilon);
          Scalar f11_M = sqrt ( ((e00_m|e00_m) + (e11_m|e11_m) + (e00_m|e11_m) + (dN11r|dN11r)) * Scalar(3.3) + epsilon);
#endif

          Scalar f00_D = sqrt ( ( (f00|f00) + (f11|f11) + (f00|f11) + (dN00|dN00) ) * Scalar(3.3) + epsilon);
          Scalar f10_D = sqrt ( ( (h00|h00) + (h11|h11) + (h00|h11) + (dN10|dN10) ) * Scalar(3.3) + epsilon);
          Scalar f01_D = sqrt ( ( (g00|g00) + (g11|g11) + (g00|g11) + (dN01|dN01) ) * Scalar(3.3) + epsilon);
          Scalar f11_D = sqrt ( ( (e00|e00) + (e11|e11) + (e00|e11) + (dN11|dN11) ) * Scalar(3.3) + epsilon);
//          Scalar f11_D = sqrt (epsilon);

#ifndef __do_depth_only__
          f00_it += (min(f00_D, depthJump) + rotWeight * min(f00_M, rotJump) ) * sij_weights[pq_ids];
          f10_it += (min(f10_D, depthJump) + rotWeight * min(f10_M, rotJump) ) * sij_weights[pq_ids];
          f01_it += (min(f01_D, depthJump) + rotWeight * min(f01_M, rotJump) ) * sij_weights[pq_ids];
          f11_it += (min(f11_D, depthJump) + rotWeight * min(f11_M, rotJump) ) * sij_weights[pq_ids];
#else
          f00_it += min(f00_D, depthJump) * sij_weights[pq_ids];
          f10_it += min(f10_D, depthJump) * sij_weights[pq_ids];
          f01_it += min(f01_D, depthJump) * sij_weights[pq_ids];
          f11_it += min(f11_D, depthJump) * sij_weights[pq_ids];
#endif
        }

        F00[Sij] = f00_it;
        F01[Sij] = f01_it;
        F10[Sij] = f10_it;
        F11[Sij] = f11_it;

//        weightMap[ Idk[Sij] ].insert( std::pair<int,Scalar> ( Idk[Sij], f10_it + f01_it - f11_it - f00_it ) );
        score += max(max(max(F00[Sij], F01[Sij]), F10[Sij]), F11[Sij]);

        if ( score != score || std::numeric_limits<Scalar>::infinity()==score )
        {
          printf ("WHAT WHY %f \n", score);
          mexEvalString("drawnow");
          int breakhere = 0;
        }  
    }
    return score;
  }
*/
  // try to apply the sqrt or thresholdig after the accumulation
  Scalar compute_score_combiDepth_sqrt( int id, std::vector<int>& dummy )
  {
    compute_InitTrial( id );
    Scalar score =0;

#pragma omp parallel for schedule (static)
    for( int Sij = 0; Sij < Idk.size(); Sij++)
    {
      Scalar f00_it(0.);
      Scalar f01_it(0.);
      Scalar f10_it(0.);
      Scalar f11_it(0.);
 
      Scalar f00_it_M(0.);
      Scalar f01_it_M(0.);
      Scalar f10_it_M(0.);
      Scalar f11_it_M(0.);

      std::vector<std::pair<int,int> >& pids = segIJ_P_ids[Sij];
      std::vector<std::pair<int,int> >& qids = segIJ_Q_ids[Sij];
      std::vector< Scalar >& sij_weights     = segIJ_W[Sij];

      const int segI = Idk[Sij];
      const int segJ = Idl[Sij];

      const std::vector< P2 >& segP2d_0l_I = segP2d_0l[segI];
      const std::vector< P2 >& segP2d_0l_J = segP2d_0l[segJ];
      const std::vector< P2 >& segP2d_0r_I = segP2d_0r[segI];
      const std::vector< P2 >& segP2d_0r_J = segP2d_0r[segJ];
      const std::vector< P2 >& segP2d_1l_I = segP2d_1l[segI];
      const std::vector< P2 >& segP2d_1l_J = segP2d_1l[segJ];
      const std::vector< P2 >& segP2d_1r_I = segP2d_1r[segI];
      const std::vector< P2 >& segP2d_1r_J = segP2d_1r[segJ];

      const std::vector< P2 >& t_segP2d_0l_I = t_segP2d_0l[segI];
      const std::vector< P2 >& t_segP2d_0l_J = t_segP2d_0l[segJ];
      const std::vector< P2 >& t_segP2d_0r_I = t_segP2d_0r[segI];
      const std::vector< P2 >& t_segP2d_0r_J = t_segP2d_0r[segJ];
      const std::vector< P2 >& t_segP2d_1l_I = t_segP2d_1l[segI];
      const std::vector< P2 >& t_segP2d_1l_J = t_segP2d_1l[segJ];
      const std::vector< P2 >& t_segP2d_1r_I = t_segP2d_1r[segI];
      const std::vector< P2 >& t_segP2d_1r_J = t_segP2d_1r[segJ];

      //P3 dN00r = roti * norm0*dn_i   - rotj * norm1*dn_j;
      //P3 dN01r = roti * norm0*dn_i   - rotk * normFix*dn_k;
      //P3 dN10r = rotk * normFix*dn_k - rotj * norm1*dn_j;
      //P3 dN11r(0.,0.,0.) ;//= rotk * normFix*dn_k - rotk * normFix*dn_k;

      P3 dN00r = gamma*(n1[ segI ] - n1[ segJ ]);
      P3 dN00  = gamma*(n0[ segI ] - n0[ segJ ]);
      P3 dN01r = gamma*(n1[ segI ] - t_n1);
      P3 dN01  = gamma*(n0[ segI ] - t_n0);
      P3 dN10r = gamma*(t_n1 - n1[ segJ ]);
      P3 dN10  = gamma*(t_n0 - n0[ segJ ]);
      P3 dN11r(0,0,0);
      P3 dN11 (0,0,0);

      for (int pq_ids = 0; pq_ids < pids.size(); pq_ids++ )
      {
        P2 f00, f00_t, f00_rt, f11, f11_t, f11_rt;
        P2 g00, g00_t, g00_rt, g11, g11_t, g11_rt;
        P2 h00, h00_t, h00_rt, h11, h11_t, h11_rt;
        P2 e00, e00_t, e00_rt, e11, e11_t, e11_rt;

          P2 pi_t0_l = segP2d_0l_I[ pids[pq_ids].first ];
          P2 pi_t1_l = segP2d_1l_I[ pids[pq_ids].first ];
          P2 pi_t0_r = segP2d_0r_I[ pids[pq_ids].first ];
          P2 pi_t1_r = segP2d_1r_I[ pids[pq_ids].first ];

          P2 qi_t0_l = segP2d_0l_I[ qids[pq_ids].first ];
          P2 qi_t1_l = segP2d_1l_I[ qids[pq_ids].first ];
          P2 qi_t0_r = segP2d_0r_I[ qids[pq_ids].first ];
          P2 qi_t1_r = segP2d_1r_I[ qids[pq_ids].first ];

          P2 pj_t0_l = segP2d_0l_J[ pids[pq_ids].second];
          P2 pj_t1_l = segP2d_1l_J[ pids[pq_ids].second];
          P2 pj_t0_r = segP2d_0r_J[ pids[pq_ids].second];
          P2 pj_t1_r = segP2d_1r_J[ pids[pq_ids].second];

          P2 qj_t0_l = segP2d_0l_J[ qids[pq_ids].second];
          P2 qj_t1_l = segP2d_1l_J[ qids[pq_ids].second];
          P2 qj_t0_r = segP2d_0r_J[ qids[pq_ids].second];
          P2 qj_t1_r = segP2d_1r_J[ qids[pq_ids].second];


          P2 pk_t0_l  = t_segP2d_0l_I[ pids[pq_ids].first ];// same as for second
          P2 pkA_t1_l = t_segP2d_1l_I[ pids[pq_ids].first ];
          P2 pk_t0_r  = t_segP2d_0r_I[ pids[pq_ids].first ];
          P2 pkA_t1_r = t_segP2d_1r_I[ pids[pq_ids].first ];
          P2 pkB_t1_l = t_segP2d_1l_J[ pids[pq_ids].second];
          P2 pkB_t1_r = t_segP2d_1r_J[ pids[pq_ids].second];

          P2 qk_t0_l  = t_segP2d_0l_I[ qids[pq_ids].first ];
          P2 qkA_t1_l = t_segP2d_1l_I[ qids[pq_ids].first ];
          P2 qk_t0_r  = t_segP2d_0r_I[ qids[pq_ids].first ];
          P2 qkA_t1_r = t_segP2d_1r_I[ qids[pq_ids].first ];
          P2 qkB_t1_l = t_segP2d_1l_J[ qids[pq_ids].second];
          P2 qkB_t1_r = t_segP2d_1r_J[ qids[pq_ids].second];


          // four times pi qj pk qk combinations
          get2dDiff_fast( pi_t0_l, pi_t0_r, pi_t1_l, pi_t1_r,
                          pj_t0_l, pj_t0_r, pj_t1_l, pj_t1_r,
                          f00, f00_t, f00_rt );
          get2dDiff_fast( qi_t0_l, qi_t0_r, qi_t1_l, qi_t1_r,
                          qj_t0_l, qj_t0_r, qj_t1_l, qj_t1_r,
                          f11, f11_t, f11_rt );

          get2dDiff_fast( pi_t0_l, pi_t0_r, pi_t1_l, pi_t1_r,
                          pk_t0_l, pk_t0_r, pkB_t1_l, pkB_t1_r,
                          g00, g00_t, g00_rt );
          get2dDiff_fast( qi_t0_l, qi_t0_r, qi_t1_l, qi_t1_r,
                          qk_t0_l, qk_t0_r, qkB_t1_l, qkB_t1_r,
                          g11, g11_t, g11_rt );

          get2dDiff_fast( pk_t0_l, pk_t0_r, pkA_t1_l, pkA_t1_r,
                          pj_t0_l, pj_t0_r, pj_t1_l, pj_t1_r,
                          h00, h00_t, h00_rt );
          get2dDiff_fast( qk_t0_l, qk_t0_r, qkA_t1_l, qkA_t1_r,
                          qj_t0_l, qj_t0_r, qj_t1_l, qj_t1_r,
                          h11, h11_t, h11_rt );

          get2dDiff_fast( pk_t0_l, pk_t0_r, pkA_t1_l, pkA_t1_r,
                          pk_t0_l, pk_t0_r, pkB_t1_l, pkB_t1_r,
                          e00, e00_t, e00_rt );
          get2dDiff_fast( qk_t0_l, qk_t0_r, qkA_t1_l, qkA_t1_r,
                          qk_t0_l, qk_t0_r, qkB_t1_l, qkB_t1_r,
                          e11, e11_t, e11_rt );

#ifdef __combi2dSmoothMotion__
          // some Hack combine both motions 
          P4 f00_m ( f00_t[0], f00_t[1], f00_rt[0], f00_rt[1] );
          P4 f11_m ( f11_t[0], f11_t[1], f11_rt[0], f11_rt[1] );
          P4 g00_m ( g00_t[0], g00_t[1], g00_rt[0], g00_rt[1] );
          P4 g11_m ( g11_t[0], g11_t[1], g11_rt[0], g11_rt[1] );
          P4 h00_m ( h00_t[0], h00_t[1], h00_rt[0], h00_rt[1] );
          P4 h11_m ( h11_t[0], h11_t[1], h11_rt[0], h11_rt[1] );
          P4 e00_m ( e00_t[0], e00_t[1], e00_rt[0], e00_rt[1] );
          P4 e11_m ( e11_t[0], e11_t[1], e11_rt[0], e11_rt[1] );

          Scalar f00_M = ( ((f00_m|f00_m) + (f11_m|f11_m) + (f00_m|f11_m) + (dN00r|dN00r)) * Scalar(3.3) + epsilon);
          Scalar f10_M = ( ((h00_m|h00_m) + (h11_m|h11_m) + (h00_m|h11_m) + (dN10r|dN10r)) * Scalar(3.3) + epsilon);
          Scalar f01_M = ( ((g00_m|g00_m) + (g11_m|g11_m) + (g00_m|g11_m) + (dN01r|dN01r)) * Scalar(3.3) + epsilon);
          Scalar f11_M = ( ((e00_m|e00_m) + (e11_m|e11_m) + (e00_m|e11_m) ) * Scalar(3.3) + epsilon);
#endif

          Scalar f00_D = ( ( (f00|f00) + (f11|f11) + (f00|f11) + (dN00|dN00) ) * Scalar(3.3) + epsilon);
          Scalar f10_D = ( ( (h00|h00) + (h11|h11) + (h00|h11) + (dN10|dN10) ) * Scalar(3.3) + epsilon);
          Scalar f01_D = ( ( (g00|g00) + (g11|g11) + (g00|g11) + (dN01|dN01) ) * Scalar(3.3) + epsilon);
          Scalar f11_D = (epsilon);

          Scalar ww2 = sij_weights[pq_ids] * sij_weights[pq_ids];

#ifndef __do_depth_only__

          f00_it += min(f00_D, depthJump*depthJump) * ww2;
          f10_it += min(f10_D, depthJump*depthJump) * ww2;
          f01_it += min(f01_D, depthJump*depthJump) * ww2;
          f11_it += min(f11_D, depthJump*depthJump) * ww2;

          f00_it_M += min(f00_M, rotJump*rotJump)   * ww2;
          f10_it_M += min(f10_M, rotJump*rotJump)   * ww2;
          f01_it_M += min(f01_M, rotJump*rotJump)   * ww2;
          f11_it_M += min(f11_M, rotJump*rotJump)   * ww2;
#else
          f00_it += min(f00_D, depthJump*depthJump) * ww2;
          f10_it += min(f10_D, depthJump*depthJump) * ww2;
          f01_it += min(f01_D, depthJump*depthJump) * ww2;
          f11_it += min(f11_D, depthJump*depthJump) * ww2;
#endif
      }

      F00[Sij] = sqrt(f00_it) + rotWeight * sqrt(f00_it_M);
      F01[Sij] = sqrt(f01_it) + rotWeight * sqrt(f01_it_M);
      F10[Sij] = sqrt(f10_it) + rotWeight * sqrt(f10_it_M);
      F11[Sij] = sqrt(f11_it) + rotWeight * sqrt(f11_it_M);

      weightMap[ Idk[Sij] ].insert( std::pair<int,Scalar> ( Idk[Sij], F10[Sij] + F01[Sij] - F11[Sij] - F00[Sij] ) );
      score += max(max(max(F00[Sij], F01[Sij]), F10[Sij]), F11[Sij]);

      if ( score != score || std::numeric_limits<Scalar>::infinity()==score )
      {
        printf ("WHAT WHY %f \n", score);
        mexEvalString("drawnow");
        int breakhere = 0;
      }  
    }
    return score;
  }

  // try to apply the thresholdig after the accumulation
  Scalar compute_score_combiDepth( int id, std::vector<int>& dummy )
  {
    compute_InitTrial( id );
    Scalar score =0;

#pragma omp parallel for schedule (static)
    for( int Sij = 0; Sij < Idk.size(); Sij++)
    {
      Scalar f00_it(0.);
      Scalar f01_it(0.);
      Scalar f10_it(0.);
      Scalar f11_it(0.);
 
      Scalar f00_it_M(0.);
      Scalar f01_it_M(0.);
      Scalar f10_it_M(0.);
      Scalar f11_it_M(0.);

      std::vector<std::pair<int,int> >& pids = segIJ_P_ids[Sij];
      std::vector<std::pair<int,int> >& qids = segIJ_Q_ids[Sij];
      std::vector< Scalar >& sij_weights     = segIJ_W[Sij];

      const int segI = Idk[Sij];
      const int segJ = Idl[Sij];

      const std::vector< P2 >& segP2d_0l_I = segP2d_0l[segI];
      const std::vector< P2 >& segP2d_0l_J = segP2d_0l[segJ];
      const std::vector< P2 >& segP2d_0r_I = segP2d_0r[segI];
      const std::vector< P2 >& segP2d_0r_J = segP2d_0r[segJ];
      const std::vector< P2 >& segP2d_1l_I = segP2d_1l[segI];
      const std::vector< P2 >& segP2d_1l_J = segP2d_1l[segJ];
      const std::vector< P2 >& segP2d_1r_I = segP2d_1r[segI];
      const std::vector< P2 >& segP2d_1r_J = segP2d_1r[segJ];

      const std::vector< P2 >& t_segP2d_0l_I = t_segP2d_0l[segI];
      const std::vector< P2 >& t_segP2d_0l_J = t_segP2d_0l[segJ];
      const std::vector< P2 >& t_segP2d_0r_I = t_segP2d_0r[segI];
      const std::vector< P2 >& t_segP2d_0r_J = t_segP2d_0r[segJ];
      const std::vector< P2 >& t_segP2d_1l_I = t_segP2d_1l[segI];
      const std::vector< P2 >& t_segP2d_1l_J = t_segP2d_1l[segJ];
      const std::vector< P2 >& t_segP2d_1r_I = t_segP2d_1r[segI];
      const std::vector< P2 >& t_segP2d_1r_J = t_segP2d_1r[segJ];

      //P3 dN00r = roti * norm0*dn_i   - rotj * norm1*dn_j;
      //P3 dN01r = roti * norm0*dn_i   - rotk * normFix*dn_k;
      //P3 dN10r = rotk * normFix*dn_k - rotj * norm1*dn_j;
      //P3 dN11r(0.,0.,0.) ;//= rotk * normFix*dn_k - rotk * normFix*dn_k;

      P3 dN00r = gamma*(n1[ segI ] - n1[ segJ ]);
      P3 dN00  = gamma*(n0[ segI ] - n0[ segJ ]);
      P3 dN01r = gamma*(n1[ segI ] - t_n1);
      P3 dN01  = gamma*(n0[ segI ] - t_n0);
      P3 dN10r = gamma*(t_n1 - n1[ segJ ]);
      P3 dN10  = gamma*(t_n0 - n0[ segJ ]);
      P3 dN11r(0,0,0);
      P3 dN11 (0,0,0);

      Scalar sumW(0);

      for (int pq_ids = 0; pq_ids < pids.size(); pq_ids++ )
      {
        P2 f00, f00_t, f00_rt, f11, f11_t, f11_rt;
        P2 g00, g00_t, g00_rt, g11, g11_t, g11_rt;
        P2 h00, h00_t, h00_rt, h11, h11_t, h11_rt;
        P2 e00, e00_t, e00_rt, e11, e11_t, e11_rt;

          P2 pi_t0_l = segP2d_0l_I[ pids[pq_ids].first ];
          P2 pi_t1_l = segP2d_1l_I[ pids[pq_ids].first ];
          P2 pi_t0_r = segP2d_0r_I[ pids[pq_ids].first ];
          P2 pi_t1_r = segP2d_1r_I[ pids[pq_ids].first ];

          P2 qi_t0_l = segP2d_0l_I[ qids[pq_ids].first ];
          P2 qi_t1_l = segP2d_1l_I[ qids[pq_ids].first ];
          P2 qi_t0_r = segP2d_0r_I[ qids[pq_ids].first ];
          P2 qi_t1_r = segP2d_1r_I[ qids[pq_ids].first ];

          P2 pj_t0_l = segP2d_0l_J[ pids[pq_ids].second];
          P2 pj_t1_l = segP2d_1l_J[ pids[pq_ids].second];
          P2 pj_t0_r = segP2d_0r_J[ pids[pq_ids].second];
          P2 pj_t1_r = segP2d_1r_J[ pids[pq_ids].second];

          P2 qj_t0_l = segP2d_0l_J[ qids[pq_ids].second];
          P2 qj_t1_l = segP2d_1l_J[ qids[pq_ids].second];
          P2 qj_t0_r = segP2d_0r_J[ qids[pq_ids].second];
          P2 qj_t1_r = segP2d_1r_J[ qids[pq_ids].second];


          P2 pk_t0_l  = t_segP2d_0l_I[ pids[pq_ids].first ];// same as for second
          P2 pkA_t1_l = t_segP2d_1l_I[ pids[pq_ids].first ];
          P2 pk_t0_r  = t_segP2d_0r_I[ pids[pq_ids].first ];
          P2 pkA_t1_r = t_segP2d_1r_I[ pids[pq_ids].first ];
          P2 pkB_t1_l = t_segP2d_1l_J[ pids[pq_ids].second];
          P2 pkB_t1_r = t_segP2d_1r_J[ pids[pq_ids].second];

          P2 qk_t0_l  = t_segP2d_0l_I[ qids[pq_ids].first ];
          P2 qkA_t1_l = t_segP2d_1l_I[ qids[pq_ids].first ];
          P2 qk_t0_r  = t_segP2d_0r_I[ qids[pq_ids].first ];
          P2 qkA_t1_r = t_segP2d_1r_I[ qids[pq_ids].first ];
          P2 qkB_t1_l = t_segP2d_1l_J[ qids[pq_ids].second];
          P2 qkB_t1_r = t_segP2d_1r_J[ qids[pq_ids].second];


          // four times pi qj pk qk combinations
          get2dDiff_fast( pi_t0_l, pi_t0_r, pi_t1_l, pi_t1_r,
                          pj_t0_l, pj_t0_r, pj_t1_l, pj_t1_r,
                          f00, f00_t, f00_rt );
          get2dDiff_fast( qi_t0_l, qi_t0_r, qi_t1_l, qi_t1_r,
                          qj_t0_l, qj_t0_r, qj_t1_l, qj_t1_r,
                          f11, f11_t, f11_rt );

          get2dDiff_fast( pi_t0_l, pi_t0_r, pi_t1_l, pi_t1_r,
                          pk_t0_l, pk_t0_r, pkB_t1_l, pkB_t1_r,
                          g00, g00_t, g00_rt );
          get2dDiff_fast( qi_t0_l, qi_t0_r, qi_t1_l, qi_t1_r,
                          qk_t0_l, qk_t0_r, qkB_t1_l, qkB_t1_r,
                          g11, g11_t, g11_rt );

          get2dDiff_fast( pk_t0_l, pk_t0_r, pkA_t1_l, pkA_t1_r,
                          pj_t0_l, pj_t0_r, pj_t1_l, pj_t1_r,
                          h00, h00_t, h00_rt );
          get2dDiff_fast( qk_t0_l, qk_t0_r, qkA_t1_l, qkA_t1_r,
                          qj_t0_l, qj_t0_r, qj_t1_l, qj_t1_r,
                          h11, h11_t, h11_rt );

          get2dDiff_fast( pk_t0_l, pk_t0_r, pkA_t1_l, pkA_t1_r,
                          pk_t0_l, pk_t0_r, pkB_t1_l, pkB_t1_r,
                          e00, e00_t, e00_rt );
          get2dDiff_fast( qk_t0_l, qk_t0_r, qkA_t1_l, qkA_t1_r,
                          qk_t0_l, qk_t0_r, qkB_t1_l, qkB_t1_r,
                          e11, e11_t, e11_rt );

#ifdef __combi2dSmoothMotion__
          // some Hack combine both motions 
          P4 f00_m ( f00_t[0], f00_t[1], f00_rt[0], f00_rt[1] );
          P4 f11_m ( f11_t[0], f11_t[1], f11_rt[0], f11_rt[1] );
          P4 g00_m ( g00_t[0], g00_t[1], g00_rt[0], g00_rt[1] );
          P4 g11_m ( g11_t[0], g11_t[1], g11_rt[0], g11_rt[1] );
          P4 h00_m ( h00_t[0], h00_t[1], h00_rt[0], h00_rt[1] );
          P4 h11_m ( h11_t[0], h11_t[1], h11_rt[0], h11_rt[1] );
          P4 e00_m ( e00_t[0], e00_t[1], e00_rt[0], e00_rt[1] );
          P4 e11_m ( e11_t[0], e11_t[1], e11_rt[0], e11_rt[1] );

          Scalar f00_M = sqrt( ((f00_m|f00_m) + (f11_m|f11_m) + (f00_m|f11_m) + (dN00r|dN00r)) * Scalar(3.3) + epsilon);
          Scalar f10_M = sqrt( ((h00_m|h00_m) + (h11_m|h11_m) + (h00_m|h11_m) + (dN10r|dN10r)) * Scalar(3.3) + epsilon);
          Scalar f01_M = sqrt( ((g00_m|g00_m) + (g11_m|g11_m) + (g00_m|g11_m) + (dN01r|dN01r)) * Scalar(3.3) + epsilon);
          Scalar f11_M = sqrt( ((e00_m|e00_m) + (e11_m|e11_m) + (e00_m|e11_m) ) * Scalar(3.3) + epsilon);
#endif

          Scalar f00_D = sqrt( ( (f00|f00) + (f11|f11) + (f00|f11) + (dN00|dN00) ) * Scalar(3.3) + epsilon);
          Scalar f10_D = sqrt( ( (h00|h00) + (h11|h11) + (h00|h11) + (dN10|dN10) ) * Scalar(3.3) + epsilon);
          Scalar f01_D = sqrt( ( (g00|g00) + (g11|g11) + (g00|g11) + (dN01|dN01) ) * Scalar(3.3) + epsilon);
          Scalar f11_D = sqrt(epsilon);

          Scalar ww2 = sij_weights[pq_ids];
          sumW += ww2;

#ifndef __do_depth_only__

          f00_it += f00_D * ww2;
          f10_it += f10_D * ww2;
          f01_it += f01_D * ww2;
          f11_it += f11_D * ww2;

          f00_it_M += f00_M * ww2;
          f10_it_M += f10_M * ww2;
          f01_it_M += f01_M * ww2;
          f11_it_M += f11_M * ww2;
#else
          f00_it += f00_D * ww2;
          f10_it += f10_D * ww2;
          f01_it += f01_D * ww2;
          f11_it += f11_D * ww2;
#endif
      }
      
      F00[Sij] = (min(f00_it/sumW, depthJump) + rotWeight * min(f00_it_M/sumW, rotJump))*sumW;
      F01[Sij] = (min(f01_it/sumW, depthJump) + rotWeight * min(f01_it_M/sumW, rotJump))*sumW;
      F10[Sij] = (min(f10_it/sumW, depthJump) + rotWeight * min(f10_it_M/sumW, rotJump))*sumW;
      F11[Sij] = (min(f11_it/sumW, depthJump) + rotWeight * min(f11_it_M/sumW, rotJump))*sumW;

      //      F00[Sij] = (min(sqrt(f00_it/sumW), depthJump) + rotWeight * min(sqrt(f00_it_M/sumW), rotJump))*sumW;
      //      F01[Sij] = (min(sqrt(f01_it/sumW), depthJump) + rotWeight * min(sqrt(f01_it_M/sumW), rotJump))*sumW;
      //      F10[Sij] = (min(sqrt(f10_it/sumW), depthJump) + rotWeight * min(sqrt(f10_it_M/sumW), rotJump))*sumW;
      //      F11[Sij] = (min(sqrt(f11_it/sumW), depthJump) + rotWeight * min(sqrt(f11_it_M/sumW), rotJump))*sumW;

      weightMap[ Idk[Sij] ].insert( std::pair<int,Scalar> ( Idk[Sij], F10[Sij] + F01[Sij] - F11[Sij] - F00[Sij] ) );
      score += max(max(max(F00[Sij], F01[Sij]), F10[Sij]), F11[Sij]);

      if ( score != score || std::numeric_limits<Scalar>::infinity()==score )
      {
        printf ("WHAT WHY %f \n", score);
        mexEvalString("drawnow");
        int breakhere = 0;
      }  
    }
    return score;
  }

  // try to apply the thresholdig after the accumulation
  Scalar compute_score_Fuse( std::vector<int>& currentSolution, std::vector<int>& trialIds )
  {
    compute_InitMappings( currentSolution );
    compute_InitTrial( trialIds );
    Scalar score =0;

//#pragma omp parallel for schedule (static)
    for( int Sij = 0; Sij < Idk.size(); Sij++)
    {
      Scalar f00_it(0.);
      Scalar f01_it(0.);
      Scalar f10_it(0.);
      Scalar f11_it(0.);
 
      Scalar f00_it_M(0.);
      Scalar f01_it_M(0.);
      Scalar f10_it_M(0.);
      Scalar f11_it_M(0.);

      std::vector<std::pair<int,int> >& pids = segIJ_P_ids[Sij];
      std::vector<std::pair<int,int> >& qids = segIJ_Q_ids[Sij];
      std::vector< Scalar >& sij_weights     = segIJ_W[Sij];

      const int segI = Idk[Sij];
      const int segJ = Idl[Sij];

      const std::vector< P2 >& segP2d_0l_I = segP2d_0l[segI];
      const std::vector< P2 >& segP2d_0l_J = segP2d_0l[segJ];
      const std::vector< P2 >& segP2d_0r_I = segP2d_0r[segI];
      const std::vector< P2 >& segP2d_0r_J = segP2d_0r[segJ];
      const std::vector< P2 >& segP2d_1l_I = segP2d_1l[segI];
      const std::vector< P2 >& segP2d_1l_J = segP2d_1l[segJ];
      const std::vector< P2 >& segP2d_1r_I = segP2d_1r[segI];
      const std::vector< P2 >& segP2d_1r_J = segP2d_1r[segJ];

      const std::vector< P2 >& t_segP2d_0l_I = t_segP2d_0l[segI];
      const std::vector< P2 >& t_segP2d_0l_J = t_segP2d_0l[segJ];
      const std::vector< P2 >& t_segP2d_0r_I = t_segP2d_0r[segI];
      const std::vector< P2 >& t_segP2d_0r_J = t_segP2d_0r[segJ];
      const std::vector< P2 >& t_segP2d_1l_I = t_segP2d_1l[segI];
      const std::vector< P2 >& t_segP2d_1l_J = t_segP2d_1l[segJ];
      const std::vector< P2 >& t_segP2d_1r_I = t_segP2d_1r[segI];
      const std::vector< P2 >& t_segP2d_1r_J = t_segP2d_1r[segJ];

      //P3 dN00r = roti * norm0*dn_i   - rotj * norm1*dn_j;
      //P3 dN01r = roti * norm0*dn_i   - rotk * normFix*dn_k;
      //P3 dN10r = rotk * normFix*dn_k - rotj * norm1*dn_j;
      //P3 dN11r(0.,0.,0.) ;//= rotk * normFix*dn_k - rotk * normFix*dn_k;

      P3 dN00r = gamma*(n1[ segI ] - n1[ segJ ]);
      P3 dN00  = gamma*(n0[ segI ] - n0[ segJ ]);
      /*
      P3 dN01r = gamma*(n1[ segI ] - t_n1);
      P3 dN01  = gamma*(n0[ segI ] - t_n0);
      P3 dN10r = gamma*(t_n1 - n1[ segJ ]);
      P3 dN10  = gamma*(t_n0 - n0[ segJ ]);
      P3 dN11r(0,0,0);
      P3 dN11 (0,0,0);
      */
      P3 dN01r = gamma*(n1[ segI ] - tv_n1[segJ]);
      P3 dN01  = gamma*(n0[ segI ] - tv_n0[segJ]);
      P3 dN10r = gamma*(tv_n1[segI] - n1[ segJ ]);
      P3 dN10  = gamma*(tv_n0[segJ] - n0[ segJ ]);
      P3 dN11r = gamma*(tv_n1[segI] - tv_n1[segJ]);
      P3 dN11  = gamma*(tv_n0[segI] - tv_n0[segJ]);


      Scalar sumW(0);

      for (int pq_ids = 0; pq_ids < pids.size(); pq_ids++ )
      {
        P2 f00, f00_t, f00_rt, f11, f11_t, f11_rt;
        P2 g00, g00_t, g00_rt, g11, g11_t, g11_rt;
        P2 h00, h00_t, h00_rt, h11, h11_t, h11_rt;
        P2 e00, e00_t, e00_rt, e11, e11_t, e11_rt;

          P2 pi_t0_l = segP2d_0l_I[ pids[pq_ids].first ];
          P2 pi_t1_l = segP2d_1l_I[ pids[pq_ids].first ];
          P2 pi_t0_r = segP2d_0r_I[ pids[pq_ids].first ];
          P2 pi_t1_r = segP2d_1r_I[ pids[pq_ids].first ];

          P2 qi_t0_l = segP2d_0l_I[ qids[pq_ids].first ];
          P2 qi_t1_l = segP2d_1l_I[ qids[pq_ids].first ];
          P2 qi_t0_r = segP2d_0r_I[ qids[pq_ids].first ];
          P2 qi_t1_r = segP2d_1r_I[ qids[pq_ids].first ];

          P2 pj_t0_l = segP2d_0l_J[ pids[pq_ids].second];
          P2 pj_t1_l = segP2d_1l_J[ pids[pq_ids].second];
          P2 pj_t0_r = segP2d_0r_J[ pids[pq_ids].second];
          P2 pj_t1_r = segP2d_1r_J[ pids[pq_ids].second];

          P2 qj_t0_l = segP2d_0l_J[ qids[pq_ids].second];
          P2 qj_t1_l = segP2d_1l_J[ qids[pq_ids].second];
          P2 qj_t0_r = segP2d_0r_J[ qids[pq_ids].second];
          P2 qj_t1_r = segP2d_1r_J[ qids[pq_ids].second];


          P2 pk_t0_l  = t_segP2d_0l_I[ pids[pq_ids].first ];// same as for second
          P2 pkA_t1_l = t_segP2d_1l_I[ pids[pq_ids].first ];
          P2 pk_t0_r  = t_segP2d_0r_I[ pids[pq_ids].first ];
          P2 pkA_t1_r = t_segP2d_1r_I[ pids[pq_ids].first ];
          P2 pkB_t1_l = t_segP2d_1l_J[ pids[pq_ids].second];
          P2 pkB_t1_r = t_segP2d_1r_J[ pids[pq_ids].second];

          P2 qk_t0_l  = t_segP2d_0l_I[ qids[pq_ids].first ];
          P2 qkA_t1_l = t_segP2d_1l_I[ qids[pq_ids].first ];
          P2 qk_t0_r  = t_segP2d_0r_I[ qids[pq_ids].first ];
          P2 qkA_t1_r = t_segP2d_1r_I[ qids[pq_ids].first ];
          P2 qkB_t1_l = t_segP2d_1l_J[ qids[pq_ids].second];
          P2 qkB_t1_r = t_segP2d_1r_J[ qids[pq_ids].second];


          // four times pi qj pk qk combinations
          get2dDiff_fast( pi_t0_l, pi_t0_r, pi_t1_l, pi_t1_r,
                          pj_t0_l, pj_t0_r, pj_t1_l, pj_t1_r,
                          f00, f00_t, f00_rt );
          get2dDiff_fast( qi_t0_l, qi_t0_r, qi_t1_l, qi_t1_r,
                          qj_t0_l, qj_t0_r, qj_t1_l, qj_t1_r,
                          f11, f11_t, f11_rt );

          get2dDiff_fast( pi_t0_l, pi_t0_r, pi_t1_l, pi_t1_r,
                          pk_t0_l, pk_t0_r, pkB_t1_l, pkB_t1_r,
                          g00, g00_t, g00_rt );
          get2dDiff_fast( qi_t0_l, qi_t0_r, qi_t1_l, qi_t1_r,
                          qk_t0_l, qk_t0_r, qkB_t1_l, qkB_t1_r,
                          g11, g11_t, g11_rt );

          get2dDiff_fast( pk_t0_l, pk_t0_r, pkA_t1_l, pkA_t1_r,
                          pj_t0_l, pj_t0_r, pj_t1_l, pj_t1_r,
                          h00, h00_t, h00_rt );
          get2dDiff_fast( qk_t0_l, qk_t0_r, qkA_t1_l, qkA_t1_r,
                          qj_t0_l, qj_t0_r, qj_t1_l, qj_t1_r,
                          h11, h11_t, h11_rt );

          get2dDiff_fast( pk_t0_l, pk_t0_r, pkA_t1_l, pkA_t1_r,
                          pk_t0_l, pk_t0_r, pkB_t1_l, pkB_t1_r,
                          e00, e00_t, e00_rt );
          get2dDiff_fast( qk_t0_l, qk_t0_r, qkA_t1_l, qkA_t1_r,
                          qk_t0_l, qk_t0_r, qkB_t1_l, qkB_t1_r,
                          e11, e11_t, e11_rt );

#ifdef __combi2dSmoothMotion__
          // some Hack combine both motions 
          P4 f00_m ( f00_t[0], f00_t[1], f00_rt[0], f00_rt[1] );
          P4 f11_m ( f11_t[0], f11_t[1], f11_rt[0], f11_rt[1] );
          P4 g00_m ( g00_t[0], g00_t[1], g00_rt[0], g00_rt[1] );
          P4 g11_m ( g11_t[0], g11_t[1], g11_rt[0], g11_rt[1] );
          P4 h00_m ( h00_t[0], h00_t[1], h00_rt[0], h00_rt[1] );
          P4 h11_m ( h11_t[0], h11_t[1], h11_rt[0], h11_rt[1] );
          P4 e00_m ( e00_t[0], e00_t[1], e00_rt[0], e00_rt[1] );
          P4 e11_m ( e11_t[0], e11_t[1], e11_rt[0], e11_rt[1] );

          Scalar f00_M = sqrt( ((f00_m|f00_m) + (f11_m|f11_m) + (f00_m|f11_m) + (dN00r|dN00r)) * Scalar(3.3) + epsilon);
          Scalar f10_M = sqrt( ((h00_m|h00_m) + (h11_m|h11_m) + (h00_m|h11_m) + (dN10r|dN10r)) * Scalar(3.3) + epsilon);
          Scalar f01_M = sqrt( ((g00_m|g00_m) + (g11_m|g11_m) + (g00_m|g11_m) + (dN01r|dN01r)) * Scalar(3.3) + epsilon);
          Scalar f11_M = sqrt( ((e00_m|e00_m) + (e11_m|e11_m) + (e00_m|e11_m) + (dN11r|dN11r)) * Scalar(3.3) + epsilon);
#endif

          Scalar f00_D = sqrt( ( (f00|f00) + (f11|f11) + (f00|f11) + (dN00|dN00) ) * Scalar(3.3) + epsilon);
          Scalar f10_D = sqrt( ( (h00|h00) + (h11|h11) + (h00|h11) + (dN10|dN10) ) * Scalar(3.3) + epsilon);
          Scalar f01_D = sqrt( ( (g00|g00) + (g11|g11) + (g00|g11) + (dN01|dN01) ) * Scalar(3.3) + epsilon);
          Scalar f11_D = sqrt( ( (e00|e00) + (e11|e11) + (e00|e11) + (dN11|dN11) ) * Scalar(3.3) + epsilon);
//          Scalar f11_D = sqrt(epsilon);

          Scalar ww2 = sij_weights[pq_ids];
          sumW += ww2;

#ifndef __do_depth_only__

          f00_it += f00_D * ww2;
          f10_it += f10_D * ww2;
          f01_it += f01_D * ww2;
          f11_it += f11_D * ww2;

          f00_it_M += f00_M * ww2;
          f10_it_M += f10_M * ww2;
          f01_it_M += f01_M * ww2;
          f11_it_M += f11_M * ww2;
#else
          f00_it += f00_D * ww2;
          f10_it += f10_D * ww2;
          f01_it += f01_D * ww2;
          f11_it += f11_D * ww2;
#endif
      }
      
      F00[Sij] = (min(f00_it/sumW, depthJump) + rotWeight * min(f00_it_M/sumW, rotJump))*sumW;
      F01[Sij] = (min(f01_it/sumW, depthJump) + rotWeight * min(f01_it_M/sumW, rotJump))*sumW;
      F10[Sij] = (min(f10_it/sumW, depthJump) + rotWeight * min(f10_it_M/sumW, rotJump))*sumW;
      F11[Sij] = (min(f11_it/sumW, depthJump) + rotWeight * min(f11_it_M/sumW, rotJump))*sumW;

      //      F00[Sij] = (min(sqrt(f00_it/sumW), depthJump) + rotWeight * min(sqrt(f00_it_M/sumW), rotJump))*sumW;
      //      F01[Sij] = (min(sqrt(f01_it/sumW), depthJump) + rotWeight * min(sqrt(f01_it_M/sumW), rotJump))*sumW;
      //      F10[Sij] = (min(sqrt(f10_it/sumW), depthJump) + rotWeight * min(sqrt(f10_it_M/sumW), rotJump))*sumW;
      //      F11[Sij] = (min(sqrt(f11_it/sumW), depthJump) + rotWeight * min(sqrt(f11_it_M/sumW), rotJump))*sumW;

//      weightMap[ Idk[Sij] ].insert( std::pair<int,Scalar> ( Idk[Sij], F10[Sij] + F01[Sij] - F11[Sij] - F00[Sij] ) );
      score += max(max(max(F00[Sij], F01[Sij]), F10[Sij]), F11[Sij]);

      if ( score != score || std::numeric_limits<Scalar>::infinity()==score )
      {
        printf ("WHAT WHY %f \n", score);
//        mexEvalString("drawnow");
        int breakhere = 0;
      }  
    }
    return score;
  }


  /// returns the map describing the submodularity of the edge between a segment and its neighbour
  std::vector< std::map< int, Scalar> >&  getWeightMap() { return weightMap; };


private:

  /*! 
  idea is to scale the distances computed accordingly: 2 points of different planes have a distance of 0.1
  but 1 plane is close to the camera so that scaleDepth becomes 5
  */
  inline void getJumpThresh( Scalar depth, Scalar& newThresh)//, Scalar& newScaleDepth, Scalar& newScaleMotion ) 
  {
    newThresh  = min( newThresh,  depth*depth/(depth - epix));
    // e.g. depth == 10, then 100/(390) = 1/4=0.25 -> factor 2 if jump = 0.5
//    newScale   = max( newScale,   depthJump/newThresh );
//    newScaleDepth   = max( newScaleDepth,  depthJump/newThresh );
//    newScaleMotion  = max( newScaleMotion, rotJump  /newThresh );
  };

  inline void getScale( Scalar thresh1, Scalar thresh2, Scalar& newScaleDepth, Scalar& newScaleMotion ) 
  {
    Scalar thresh   = min( thresh1, thresh2);
    newScaleDepth   = max( Scalar(1), depthJump/thresh );
    newScaleMotion  = max( Scalar(1), rotJump  /thresh );
  }

  /// for a general camera setup: p and q are the 3d points at timstep t0 (of the same 2d point)
  inline void get2dDiff( P3 p_t0, P3 q_t0, P3 p_t1, P3 q_t1, P2& dispDifft0, P2& dispDifft1, P2& flowDiff )
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

  // do that for all 8 points
  inline void get2dprojections( const P3& p_t0, P2& pixP_t0_l, P2& pixP_t0_r )
  {
    // pixel directly:
    P3 pixP_t0_lT = Pl * p_t0;
    P3 pixP_t0_rT = pixP_t0_lT +pr;
    pixP_t0_lT /= pixP_t0_lT[2];
    pixP_t0_rT /= pixP_t0_rT[2];
    pixP_t0_l = P2(pixP_t0_lT[0], pixP_t0_lT[1]);
    pixP_t0_r = P2(pixP_t0_rT[0], pixP_t0_rT[1]);
  }

  inline void get2dDiff_fast( const P2& pixP_t0_l, const P2& pixP_t0_r, const P2& pixP_t1_l, const P2& pixP_t1_r, 
                              const P2& pixQ_t0_l, const P2& pixQ_t0_r, const P2& pixQ_t1_l, const P2& pixQ_t1_r,
                              P2& dispDifft0, P2& dispDifft1, P2& flowDiff )
  {
    P2 difft0 = (pixP_t0_r-pixP_t0_l) - (pixQ_t0_r-pixQ_t0_l);
    P2 difft1 = (pixP_t1_r-pixP_t1_l) - (pixQ_t1_r-pixQ_t1_l);
    flowDiff  = (pixP_t1_l-pixP_t0_l) - (pixQ_t1_l-pixQ_t0_l);
    ///////////////////
    dispDifft0 = P2(difft0[0],difft0[1]);
    dispDifft1 = P2(difft1[0]-difft0[0], difft1[1]-difft0[1]); // correct but non sense different difference to cam so also different disp difference or not ?
  }



  /// for a non-general camera setup: p and q are the 3d points at timstep t0 (of the same 2d point), here special case Pl=Pr, p_t0 and q_t0 must be in pixel coords, pr similarly
  inline void get2dDiff_fake( const P3& p_t0, const P3& q_t0, const P3& p_t1, const P3& q_t1, P2& dispDifft0, P2& dispDifft1, P2& flowDiff )
  {
    // pixel directly:
    P3 pixP_t0_l = Pl * p_t0;
    P3 pixQ_t0_l = Pl * q_t0;
    // assume last row in pl = 0,0,1, so: ok
    P3 pixP_t0_r = pixP_t0_l +pr;
    P3 pixQ_t0_r = pixQ_t0_l +pr;
    pixP_t0_l /= pixP_t0_l[2];
    pixQ_t0_l /= pixQ_t0_l[2];
    pixP_t0_r /= pixP_t0_r[2];
    pixQ_t0_r /= pixQ_t0_r[2];

    P3 pixP_t1_l = Pl * p_t1;
    P3 pixP_t1_r = pixP_t1_l + pr; 
    pixP_t1_r /= pixP_t1_r[2];
    pixP_t1_l /= pixP_t1_l[2];

    P3 pixQ_t1_l = Pl * q_t1;
    P3 pixQ_t1_r = pixQ_t1_l + pr;
    pixQ_t1_l /= pixQ_t1_l[2];
    pixQ_t1_r /= pixQ_t1_r[2];


    P3 difft0 = (pixP_t0_r-pixP_t0_l) - (pixQ_t0_r-pixQ_t0_l);
    P3 difft1 = (pixP_t1_r-pixP_t1_l) - (pixQ_t1_r-pixQ_t1_l);
    P3 diffF  = (pixP_t1_l-pixP_t0_l) - (pixQ_t1_l-pixQ_t0_l);
    ///////////////////

    flowDiff   = P2(diffF[0],  diffF[1]);
    dispDifft0 = P2(difft0[0], difft0[1]);
    dispDifft1 = P2(difft1[0]-difft0[0], difft1[1]-difft0[1]); // correct but non sense different difference to cam so also different disp difference or not ?
  }


  /// for a general camera setup: p and q are the 3d points at timstep t0 (of the same 2d point)
  inline void get2d_3dDiff( const P3& p_t0, const P3& q_t0, const P3& p_t1, const P3& q_t1, P2& dispDifft0, P3& motionDiff )
  {
    P3 pixP_t0_l = Pl * p_t0;pixP_t0_l /= pixP_t0_l[2];
    P3 pixQ_t0_l = Pl * q_t0;pixQ_t0_l /= pixQ_t0_l[2];

//    P3 pixP_t1_l = Pl * p_t1;pixP_t1_l /= pixP_t1_l[2];
//    P3 pixQ_t1_l = Pl * q_t1;pixQ_t1_l /= pixQ_t1_l[2];

    P3 pixP_t0_r = Pr * p_t0 +pr;pixP_t0_r /= pixP_t0_r[2];
    P3 pixQ_t0_r = Pr * q_t0 +pr;pixQ_t0_r /= pixQ_t0_r[2];

//    P3 pixP_t1_r = Pr * p_t1 +pr; pixP_t1_r /= pixP_t1_r[2];
//    P3 pixQ_t1_r = Pr * q_t1 +pr; pixQ_t1_r /= pixQ_t1_r[2];

    P3 difft0 = (pixP_t0_r-pixP_t0_l) - (pixQ_t0_r-pixQ_t0_l);
//    P3 difft1 = (pixP_t1_r-pixP_t1_l) - (pixQ_t1_r-pixQ_t1_l);
//    P3 diffF  = (pixP_t1_l-pixP_t0_l) - (pixQ_t1_l-pixQ_t0_l);
    ///////////////////

    motionDiff = (p_t1-p_t0) - (q_t1-q_t0);// vector a minus vector b

    // 5 too much ? before *2: decent
//    motionDiff[0] *= 1;
//    motionDiff[1] *= 1;

    dispDifft0 = P2(difft0[0], difft0[1]);
  }


  /// simplified special case assumes rectified cameras
  inline Scalar getDisparity( Scalar d_i, Scalar d_j )
  {
    // no absolute value!
    return ( epix* (1/d_i-1/d_j) );
  }

  /// 4 views disparity, optical flow and disparity distance
  inline void get4D_Disparity( P3 p_i, P3 p_j ) {}

  inline void scaleMotionXY( P3& f00, P3& f01, P3& f10, P3& f11 )
  {
    f00[0] *= 1;
    f00[1] *= 1;
    f10[0] *= 1;
    f10[1] *= 1;
    f01[0] *= 1;
    f01[1] *= 1;
    f11[0] *= 1;
    f11[1] *= 1;
  }

  inline void getScales( Scalar dp_i, Scalar dq_i, Scalar dp_j, Scalar dq_j, Scalar dp_k, Scalar dq_k,
    Scalar& newScaleDepth_ij, Scalar& newScaleMotion_ij,  Scalar& newScaleDepth_ik, Scalar& newScaleMotion_ik, 
    Scalar& newScaleDepth_jk, Scalar& newScaleMotion_jk, Scalar& newScaleDepth_k, Scalar& newScaleMotion_k )
  {
    Scalar threshi = depthJump;
    Scalar threshj = depthJump;
    Scalar threshk = depthJump;

    getJumpThresh( dp_i, threshi );
    getJumpThresh( dq_i, threshi );

    getJumpThresh( dp_j, threshj );
    getJumpThresh( dq_j, threshj );

    getJumpThresh( dp_k, threshk );
    getJumpThresh( dq_k, threshk );

    getScale( threshi, threshj, newScaleDepth_ij, newScaleMotion_ij);
    getScale( threshi, threshk, newScaleDepth_ik, newScaleMotion_ik);
    getScale( threshj, threshk, newScaleDepth_jk, newScaleMotion_jk);
    getScale( threshk, threshk, newScaleDepth_k, newScaleMotion_k);
  }

  //////////////////////////////////////////////////

  std::vector< int > segIJMap;

  /// weights per i,j seg pair pixel segIJ_W[ Seg(i,j)-Id ] [1..N]
  std::vector< std::vector< Scalar > > segIJ_W;

  /// store p2d per segment - to be tranformed per motion and store the current situation
  std::vector< std::vector< P3 > > segP2d;
//  std::vector< std::vector< P3 > > segQ2d;

  std::vector< std::vector<std::pair<int,int> > > segIJ_P_ids;
  std::vector< std::vector<std::pair<int,int> > > segIJ_Q_ids;

  // projected and transofrmed points by the current solution
  std::vector< std::vector< P2 > > segP2d_0l;
  std::vector< std::vector< P2 > > segP2d_1l;
  std::vector< std::vector< P2 > > segP2d_0r;
  std::vector< std::vector< P2 > > segP2d_1r;

  // projected and transofrmed points by the trial solution
  std::vector< std::vector< P2 > > t_segP2d_0l;
  std::vector< std::vector< P2 > > t_segP2d_1l;
  std::vector< std::vector< P2 > > t_segP2d_0r;
  std::vector< std::vector< P2 > > t_segP2d_1r;

  std::vector< P3 > n0;
  std::vector< P3 > n1;

  P3 t_n0;
  P3 t_n1;

  // trial normals, etc.
  std::vector< P3 > tv_n0;
  std::vector< P3 > tv_n1;

  M3 Pl;//R=I t=0

  M3 Pr;//K*R
  P3 pr;//K*t
  std::vector<std::map<int, Scalar> > edgeWeightMap;

  /// the edges 0,1,2 are 2d coordinates of point p 
  std::vector<P3> edgesP;
  /// the edges 0,1,2 are 2d coordinates of point q
  std::vector<P3> edgesQ;

  /// the centers of the pathces in 2d
//  std::vector<P3> centers;
  std::vector<P3> centersA;
  std::vector<P3> centersB;

  /// the weights of the edges
  std::vector<Scalar> weights;

  /// ids of the patches id k has a corresponding neighbour l in Idl
  std::vector<int> Idk;

  /// ids of the patches id l has a corresponding neighbour k in Idl
  std::vector<int> Idl;

  const std::vector<P3>* normals;

  /// translations
  const std::vector<P3>* tra;
  /// rotations
  const std::vector<M3>* rot;


  std::vector<Scalar> F00;
  std::vector<Scalar> F01;
  std::vector<Scalar> F10;
  std::vector<Scalar> F11;
  
  /// weights on the edges per pixel (4/8-neighbourhood)
  Scalar* halfEdgeX;
  Scalar* halfEdgeY;

  /// weight on the cross edges (8-neighbourhood)
  Scalar* halfEdgeXY;
  Scalar* halfEdgeiXY;

  Scalar gamma;
  double epsilon;
  Scalar rotJump;
  Scalar depthJump;
  Scalar rotWeight;

  int nSegments;

  /// used to adjust the threshold in case of 2d jumps of more than e.g. 1 pixel (given as parameter)
  Scalar epix;

  /// keeps track of the weights at the edges, test: use for non-submodular edge selection for occlusion mapping
  std::vector< std::map< int, Scalar> > weightMap;
  
  /// simplified camera configuration: Pl == Pr
  bool simpleProjection;

  M3 iK;
};
#endif
#undef __do_strict_min__
#undef __do_depth_only__
