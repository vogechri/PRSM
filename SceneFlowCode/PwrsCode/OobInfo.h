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

#ifndef _OOB_INFO_H
#define _OOB_INFO_H

//////////////////////////////////////////////////////////////////////////////////////
/// informs whether a segment is occluded or not -> drop the truncation here
template<typename Scalar> class oobInfo
{
public:

  typedef Math::Mat3x3T<Scalar>     M3;
  typedef Math::VectorT<Scalar, 3>  P3;

  oobInfo( int nSegments_, int nAssignments_ ) 
    : nAssignments(nAssignments_), nSegments(nSegments_) 
  { oob.resize(nAssignments_ * nSegments_, false); };

  ~oobInfo( ){};

  inline bool occluded(int segment, int assignment)
  {
    return oob[ assignment*nSegments + segment ] ;
  }

  void getOcclusionVec( std::vector<int> solution, int proposal, std::vector<bool>& oobs)
  {
    std::copy( oob.begin()+proposal*nSegments, oob.begin()+(proposal+1)*nSegments, oobs.begin() );
    for(int i=0;i<solution.size();i++)
      if( occluded(i, solution[i]) ) 
        oobs[i] = true;

    return;// oobs;
  }

  void getOcclusionVec( int proposal, std::vector<bool>& oobs)
  {
    std::copy( oob.begin()+proposal*nSegments, oob.begin()+(proposal+1)*nSegments, oobs.begin() );
    return;
  }

  void getOcclusionVec( std::vector<int> solution, std::vector<bool>& oobs)
  {
    for(int i=0;i<solution.size();i++)
      if( occluded(i, solution[i]) ) 
        oobs[i] = true;
    return;
  }
  

  inline void setOccluded(std::vector<bool> occlusions, int assignment)
  {
    std::copy( occlusions.begin(), occlusions.end(), oob.begin() + assignment*nSegments );
  }

  inline void setOccluded(int segment, int assignment)
  {
    return oob[ assignment*nSegments + segment ] = true;
  }
   
  /// not needed
  inline void unsetOccluded(int segment, int assignment)
  {
    return oob[ assignment*nSegments + segment ] = false;
  }

  //////////////////////////////////////////////////
  std::vector<bool> oob;
  int nSegments;
  int nAssignments;
};


#endif