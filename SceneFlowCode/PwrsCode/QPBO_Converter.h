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

/////////////////////////////////////////////////////////////////////////////
// Converts a general pseudo boolean polynomial to one with quadratic form //
/////////////////////////////////////////////////////////////////////////////

#ifndef __QPBO_Converter_
#define __QPBO_Converter_

#include <math.h>

using namespace std;
using namespace Math;

//#define round(x)          floor(x+0.5)
#define round_matlab(x)   floor(x-0.5)

//////////////////////////////////////////////////////////////////////////////////////
struct occStorage
{
public:
  occStorage( int _seg_i, int _seg_j, double _f00, double _f01, double _f10, double _f11 )
    : seg_i(_seg_i), seg_j(_seg_j), f00(_f00), f01( _f01), f10(_f10), f11(_f11)
  {};

  ~occStorage() {};
  ///////////////////////////////////////
  int seg_i;
  int seg_j;

  double f00;
  double f10;
  double f01;
  double f11;
};

/////////////////////////////////////////////////////////////////////////

/// stores new edges
template<typename Scalar>
struct edgePenalty
{
  edgePenalty(int _seg_i, int _seg_j, Scalar _penalty)
    : seg_i(_seg_i),  seg_j(_seg_j),  f11(_penalty)
  {};

  ~edgePenalty(){};

  int seg_i;
  int seg_j;
  Scalar f11;
};
#endif
