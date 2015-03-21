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
License along with this software.
*/

#include <iostream>  // iostream must be before mex -> else redefinition
#include "mex.h"

#include "FullStepMatlab.h"

void mexFunction ( int nlhs, mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{
  growGrid ( nlhs, plhs, nrhs, prhs );
}
