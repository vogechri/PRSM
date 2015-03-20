// StereoTest.cpp : Defines the entry point for the console application.
//

#include <opencv2/opencv.hpp>
#include <stdio.h>

#define _matlab_output

#ifdef _matlab_output
#include "mex.h"
#endif

#define opencv243

using namespace cv;
void print_help()
{
    printf("Usage: stereo_match <left_image> <right_image> [--algorithm=bm|sgbm|hh] [--blocksize=<block_size>]\n"
           "[--max-disparity=<max_disparity>] [-i <intrinsic_filename>] [-e <extrinsic_filename>]\n"
           "[--no-display] [-o <disparity_image>] [-p <point_cloud_file>]\n");
}


// call (N,M, img1, img2, M1,M2,D1,D2,R,T, method, maxDisp, blockSize)
void get_Stereo ( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

const int invertMap = 0;

int N, M, minDisp, maxDisp, blockSize;
int doRectification = 1;

int elements, nsubs;
double *output_Points, *output_Q;

char buffer[256];

N         = (int) *mxGetPr(prhs[0]);
M         = (int) *mxGetPr(prhs[1]);
minDisp   = (int) *mxGetPr(prhs[10]);
maxDisp   = (int) *mxGetPr(prhs[11]);
blockSize = (int) *mxGetPr(prhs[12]);

int p1(50),p2(800);

if (nrhs > 13) 
  doRectification = (int) *mxGetPr(prhs[13]);

if (nrhs > 14) 
  p1 = (int) *mxGetPr(prhs[14]);
if (nrhs > 15) 
  p2 = (int) *mxGetPr(prhs[15]);

#ifdef _matlab_output
if (! (nrhs == 13 || nrhs == 14 || nrhs == 15 || nrhs == 16))
  mexErrMsgTxt("Must have 13, 14, 15 or 16 input argument");

if (nlhs < 2)
  mexErrMsgTxt("Must have at least two output arguments");
#endif

/* Get the number of elements in the input argument */
#ifdef _matlab_output

elements=mxGetNumberOfElements(prhs[3]);

sprintf (buffer, "The amount of elements in the matrix is not equal to the amount given by the specified dimensions!\n all: %d product %d   N: %d, M: %d\n", elements, N*M, N, M);
if (elements != N*M)
  mexErrMsgTxt(buffer);
#endif
/* Get the number of dimensions in array */
//nsubs=mxGetNumberOfDimensions(prhs[9]);
//msubs=mxGetNumberOfDimensions(prhs[6]);

unsigned char* ImgA, *ImgB;
double *MA, *MB, *DA, *DB, *Ri, *Ti;

ImgA = (unsigned char*) mxGetPr(prhs[2]);
ImgB = (unsigned char*) mxGetPr(prhs[3]);
MA   = mxGetPr(prhs[4]);
MB   = mxGetPr(prhs[5]);
DA   = mxGetPr(prhs[6]);
DB   = mxGetPr(prhs[7]);
Ri   = mxGetPr(prhs[8]);
Ti   = mxGetPr(prhs[9]);

Mat img1(N,M, CV_8UC1), img2(N,M, CV_8UC1);

for(int m=0;m<M;m++)//M
{
  for(int n=0;n<N;n++)//N row index
  {
    int idv = m*N+n;
    img1.at<unsigned char>(n,m) = ImgA[idv];
 		img2.at<unsigned char>(n,m) = ImgB[idv];
  }
}

// first

Size img_size = img1.size();

//Mat M1, D1, M2, D2, R, T, R1, R2, P1, P2, Q;

Mat M1t(3,3, CV_64F, MA), D1(4,1, CV_64F, DA), M2t(3,3, CV_64F, MB), D2(1,4, CV_64F, DB);
Mat Rt(3,3, CV_64F, Ri), T(3,1, CV_64F, Ti);
Mat R1, R2, P1, P2, Q = Mat::eye( 4, 4, CV_64F);
Mat M1, M2, R;
Rect roi1, roi2;

M1 = M1t.t();M2 = M2t.t();R = Rt.t();

if (doRectification)
{
// rectification ok but the images are cut A LOT TO EARLY
#ifdef opencv243
  stereoRectify( M1, D1, M2, D2, img_size, R, T, R1, R2, P1, P2, Q, 1, -1, img_size, &roi1, &roi2);//CALIB_ZERO_DISPARITY);//&roi1, &roi2 );
#else
  stereoRectify( M1, D1, M2, D2, img_size, R, T, R1, R2, P1, P2, Q, 1, img_size, &roi1, &roi2, 0);//CALIB_ZERO_DISPARITY);//&roi1, &roi2 );
#endif

Mat map11, map12, map21, map22;
initUndistortRectifyMap(M1, D1, R1, P1, img_size, CV_16SC2, map11, map12);
initUndistortRectifyMap(M2, D2, R2, P2, img_size, CV_16SC2, map21, map22);

Mat img1r, img2r;
remap(img1, img1r, map11, map12, INTER_LINEAR);
remap(img2, img2r, map21, map22, INTER_LINEAR);

img1 = img1r;
img2 = img2r;
}

Mat disp, disp8;

int method = 0; // BM stereo sucks
if (method == 1)
{

    StereoBM bm;
//    SADWindowSize = 9;
    bm.state->roi1 = roi1;
    bm.state->roi2 = roi2;
    bm.state->preFilterCap = 31;
    bm.state->SADWindowSize = blockSize > 0 ? blockSize : 9;
    bm.state->minDisparity = 0;
    bm.state->numberOfDisparities = maxDisp;
    bm.state->textureThreshold = 10;
    bm.state->uniquenessRatio = 15;
    bm.state->speckleWindowSize = 100;
    bm.state->speckleRange = 32;
    bm.state->disp12MaxDiff = 1;

    // NOTE THAT THIS IS WEIRD BUT I NEED THIS FLIP!
    bm(img2, img1, disp);
//    bm(img1, img2, disp);

}
else
{
    StereoSGBM sgbm;

    // changed some parameters according to website KITTI

    sgbm.preFilterCap = 15;//63;
    sgbm.SADWindowSize = blockSize > 0 ? blockSize : 5;
    
    int cn = img1.channels();
    if (!p1)
      sgbm.P1 = 50;//4*cn*sgbm.SADWindowSize*sgbm.SADWindowSize; // smoothness == change of -1 // before used 8 not 4
    else
      sgbm.P1 = p1;
    if (!p2)
      sgbm.P2 = 800;//8*p1;//32*cn*sgbm.SADWindowSize*sgbm.SADWindowSize;// larger changes
    else
      sgbm.P2 = p2;

    sgbm.minDisparity = minDisp;
    sgbm.numberOfDisparities = maxDisp;
    sgbm.uniquenessRatio = 0;//10;
    sgbm.speckleWindowSize = 100;
//  less but more reliable estimate - or the opposite? : more pixel defined
//    sgbm.uniquenessRatio = 60;
//    sgbm.speckleWindowSize = 200;
    sgbm.speckleRange = 32;
    sgbm.disp12MaxDiff = 1; // difference for left right check
    sgbm.fullDP  = 1;//= alg == STEREO_HH;

    // NOTE THAT THIS IS WEIRD BUT I NEED THIS !
    if (invertMap)
      sgbm(img2, img1, disp);
    else
      sgbm(img1, img2, disp);

}

Size disp_size = img1.size();//disp.size();

int width = disp_size.width;
int height = disp_size.height;

// Create the output array 

plhs[0] = mxCreateDoubleMatrix(height, width, mxREAL);
output_Points  = mxGetPr(plhs[0]);
for (int i=0;i<height;i++)
  for(int j=0;j<width;j++)
    output_Points[j*height+i] = float( disp.at<short>(i,j) ) /16.0;
//    output_Points[j*height+i] = img1.at<unsigned char>(i,j);
    
#ifdef _matlab_output
plhs[1] = mxCreateDoubleMatrix( 4, 4, mxREAL);
output_Q  = mxGetPr(plhs[1]);
#else
output_Q = (double*) malloc(sizeof(double) * 4*4);
#endif
for (int i=0;i<4;i++)
  for(int j=0;j<4;j++)
{
   output_Q[j*4+i] = Q.at<double>(i,j);
}

#ifdef _matlab_output
plhs[2] = mxCreateDoubleMatrix( 3, 3, mxREAL);
double *output_H1  = mxGetPr(plhs[2]);
plhs[3] = mxCreateDoubleMatrix( 3, 3, mxREAL);
double *output_H2  = mxGetPr(plhs[3]);
#else
double *output_H1 = (double*) malloc(sizeof(double) * 3*3);
double *output_H2 = (double*) malloc(sizeof(double) * 3*3);
#endif

Mat H1 = Mat::eye( 3, 3, CV_64F);
Mat H2 = Mat::eye( 3, 3, CV_64F);

if (doRectification)
{
  H1 = M1 * (P1.colRange(0,3) * R1).inv();
  H2 = M2 * (P2.colRange(0,3) * R2).inv();
}
#ifdef _matlab_output
for (int i=0;i<3;i++)
  for(int j=0;j<3;j++)
  {
    output_H2[j*3+i] = H2.at<double>(i,j);
    output_H1[j*3+i] = H1.at<double>(i,j);
  }
#endif

#ifndef _matlab_output
free (output_H1);
free (output_H2);
free (output_Q);
//free (output_Points);

printf("The end\n");

printf("Not yet\n");

#endif

}
