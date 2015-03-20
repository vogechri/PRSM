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

#ifndef _BILINEAR_IP_H
#define _BILINEAR_IP_H

#include <math.h>

#define _mm_extract_pd(R,I) (*((double*)(&R)+I))

template <typename TScalar>
class BiLinearIP
{
public:

  typedef TScalar Scalar;

  BiLinearIP(int _N, int _M, Scalar*& _Img ) 
    : N(_N), M(_M), I(_Img), maxConstN(Scalar(_N)-1.000001), maxConstM(Scalar(_M)-1.000001)
  {
    mm = _mm_set1_epi32(M);
    ones = _mm_set1_pd( 1.0);
    maxConstM_128d = _mm_set1_pd(maxConstM);
    maxConstN_128d = _mm_set1_pd(maxConstN);
    minConst_128d  = _mm_set1_pd(0.000001);
  };

  ~BiLinearIP() {};

  void interpolate_noOmp_short( int elements, short *&inputF, Scalar*& X, Scalar *&Y, short *&outputF )
  {
    for(int id=0;id<elements;id++)
    {
        // c++ style indices not matlab
        Scalar idx = X[id] - Scalar(1.0);
        Scalar idy = Y[id] - Scalar(1.0);

        idx = std::max( Scalar(0.000001), std::min(idx, maxConstN) );
        idy = std::max( Scalar(0.000001), std::min(idy, maxConstM) );

        outputF[id] = GetPixel_16(inputF, idx, idy);
    }
  }

  void interpolate_short( short *&inputF, int elements, Scalar*& X, Scalar *&Y, short *&outputF )
  {
    for(int id=0;id<elements;id++)
    {
        // c++ style indices not matlab
        Scalar idx = X[id] -Scalar(1.0);
        Scalar idy = Y[id] -Scalar(1.0);

        idx = std::max( Scalar(0.00000001), std::min(idx, Scalar(N-1)) );
        idy = std::max( Scalar(0.00000001), std::min(idy, Scalar(M-1)) );

        outputF[id] = GetPixel_16(inputF, idx, idy);
    }
  }
 

  void interpolate_box_short( short *&inputF, Scalar*& X, Scalar *&Y, short *&outputF, int xmin, int xmax, int ymin, int ymax )
  {
    for(int i=xmin;i<xmax;i++)
      for(int j=ymin;j<ymax;j++)
      {
        int id = i*M + j;

        // c++ style indices not matlab
        Scalar idx = X[id] -Scalar(1.0);
        Scalar idy = Y[id] -Scalar(1.0);

        idx = std::max( Scalar(0.00000001), std::min(idx, Scalar(N-1)) );
        idy = std::max( Scalar(0.00000001), std::min(idy, Scalar(M-1)) );

        outputF[id] = GetPixel_16(inputF, idx, idy);
      }
  }

  inline short GetPixel_16 (short*& img, Scalar &x, Scalar &y)
  {
    unsigned int Fx = (unsigned int)(x * 65536.0); // convert to Fixed8
    unsigned int Fy = (unsigned int)(y * 65536.0); // convert to Fixed8
    unsigned int px = (Fx & 0xFFFF0000) >> 16; // floor function
    unsigned int py = (Fy & 0xFFFF0000) >> 16; // floor function

    // the point: out of bounds pixel have weight 0 - 
    // see caller above; x is N-1 at most - y M-1
    // therefore the images are expanded !
    //if (px * M + py + M > N*M-1 || px * M + py + M < 0) can happen
    const short* p0 = img + px * M + py; // pointer to first pixel

    const unsigned int fx = Fx & 0x0000FFFF; // frac function
    const unsigned int fy = Fy & 0x0000FFFF; // frac function

    // load the four neighboring pixels
/*    const int p1 = p0[0];
    const int p3 = p0[1];
    const int p2 = p0[0 + M];
    const int p4 = p0[1 + M];
*/
    int p1 = *((int*) p0);
    int p2 = *((int*) (&p0[M]));
    int p3 = p1>>16;p1 = p1&0x0000FFFF;
    int p4 = p2>>16;p2 = p2&0x0000FFFF;

    int w2 = (fx * (0x00010000 - fy))  >> 16;
    int w3 = ((0x00010000 - fx) * fy ) >> 16;
    int w4 = (fx * fy )  >> 16;
//    int w1 = ((0x00010000 - fx) * (0x00010000 - fy) ) >> 16; // EVIL LINE: 1*1 = overflow
    int w1 = 0x00010000 - (w4+w2+w3); // no overflow

    // convert to short here: there are overflows for some reason ?!
    short outr = (short) (0x00007FFF & ((p1 * w1 + p2 * w2 + p3 * w3 + p4 * w4)>>9));
    return outr;
  }

  void approximateImageGradients() {};

private:

  int N;
  int M;

  const Scalar maxConstN;
  const Scalar maxConstM;

  const Scalar* I;
//  __declspec(align(16)) const Scalar* I;
  __m128d ones;
  __m128i mm;

  __m128d maxConstM_128d;
  __m128d maxConstN_128d;
  __m128d minConst_128d;
};

#endif
