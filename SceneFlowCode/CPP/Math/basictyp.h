/*----------------------------------------------------------------------*/
/* modul      : basictyp.h						*/
/*                                                                      */
/* description: Some basic data types which should be common to several	*/
/*		moduls.							*/
/*                                                                      */
/* author     : P. J. Neugebauer                                        */
/* date       : 05.12.96                                                */
/*----------------------------------------------------------------------*/
/* changed by : Matthias Block                                          */
/* date       : 22.04.1999                                              */
/* description: added the data types used for the preview               */
/*----------------------------------------------------------------------*/
/* changed by : Matthias Block                                          */
/* date       : 01.02.2000                                              */
/* description: removed flag in Point3D_T -> was obsolete               */
/*----------------------------------------------------------------------*/

#ifndef __BASICTYP
#define __BASICTYP

#include "math.h"

/*----------------------------------------------------------------------*/
/* datatypes                                                            */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* name       : Point3D_t						*/
/* description: Represents floating point coordinates of a point in     */
/*              three dimensions.                                       */
/*----------------------------------------------------------------------*/

typedef struct Point3DType
{
  public:
    Point3DType() {}
    Point3DType(double x, double y, double z)
    {
      Point3DType::x = x;
      Point3DType::y = y;
      Point3DType::z = z;
    }
    
    inline bool normalize()
    {
      double VecNorm = sqrt(x*x + y*y + z*z);
      if(VecNorm == 0.0)
        return false;
        
      VecNorm = 1.0/VecNorm;     
      x *= VecNorm;
      y *= VecNorm;
      z *= VecNorm;
      
      return true;
    }
    
    inline double length()
    {
      return sqrt(x*x + y*y + z*z);
    }
    
    inline double squaredLength()
    {
      return (x*x + y*y + z*z);
    }
  
  public:
    double x;
    double y;
    double z;
} Point3D_t;

/* 
// schoener ware sicherlich:
// um z.B. beim rendern glVertex3dv zu benutzen, etc
// jedoch werden zum betrachten die y werte negiert !!
typedef union Point3DType 
{
  struct
  {
    double x;
    double y;
    double z;
  };

  double vec[3];
} Point3D_t;
*/

/*----------------------------------------------------------------------*/
/* name       : IPoint3D_t						*/
/* description: Represents floating point coordinates of a point in     */
/*              three dimensions.                                       */
/*----------------------------------------------------------------------*/
typedef struct IPoint3DType {
  long x;
  long y;
  long z;
} IPoint3D_t;


/*----------------------------------------------------------------------*/
/* name       : FPoint3D_t						*/
/* description: a type used to store point coordinates in 3D-space      */
/*		as float-values						*/
/*----------------------------------------------------------------------*/
typedef struct FPoint3DType {
  float x;
  float y;
  float z;
} FPoint3D_t;


/*----------------------------------------------------------------------*/
/* name       : Point2D_t						*/
/* description: Represents floating point coordinates of a point in two */
/*              dimensions.                                             */
/*----------------------------------------------------------------------*/
typedef struct Point2DType {
  double u;
  double v;
} Point2D_t;


/*----------------------------------------------------------------------*/
/* name       : IPoint2D_t						*/
/* description: Represents floating point coordinates of a point in two	*/
/*              dimensions.						*/
/*----------------------------------------------------------------------*/
typedef struct IPoint2DType {
  long u;
  long v;
} IPoint2D_t;

/*----------------------------------------------------------------------*/
/* name       : FPoint2D_t						*/
/* description: Represents floating point coordinates of a point in two */
/*              dimensions.                                             */
/*----------------------------------------------------------------------*/
typedef struct FPoint2DType {
  float u;
  float v;
} FPoint2D_t;

/*----------------------------------------------------------------------*/
/* name       : Point4D_t						*/
/* description: Represents floating point coordinates of a point in four*/
/*              dimensions (describes a point in homogenous coords or   */
/*              a plane in n x + w == 0                                 */
/*              where n == {x,y,z} and x \in R^3                        */
/*----------------------------------------------------------------------*/
typedef struct Point4DType {
  double x;
  double y;
  double z;
  double w;
} Point4D_t;


#endif /* __BASICTYP */





