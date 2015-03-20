#ifndef _EASY_SPARSE_MATRIX_T
#define _EASY_SPARSE_MATRIX_T

//== INCLUDES =======================================================

#include "VectorT.h"
#include "GenProg.h"

//== NAMESPACE ================================================================

namespace Math {


  //== CLASS DEFINITION =========================================================


  /** /class EasySparseMatrix

  Stores a sparse symmetric! matrix efficiently.
  **/

  template <class Scalar>
  class EasySparseMatrix
  {
  public:
    typedef Scalar                   value_type;
    typedef EasySparseMatrix<Scalar> type;
    typedef EasySparseMatrix<Scalar> Self;


    EasySparseMatrix(): _init(false)  { };

    ~EasySparseMatrix() { };

    bool initialized() { return _init; }

    /// @param rSize: in: the id of the row, out the size of the row
    void getRowPtr(int &rSize, Scalar*& valPtr, int*& idPtr)
    {
      if (rSize > static_cast<int>(colptr.size()) -1 )
      {
        valPtr = NULL;
        rSize = 0;
        idPtr = 0;
      }

      valPtr = static_cast<Scalar*> (&(values[ colptr[rSize] ]));
      idPtr  = static_cast<int*> (&(rowind[ colptr[rSize] ]));
      rSize  = colptr[rSize+1] - colptr[rSize]; 
    }

    /// with a map which resolves several adds to an already set value with an addition or a set else, and does not requiere pre-sorting
    void begin_map_row()
    {
      if (colptr.empty() || colptr.back() != (int)values.size())
      {
        colptr.push_back(values.size());
//        n_rows = colptr.size()-1;
      }

      map_row.clear();
    };

    /// if unsure wether existing use this for adding the value to the possibly existing value, the entries do not need to be sorted!
    void add_map_row_value(int _i, Scalar _val)
    {
      if (map_row.find(_i) != map_row.end() )
        map_row[_i] += _val;
      else 
        map_row[_i]  = _val;
    };

    /// if sure that the entry is not existing, use this for setting the value, the entries do not need to be sorted!
    void set_map_row_value(int _i, Scalar _val)
    {
      map_row[_i] = _val;
    };

    /// finalize the row started and defined earlier, @see begin_map_row
    void end_map_row()
    {
      typename std::map<int , Scalar>::const_iterator row_it;

      begin_row();

      for (row_it = map_row.begin();row_it != map_row.end();row_it++)
      {
        // existiert index überhaupt , also steht der eintrag in der 
        // symmetrisch gemachten matrix ?
        if (row_it->first >= 0)
        {
          set_row_value(row_it->first , row_it->second);
          if (isnan(row_it->second)) std::cerr << "isNAN: " << row_it->second << "\n";
          if (isinf(row_it->second)) std::cerr << "isINF: " << row_it->second << "\n";
        }
      }

      end_row();

      map_row.clear();
    };

    //void newRow()
    //{
    //  colptr.push_back( values.size() );
    //}

    /*
    void addValue(int id, Scalar val)
    {
    if (colptr.empty() || colptr.back() != (int)values.size())
    {
    colptr.push_back(values.size());
    n_rows = colptr.size()-1;
    }

    map_row.clear();
    }
    */

    /// preallocate memory ( saves a lot of time when setting up the data)
    void reserve(int nRows, int nElements)
    {
      colptr.clear();
      colptr.reserve(nRows+1);
      values.clear();
      values.reserve(nElements);
      rowind.clear();
      rowind.reserve(nElements);
    }

    /// reset the memory use for efficiency @see reserve
    void finalize()
    {
      //      colptr.push_back(values.size());
      colptr.resize(colptr.size());
      rowind.resize(rowind.size());
      values.resize(values.size());
      _init = true;
    }

  private:

    /// internally used for ordering of the elements
    void begin_row()
    {
      if (colptr.empty() || colptr.back() != (int)values.size())
      {
        colptr.push_back(values.size());
//        n_rows = colptr.size()-1;
      }
    }

    /// internally used for ordering of the elements
    void set_row_value(int id, Scalar val)
    {
      // needed for rhs computation!
//      if (id > static_cast<int>(colptr.size()) -1 ) return;

      rowind.push_back(id);
      values.push_back(val);
    }

    /// internally used for ordering of the elements
    void end_row()
    {
      if (colptr.empty() || colptr.back() != (int)values.size())
      {
        colptr.push_back(values.size());
//        n_rows = colptr.size()-1;
      }
    }

    bool _init;
    std::map <int , Scalar>    map_row;

    std::vector<Scalar>        values;
    std::vector<int>           colptr;
    std::vector<int>           rowind;

  };
};
#endif
