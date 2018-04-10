/*
   Copyright (C) 1997,1998,1999,2000,2002
   Kenji Hiranabe, Eiwa System Management, Inc.
   Katsuaki Kawachi, The University of Tokyo

   This program is free software.
   Implemented by Kenji Hiranabe(hiranabe@esm.co.jp) and
                  Katsuaki Kawachi(kawachi@cim.pe.u-tokyo.ac.jp)
   conforming to the Java(TM) 3D API specification by Sun Microsystems.

   Permission to use, copy, modify, distribute and sell this software
   and its documentation for any purpose is hereby granted without fee,
   provided that the above copyright notice appear in all copies and
   that both that copyright notice and this permission notice appear
   in supporting documentation. Kenji Hiranabe and Eiwa System Management,Inc.
   makes no representations about the suitability of this software for any
   purpose.  It is provided "AS IS" with NO WARRANTY.
*/

//
// $Id:$
//

#ifndef GMATRIX_H
#define GMATRIX_H

#include <VmUtil.h>
#include <GMatrix_.h>
#include <GVector.h>

VM_BEGIN_NS

template<class T> inline void
GMatrixT<T>::getColumn(int col, GVectorT<T> &v) const
{
    assert(numCol > col);
    assert(col >= 0);
    assert(v.getSize() >= numRow);
    T *src = elem + col;
    T *dst = v.elementData();
    for (int i = 0; i < numRow; ++i, ++dst, src += numCol) {
        *dst = *src;
    }
}

template<class T> inline void
GMatrixT<T>::getRow(int row, GVectorT<T> &v) const
{
    assert(numRow > row);
    assert(row >= 0);
    assert(v.getSize() >= numCol);
    T *src = (*this)[row];
    T *dst = v.elementData();
    memcpy(dst, src, sizeof(T) * numCol);
}

template<class T> inline void
GMatrixT<T>::invert(void)
{
    assert(numRow == numCol);
    const int n = numRow;
    GMatrixT<T> LU(n, n);
    GVectorT<T> permutation(n);
    GVectorT<T> column(n);
    GVectorT<T> unit(n);
    LUD(LU, permutation);
    for (int j = 0; j < n; ++j) {
        unit.zero();
        unit.setElement(j, 1.0);
        column.LUDBackSolve(LU, unit, permutation);
        setColumn(j, column);
    }
}

template<class T> inline int
GMatrixT<T>::LUD(GMatrixT<T> &LU, GVectorT<T> &permutation) const
{
    assert(numRow == numCol);
    assert(numRow == LU.numRow);
    assert(numRow == LU.numCol);
    assert(permutation.getSize() >= numRow);
            
    const int n = numRow;
    if (this != &LU) {
        LU.set(*this);
    }

    int even = 1;   // permutation Odd/Even

    // initialize index
    for (int i = 0; i < n; ++i) {
        permutation[i] = T(i);
    }

    // start Crout's method
    for (int j = 0; j < n; ++j) {
        T big, dum, sum;
        int imax;           // the pivot row number

        // upper portion (U)
        int i;
        for (i = 0; i < j; ++i) {
            sum = LU[i][j];
            for (int k = 0; k < i; ++k) {
                if (LU[i][k] != 0.0 && LU[k][j] != 0.0) {
                    sum -= LU[i][k] * LU[k][j];
                }
            }
            LU[i][j] = sum;
        }
        big = 0.0;
        imax = j;

        // lower part (L)
        for (i = j; i < n; ++i) {
            sum = LU[i][j];
            for (int k = 0; k < j; ++k) {
                if (LU[i][k] != 0.0 && LU[k][j] != 0.0)
                    sum -= LU[i][k] * LU[k][j];
            }
            LU[i][j] = sum;
            dum = abs(sum);
            if (dum >= big) {
                big = dum;
                imax = i;   // imax is the pivot
            }
        }

        if (j != imax) {    // if pivot is not on the diagonal
            // swap rows
            LU.swapRows(imax, j);
            T tmp = permutation[imax];
            permutation[imax] = permutation[j];
            permutation[j] = tmp;
            even = -even;
        }

        // zero-div occurs.
        // if (a[j][j] == 0.0) 

        if (j != n - 1) {
            dum = 1.0 / LU[j][j];
            for (int i = j + 1; i < n; ++i) {
                LU[i][j] *= dum;
            }
        }

    } // end of for j
    return even;
}
    
template<class T> inline void
GMatrixT<T>::mul(const GVectorT<T> &v1, const GVectorT<T> &v2)
{
    assert(numRow >= v1.getSize());
    assert(numCol >= v2.getSize());
    
    T *dst = elem;
    const T *src1 = v1.elementData();
    for (int i = 0; i < numRow; ++i, ++src1) {
        const T *src2 = v2.elementData();
        for (int j = 0; j < numCol; ++j, ++dst, ++src2) {
            *dst = (*src1) * (*src2);
        }
    }
}

template<class T> inline void
GMatrixT<T>::setColumn(int col, const GVectorT<T> &v)
{
    assert(numCol > col);
    assert(col >= 0);
    assert(v.getSize() >= numRow);
    T *dst = elem + col;
    const T *src = v.elementData();
    for (int i = 0; i < numRow; ++i, ++src, dst += numCol) {
        *dst = *src;
    }
}

template<class T> inline void
GMatrixT<T>::setRow(int row, const GVectorT<T> &v)
{
    assert(numRow > row);
    assert(row >= 0);
    assert(v.getSize() >= numCol);
    T *dst = (*this)[row];
    const T *src = v.elementData();
    memcpy(dst, src, sizeof(T) * numCol);
}

VM_END_NS


#ifdef VM_INCLUDE_IO
template<class T> inline
VM_IOSTREAM_STD::ostream& operator<<(VM_IOSTREAM_STD::ostream& o,
                                     const VM_VECMATH_NS::GMatrixT<T> &m)
{
    o << "[";
    const int numRow = m.getNumRow();
    for (int i = 0; i < numRow; ++i) {
        o << "[";
        const int numCol = m.getNumCol() - 1;
        for (int j = 0; j < numCol; ++j) {
            o << m[i][j] << ",";
        }
        o << m[i][numCol] << "]";
        if (i < numRow - 1) {
            o << "," << VM_IOSTREAM_STD::endl;
        }
    }
    return o << "]";
}
#endif /* VM_INCLUDE_IO */


VM_BEGIN_NS

typedef GMatrixT<double> GMatrix;
typedef GMatrixT<float> GMatrixf; // not defined in standard java3d.vecmath


#ifdef VM_INCLUDE_TOSTRING
template<class T> inline VM_STRING_STD::string GMatrixT<T>::toString() const {
    VM_TOSTRING
}
#endif

VM_END_NS

#endif /* GMATRIX_H */

// Local Variables:
// mode: c++
// End:
