/*
   Copyright (C) 1997,1998,1999,2000
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

#ifndef GVECTOR_H
#define GVECTOR_H

#include <VmUtil.h>
#include <GVector_.h>
#include <GMatrix_.h>

VM_BEGIN_NS

template<class T> inline void
GVectorT<T>::LUDBackSolve(const GMatrixT<T> &LU, const GVectorT<T> &b,
                          const GVectorT<T> &permutation)
{
    // not alias-safe with b and this!
    // note: this is from William H. Press et.al Numerical Recipes in C.
    // hiranabe modified 1-n indexing to 0-(n-1).
    // I fixed some bugs in NRC, which are bad permutation handling.
    assert(numRow == b.numRow);
    assert(numRow == LU.getNumRow());
    assert(numRow == LU.getNumCol());
        
    int n = numRow;
    T *indx = permutation.elementData();
    T *x = elem;           // to put answer.
    T *bdata = b.elementData();

    /* make permutated b (b'=Pb)*/
    int i;
    for (i = 0; i < n; ++i) {
        // not alias-safe!
        x[i] = bdata[(int)indx[i]];
    }

    /* forward substitution Ux' = b' */
    int ii = -1;
    for (i = 0; i < n; ++i) {
        T sum = x[i];
        if (0 <= ii) {
            for (int j = ii; j <= i - 1; ++j) {
                sum -= LU[i][j] * x[j];
            }
        } else if (sum != 0.0)  {
            /* found the first non-zero x */
            ii = i;
        }
        x[i] = sum;
    }

    /* backward substitution, solve x' */
    for (i = n - 1; i >= 0 ; --i) {
        T sum = x[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= LU[i][j] * x[j];
        }
        // zero-div may occur
        x[i] = sum / LU.getElement(i, i);
    }
}

template<class T> inline void
GVectorT<T>::SVDBackSolve(const GMatrixT<T> &U,
                          const GMatrixT<T> &W,
                          const GMatrixT<T> &V,
                          const GVectorT<T> &b)
{
    assert(numRow == U.getNumRow());
    assert(numRow == U.getNumCol());
    assert(numRow == W.getNumRow());
    assert(b.numRow == W.getNumCol());
    assert(b.numRow == V.getNumRow());
    assert(b.numRow == V.getNumCol());
    
    const int m = U.getNumRow();  // this.numRow
    const int n = V.getNumRow();  // b.numRow
    T *tmp = new T[n];

    int j;
    for (j = 0; j < n; ++j) {
        T s = 0.0; 
        T wj = W[j][j];
        if (wj != 0.0) {
            for (int i = 0; i < m; ++i) {
                s += U[i][j] * b[i];
            }
            s /= wj;
        }
        tmp[j] = s;
    }
    for (j = 0; j < n; ++j) {
        T s = 0.0;
        for (int jj = 0; jj < n; ++jj) {
            s += V[j][jj] * tmp[jj];
        }
        elem[j] = s;
    }
    delete[] tmp;
}

template<class T> inline bool
GVectorT<T>::equals(const GMatrixT<T> &m) const
{
    return equals(m[0]);
}

template<class T> inline void
GVectorT<T>::mul(const GMatrixT<T> &m, const GVectorT<T> &v)
{
    assert(v.numRow == m.getNumCol());
    assert(numRow == m.getNumRow());
        
    GVectorT<T> r(m.getNumRow());
    const int mRow = m.getNumRow();
    const int mCol = m.getNumCol();
    const T *mSrc = m[0];
    for (int i = 0; i < mRow; ++i) {
        const T *vSrc = v.elem;
        T &sum = r[i];
        for (int j = 0; j < mCol; ++j, ++mSrc, ++vSrc) {
            sum += (*mSrc) * (*vSrc);
        }
    }
    replace(r);
}

template<class T> inline void
GVectorT<T>::mul(const GVectorT<T> &v, const GMatrixT<T> &m)
{
    assert(v.numRow == m.getNumRow());
    assert(numRow == m.getNumCol());
    GVectorT<T> r(v.getSize());
    const int mRow = m.getNumRow();
    const int mCol = m.getNumCol();
    const T *mSrc = m[0];
    const T *vSrc = v.elem;
    for (int i = 0; i < mRow; ++i, ++vSrc) {
        const T vi = *vSrc;
        T &sum = r[i];
        for (int j = 0; j < mCol; ++j, ++mSrc) {
            sum += (*mSrc) * vi;
        }
    }
    replace(r);
}

VM_END_NS


#ifdef VM_INCLUDE_IO
template<class T> inline
VM_IOSTREAM_STD::ostream& operator<<(VM_IOSTREAM_STD::ostream& o,
                                     const VM_VECMATH_NS::GVectorT<T> &v)
{
    o << "(";
    const int len = v.getSize() - 1;
    for (int i = 0; i < len; ++i) {
        o << v[i] << ",";
    }
    return o << v[len] << ")";
}
#endif /* VM_INCLUDE_IO */


VM_BEGIN_NS

typedef GVectorT<double> GVector;
typedef GVectorT<float> GVectorf; // not defined in standard java3d.vecmath


#ifdef VM_INCLUDE_TOSTRING
template<class T> inline
VM_STRING_STD::string GVectorT<T>::toString(void) const {
    VM_TOSTRING
}
#endif

VM_END_NS

#endif /* GVECTOR_H */

// Local Variables:
// mode: c++
// End:
