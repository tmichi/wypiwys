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

#ifndef GVECTOR__H
#define GVECTOR__H

#include <Tuple2.h>
#include <Tuple3.h>
#include <Tuple4.h>

VM_BEGIN_NS

/**
 * A double/single precision, general, and dynamically resizeable one
 * dimensional vector class.  Index numbering begins with zero.
 */

template<class T> class GMatrixT;


template<class T>
class GVectorT {
protected:
    static T abs(T t) {return VmUtil<T>::abs(t);}
    int numRow;
    T *elem;

    void resize(int len) {
        static GMVallocator<T> galloc;
        if (len > 0) {
            T *newElem = galloc.allocate(len);
            const int minRow = numRow < len ? numRow : len;
            memcpy(newElem, elem, sizeof(T) * minRow);
            galloc.deallocate(elem, numRow);
            elem = newElem;
        } else {
            galloc.deallocate(elem, numRow);
            elem = NULL;
        }
        numRow = len;
    }

    void replace(GVectorT<T> &r) {
        static GMVallocator<T> galloc;
        galloc.deallocate(elem, numRow);
        elem = r.elem;
        numRow = r.numRow;
        r.elem = NULL;
        r.numRow = 0;
    }
    
    void create(int s) {
        static GMVallocator<T> galloc;
        assert(s >= 0);
        numRow = s;
        elem = (s > 0) ? galloc.allocate(s) : 0;
    }
    
public:
    T operator[](int i) const {
        assert(0 <= i);
        assert(i < numRow);
        return elem[i];
    }
    T& operator[](int i) {
        assert(0 <= i);
        assert(i < numRow);
        return elem[i];
    }
    T *elementData(void) const {return elem;}

    GVectorT<T>& operator=(const GVectorT<T> &v) {
        set(v);
        return *this;
    }

    ~GVectorT<T>() {
        static GMVallocator<T> galloc;
        galloc.deallocate(elem, numRow);
    }

    GVectorT<T>(void) : numRow(0) {create(0);}
    
    /**
     * Constructs a new generalized mathematic Vector with zero 
     * elements; length reprents the number of elements in the 
     * vector. !! this comment is a bug in Sun's API !!
     * @param v the values for the new vector. 
     * @param l number of elements in the vector v
     */
    GVectorT<T>(const T v[], int l) {
        create(l);
        for (int i = 0; i < l; ++i) {elem[i] = v[i];}
    }
    /**
     * Constructs a new GVector and copies the initial values from 
     * the parameter vector.  
     * @param   v the source for the new GVector's initial values 
     */
    GVectorT<T>(const GVectorT<T> &v) {create(v.getSize()); set(v);}
    
    /**
     * Constructs a new generalized mathematic Vector with zero
     * elements; length reprents the number of elements in the
     * vector.
     * @param l number of elements in this vector. 
     */
    GVectorT<T>(int l) {create(l); zero();}
    
    /**
     * Constructs a new GVector and copies the initial values from 
     * the Tuple 
     * @param   t the source for the new GVector's initial values 
     */
    GVectorT<T>(const Tuple2f &t) {create(2); set(t);}
    
    /**
     * Constructs a new GVector and copies the initial values from 
     * the Tuple 
     * @param   t the source for the new GVector's initial values 
     */
    GVectorT<T>(const Tuple2d &t) {create(2); set(t);}
    
    /**
     * Constructs a new GVector and copies the initial values from 
     * the Tuple 
     * @param   t the source for the new GVector's initial values 
     */
    GVectorT<T>(const Tuple3f &t) {create(3); set(t);}
    
    /**
     * Constructs a new GVector and copies the initial values from 
     * the Tuple 
     * @param   t the source for the new GVector's initial values 
     */
    GVectorT<T>(const Tuple3d &t) {create(3); set(t);}

    /**
     * Constructs a new GVector and copies the initial values from 
     * the Tuple 
     * @param   t the source for the new GVector's initial values 
     */
    GVectorT<T>(const Tuple4f &t) {create(4); set(t);}

    /**
     * Constructs a new GVector and copies the initial values from 
     * the Tuple 
     * @param   t the source for the new GVector's initial values 
     */
    GVectorT<T>(const Tuple4d &t) {create(4); set(t);}

    /**
     * Sets the value of this vector to sum of itself and the 
     * specified vector 
     * @param   v the second vector 
     */
    void add(const GVectorT<T> &v) {
        assert(numRow == v.numRow);
        const T *src = v.elementData();
        T *dst = elem;
        for (int i = 0; i < numRow; ++i, ++dst, ++src) {
            *dst += *src;
        }
    }

    /**
     * Sets the value of this vector to the vector sum of vectors 
     * vector1 and vector2.
     * @param   v1 the first vector 
     * @param   v2 the second vector 
     */
    void add(const GVectorT<T> &v1, const GVectorT<T> &v2) {
        assert(v1.numRow == numRow);
        assert(v1.numRow == v2.numRow);
        const T *src1 = v1.elementData();
        const T *src2 = v2.elementData();
        T *dst = elem;
        for (int i = 0; i < numRow; ++i, ++dst, ++src1, ++src2) {
            *dst = (*src1) + (*src2);
        }
    }

    /**
     * Returns the (n-space) angle in radians between this vector 
     * and the vector parameter; the return value is constrained to 
     * the range [0,PI].
     * @param   v The other vector 
     * @return  The angle in radians in the range [0,PI] 
     */
    T angle(const GVectorT<T> &v) const {
        return acos(dot(v) / norm() / v.norm());
    }

    /**
     * Returns the dot product of this vector and vector v.
     * @param   v the other vector 
     * @return  the dot product of this and v 
     */
    T dot(const GVectorT<T> &v) const {
        assert(numRow == v.numRow);
        const T *src1 = elem;
        const T *src2 = v.elementData();
        T sum = 0.0;
        for (int i = 0; i < numRow; ++i, ++src1, ++src2) {
            sum += (*src1) * (*src2);
        }
        return sum;
    }

    /**
     * Returns true if the L-infinite distance between this vector 
     * and vector v is less than or equal to the epsilon parameter, 
     * otherwise returns false. The L-infinite distance is equal 
     * to MAX[abs(x1-x2), abs(y1-y2), . . . ]. <p>
     * @param   v The vector to be compared to this vector 
     * @param   epsilon the threshold value 
     */
    bool epsilonEquals(const GVectorT<T> &v, T eps) const {
        if (numRow != v.getSize()) {
            return false;
        }
        const T *src1 = elem;
        const T *src2 = v.elementData();
        for (int i = 0; i < numRow; ++i, ++src1, ++src2) {
            if (abs((*src1) - (*src2)) > eps) {
                return false;
            }
        }
        return true;
    }

    /**
     * Returns true if all of the data members of GVector v 
     * are equal to the corresponding data members in this GVector.
     * @param   v The vector with which the comparison is made. 
     * @return  true or false 
     */
    bool equlas(const GVectorT<T> &v) const {
        if (numRow != v.getSize()) {
            return false;
        }
        const T *src1 = elem;
        const T *src2 = v.elementData();
        for (int i = 0; i < numRow; ++i, ++src1, ++src2) {
            if (*src1 != *src2) {
                return false;
            }
        }
        return true;
    }

    /**
      * Returns true if the first column of matrix m are equal to
      * the corresponding data members in this GVector.
      * @param m the matrix with which the comparison is made.
      */
    bool equals(const GMatrixT<T> &m) const;

    /**
      * Returns true if all of the data members of m are equal to the
      * corresponding data members in this GVector.
      * @param m the array with which the comparison is made.
      */
    bool equals(const T *m) const {
        const T *src1 = elem;
        const T *src2 = m;
        for (int i = 0; i < numRow; ++i, ++src1, ++src2) {
            if (*src1 != *src2) {
                return false;
            }
        }
        return true;
    }

    /**
     * Retrieves the value at the specified index value of this 
     * vector.
     * @param   i the index of the element to retrieve (zero indexed) 
     * @return  the value at the indexed element 
     */
    T getElement(int i) const {
        assert(0 <= i);
        assert(i < numRow);
        return elem[i];
    }
    
    /**
     * Returns the number of elements in this vector. 
     * @return  number of elements in this vector 
     */
    int getSize(void) const {return numRow;}
    
    /**
     * Returns a hash number based on the data values in this object. 
     * @return the integer hash value
     */
    int hashCode(void) const {
        return VmUtil<T>::hashCode(sizeof *this, this);
    }

    /**
     * Linearly interpolates between this vector and vector v1 and 
     * places the result into this tuple: 
     * this = (1-alpha)*this + alpha*v1. 
     * @deprecated the double version of this method should be used.
     * @param   v1 the first vector 
     * @param   a the alpha interpolation parameter 
     */
    void interpolate(const GVectorT<T> &v, T a) {
        assert(numRow == v.numRow);
        T *dst = elem;
        const T *src = v.elementData();
        T b = 1.0 - a;
        for (int i = 0; i < numRow; ++i, ++dst, ++src) {
            *dst = b * (*dst) + a * (*src);
        }
    }

    /**
     * Linearly interpolates between vectors v1 and v2 and places 
     * the result into this tuple: this = (1-alpha)*v1 + alpha*v2.
     * @deprecated the double version of this method should be used.
     * @param   v1 the first vector 
     * @param   v2 the second vector 
     * @param   a the alpha interpolation parameter 
     */
    void interpolate(const GVectorT<T> &v1, const GVectorT<T> &v2, T a) {
        set(v1);
        interpolate(v2, a);
    }

    /**
     * LU Decomposition Back Solve; this method takes the LU matrix 
     * and the permutation vector produced by the GMatrix method LUD 
     * and solves the equation (LU)*x = b by placing the solution 
     * vector x into this vector. This vector should be the same 
     * length or longer than b. 
     *
     * @param   LU The matrix into which the lower and upper 
     *          decompositions have been placed 
     * @param   b The b vector in the equation (LU)*x = b 
     * @param   permutation The row permuations that were necessary 
     *          to produce the LU matrix parameter 
     */
    inline void LUDBackSolve(const GMatrixT<T> &LU, const GVectorT<T> &b,
                             const GVectorT<T> &permutation);
    
    /**
     * Multiplies matrix m times Vector v and places the result 
     * into this vector (this = m*v).
     * @param   m The matrix in the multiplication 
     * @param   v The vector that is multiplied 
     */
    inline void mul(const GMatrixT<T> &m, const GVectorT<T> &v);
    
    /**
     * Multiplies the transpose of vector v (ie, v becomes a row 
     * vector with respect to the multiplication) times matrix m 
     * and places the result into this vector 
     * (this = transpose(v)*m). The result is technically a row 
     * vector, but the GVector class only knows about column vectors, 
     * and so the result is stored as a column vector. 
     * @param   m The matrix in the multiplication 
     * @param   v The vector that is temporarily transposed 
     */
    inline void mul(const GVectorT<T> &v, const GMatrixT<T> &m);
    
    /**
     * Negates the value of this vector: this = -this. 
     */
    void negate(void) {
        T *dst = elem;
        for (int i = 0; i < numRow; ++i, ++dst) {
            *dst = -(*dst);
        }
    }

    /**
     * Returns the square root of the sum of the squares of this 
     * vector (its length in n-dimensional space). 
     *
     * @return  length of this vector 
     */
    T norm(void) const {
        return sqrt(normSquared());
    }

    /**
     * Normalizes this vector in place. 
     */
    void normalize(void) {
        T *dst = elem;
        const T n = norm(); // zero-div may happen
        for (int i = 0; i < numRow; ++i, ++dst) {
            *dst /= n;
        }
    }

    /**
     * Sets the value of this vector to the normalization of 
     * vector v. 
     * @param   v the un-normalized vector 
     */
    void normalize(const GVectorT<T> &v) {
        T *dst = elem;
        const T *src = v.elem;
        const T n = norm(); // zero-div may happen
        for (int i = 0; i < numRow; ++i, ++dst, ++src) {
            *dst = (*src) / n;
        }
    }

    /**
     * Returns the sum of the squares of this vector (its length 
     * sqaured in n-dimensional space). <p>
     * @return  length squared of this vector 
     */
    T normSquared(void) const {
        const T *src = elem;
        T sum = 0.0;
        for (int i = 0; i < numRow; ++i, ++src) {
            sum += (*src) * (*src);
        }
        return sum;
    }

    /**
     * Scales this vector by the scale factor s. 
     * @param   s the scalar value 
     */
    void scale(T s) {
        T *src = elem;
        for (int i = 0; i < numRow; ++i, ++src) {
            *src *= s;
        }
    }

    /**
     * Sets the value of this vector to the scalar multiplication of 
     * the scale factor with the vector v1.
     * @param   s the scalar value 
     * @param   v the source vector 
     */
    void scale(T s, const GVectorT<T> &v) {
        assert(numRow == v.numRow);
        const T *src = v.elementData();
        T *dst = elem;
        for (int i = 0; i < numRow; ++i, ++src) {
            *dst = s * (*src);
        }
    }

    /** 
     * Sets the value of this vector to the scalar multiplication by 
     * s of vector v1 plus vector v2 (this = s*v1 + v2).
     * @param   s the scalar value 
     * @param   v1 the vector to be multiplied 
     * @param   v2 the vector to be added 
     */
    void scaleAdd(T s, const GVectorT<T> &v1, const GVectorT<T> &v2) {
        assert(numRow == v1.numRow);
        assert(numRow == v2.numRow);
        const T *src1 = v1.elementData();
        const T *src2 = v2.elementData();
        T *dst = elem;
        for (int i = 0; i < numRow; ++i, ++dst, ++src1, ++src2) {
            *dst = s * (*src1) + (*src2);
        }
    }

    /**
     * Sets the value of this vector to the values found in the 
     * array parameter. The array should be at least equal in 
     * length to the number of elements in the vector. 
     * @param   v the source array 
     */
    void set(const GVectorT<T> &v) {
        assert(numRow == v.numRow);
        const T *src = v.elementData();
        T *dst = elem;
        for (int i = 0; i < numRow; ++i, ++dst, ++src) {
            *dst = *src;
        }
    }

    /**
     * Sets the value of this vector to the values in tuple.
     * @param   t the source for the new GVector's new values 
     */
    void set(const Tuple2f &t) {
        elem[0] = T(t.x);
        elem[1] = T(t.y);
    }
    
    /**
     * Sets the value of this vector to the values in tuple.
     * @param   t the source for the new GVector's new values 
     */
    void set(const Tuple2d &t) {
        elem[0] = T(t.x);
        elem[1] = T(t.y);
    }
    
    /**
     * Sets the value of this vector to the values in tuple.
     * @param   t the source for the new GVector's new values 
     */
    void set(const Tuple3f &t) {
        elem[0] = T(t.x);
        elem[1] = T(t.y);
        elem[2] = T(t.z);
    }
    
    /**
     * Sets the value of this vector to the values in tuple.
     * @param   t the source for the new GVector's new values 
     */
    void set(const Tuple3d &t) {
        elem[0] = T(t.x);
        elem[1] = T(t.y);
        elem[2] = T(t.z);
    }
    
    /**
     * Sets the value of this vector to the values in tuple.
     * @param   t the source for the new GVector's new values 
     */
    void set(const Tuple4f &t) {
        elem[0] = T(t.x);
        elem[1] = T(t.y);
        elem[2] = T(t.z);
        elem[3] = T(t.w);
    }
    
    /**
     * Sets the value of this vector to the values in tuple.
     * @param   t the source for the new GVector's new values 
     */
    void set(const Tuple4d &t) {
        elem[0] = T(t.x);
        elem[1] = T(t.y);
        elem[2] = T(t.z);
        elem[3] = T(t.w);
    }

    /**
     * Modifies the value at the specified index of this vector. 
     * @param   i the index if the element to modify (zero indexed) 
     * @param   val the new vector element value 
     */
    void setElement(int i, T val) {
        assert(0 <= i);
        assert(i < numRow);
        elem[i] = val;
    }

    /**
     * Changes the size of this vector dynamically. If the size is 
     * increased no data values will be lost. If the size is 
     * decreased, only those data values whose vector positions 
     * were eliminated will be lost.
     * @param   l number of desired elements in this vector 
     */
    void setSize(int l) {
        assert(l >= 0);
        const int oldSiz = numRow;
        resize(l);
        for (int i = oldSiz; i < l; ++i) {elem[i] = 0.0;}
    }

    /**
     * Sets the value of this vector to the vector difference of 
     * itself and vector (this = this - v).
     * @param   v the other vector 
     */
    void sub(const GVectorT<T> &v) {
        assert(numRow == v.numRow);
        const T *src = v.elementData();
        T *dst = elem;
        for (int i = 0; i < numRow; ++i, ++dst, ++src) {
            *elem -= (*src);
        }
    }

    /**
     * Sets the value of this vector to the vector sum of vectors 
     * vector1 and vector2.
     * @param   v1 the first vector 
     * @param   v2 the second vector 
     */
    void sub(const GVectorT<T> &v1, const GVectorT<T> &v2) {
        assert(v1.numRow == numRow);
        assert(v1.numRow == v2.numRow);
        const T *src1 = v1.elementData();
        const T *src2 = v2.elementData();
        T *dst = elem;
        for (int i = 0; i < numRow; ++i, ++dst, ++src1, ++src2) {
            *elem = (*src1) - (*src2);
        }
    }

    /**
     * Solves for x in Ax = b, where x is this vector (nx1), 
     * A is mxn, b is mx1, and A = U*W*transpose(V); 
     * U,W,V must be precomputed and can be found by taking the 
     * singular value decomposition (SVD) of A using the method 
     * SVD found in the GMatrix class.
     * @param   U The U matrix produced by the GMatrix method SVD 
     * @param   W The W matrix produced by the GMatrix method SVD 
     * @param   V The V matrix produced by the GMatrix method SVD 
     * @param   b The b vector in the linear equation Ax = b 
     */
    inline void SVDBackSolve(const GMatrixT<T> &U, const GMatrixT<T> &W,
                             const GMatrixT<T> &V, const GVectorT<T> &b);

    /**
     * Sets all the values in this vector to zero. 
     */
    void zero(void) {
        T *dst = elem;
        for (int i = 0; i < numRow; ++i, ++dst) {
            *dst = 0.0;
        }
    }
    
#ifdef VM_INCLUDE_TOSTRING
    inline VM_STRING_STD::string toString(void) const;
#endif
};

VM_END_NS

#endif /* GVECTOR__H */

// Local Variables:
// mode: c++
// End:
