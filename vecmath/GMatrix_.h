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

#ifndef GMATRIX__H
#define GMATRIX__H

#include <Matrix3.h>
#include <Matrix4.h>
#ifdef VM_STD_C_HEADERS
# include <cfloat>
#else
# include <float.h>
#endif

VM_BEGIN_NS

/**
 * A double/single precision, general, real, and dynamically
 * resizeable two dimensional N x M matrix class.  Row and column
 * numbering begins with zero.  The representation is row major.
 */

template<class T> class GVectorT;


template<class T>
class GMatrixT {
protected:
    static T abs(T t) {return VmUtil<T>::abs(t);}
    
    /**
     * Number of rows in this matrix.
     */
    int numRow;

    /**
     * Number of columns in this matrix.
     */
    int numCol;

    /**
     * Number of rows x cols (i.e. number of elements)
     */
    int numRowCol;

    /**
     * Data of the GMatrix.
     * (1D array. The (i,j) element is stored in elementData[i*nCol + j])
     */
    T *elem;

    /**
     * Pointers to beginnings of rows
     */
    T **elemArray;

    void resize(int r, int c) {
        static GMVallocator<T> galloc;
        const int rc = r * c;
        if (r > 0) {
            const T *src = elem;
            T *newElem = galloc.allocate(rc);
            const int minCol = numCol < c ? numCol : c;
            const int minRow = numRow < r ? numRow : r;
            T *dst = newElem;
            const int len = sizeof(T) * minCol;
            if (len > 0) {
                int i;
                for (i = 0; i < minRow; ++i, src += numCol, dst += c) {
                    memcpy(dst, src, len);
                }
            }
            galloc.deallocate(elem, numRowCol);
            elem = newElem;

            galloc.p_deallocate(elemArray, numRow);
            elemArray = galloc.p_allocate(r);
            int i;
            if (c > 0) {
                T *p;
                for (i = 0, p = elem; i < r; ++i, p += c) {elemArray[i] = p;}
            } else {
                for (i = 0; i < r; ++i) {elemArray[i] = NULL;}
            }
        } else {
            galloc.deallocate(elem, numRowCol);
            elem = 0;
            galloc.p_deallocate(elemArray, numRow);
            elemArray = 0;
        }
        numRowCol = rc;
        numRow = r;
        numCol = c;
    }

    void create(int r, int c, const T *const m = 0) {
        static GMVallocator<T> galloc;
        const int rc = r * c;
        if (r > 0) {
            elem = galloc.allocate(rc);
            elemArray = galloc.p_allocate(r);
            int i;
            if (c > 0) {
                T *p;
                for (i = 0, p = elem; i < r; ++i, p += c) {elemArray[i] = p;}
            } else {
                for (i = 0; i < r; ++i) {elemArray[i] = NULL;}
            }
        } else {
            elem = NULL;
            elemArray = NULL;
        }
        if (m) {memcpy(elem, m, sizeof(T) * rc);}
        
        numRow = r;
        numCol = c;
        numRowCol = rc;
    }

    void capture(GMatrixT<T> &m) {
        static GMVallocator<T> galloc;
        galloc.deallocate(elem, numRowCol);
        galloc.p_deallocate(elemArray, numRow);
        
        elem = m.elem;
        elemArray = m.elemArray;
        numRow = m.numRow;
        numCol = m.numCol;
        
        m.elem = NULL;
        m.elemArray = NULL;
        m.numRow = m.numCol = m.numRowCol = 0;
    }
    
public:
    
    T *operator[](int i) const {return elemArray[i];}
    GMatrixT<T>& operator=(const GMatrixT<T> &m) {
        set(m);
        return *this;
    }

    ~GMatrixT<T>() {
        static GMVallocator<T> galloc;
        galloc.deallocate(elem, numRowCol);
        galloc.p_deallocate(elemArray, numRow);
    }

    /**
     * Constructs an nRow by nCol identity matrix. 
     * Note that even though row and column numbering begins with
     * zero, nRow and nCol will be one larger than the maximum
     * possible matrix index values.
     * @param nrow number of rows in this matrix.
     * @param ncol number of columns in this matrix.
     */
    GMatrixT<T>(int nrow, int ncol) {
        assert(nrow >= 0);
        assert(ncol >= 0);
        create(nrow, ncol);
        setIdentity();
    }

    /**
     * Constructs an nRow by nCol matrix initialized to the values 
     * in the matrix array.  The array values are copied in one row at
     * a time in row major fashion.  The array should be at least 
     * nRow*nCol in length.
     * Note that even though row and column numbering begins with 
     * zero, nRow and nCol will be one larger than the maximum
     * possible matrix index values.
     * @param nRow number of rows in this matrix.
     * @param nCol number of columns in this matrix.
     * @param matrix a 1D array that specifies a matrix in row major fashion
     */
    GMatrixT<T>(int nrow, int ncol, const T *matrix) {
        assert(nrow >= 0);
        assert(ncol >= 0);
        create(nrow, ncol, matrix);
    }
    
    /**
     * Constructs a new GMatrix and copies the initial values
     * from the parameter matrix.
     * @param m the source of the initial values of the new GMatrix
     */
    GMatrixT<T>(const GMatrixT<T> &m) {create(m.numRow, m.numCol, m.elem);}
    
    GMatrixT<T>(void) {create(0, 0);}
    
    /**
     * Sets the value of this matrix to sum of itself and matrix m.
     * @param m the other matrix
     */
    void add(const GMatrixT<T> &m) {
        assert(numRow == m.numRow);
        assert(numCol == m.numCol);
        T *src = m.elem;
        T *dst = elem;
        for (int i = 0; i < numRowCol; ++i, ++dst, ++src) {
            *dst += *src;
        }
    }
    
    /**
     * Sets the value of this matrix to the matrix sum of matrices m1 and m2.
     * @param m1 the first matrix
     * @param m2 the second matrix
     */
    void add(const GMatrixT<T> &m1, const GMatrixT<T> &m2) {
        assert(numRow == m1.numRow);
        assert(numCol == m1.numCol);
        assert(numRow == m2.numRow);
        assert(numCol == m2.numCol);
        T *src1 = m1.elem;
        T *src2 = m2.elem;
        T *dst = elem;
        for (int i = 0; i < numRowCol; ++i, ++dst, ++src1, ++src2) {
            *dst = *src1 + *src2;
        }
    }


    /**
     * Copies a sub-matrix derived from this matrix into the target matrix.
     * The upper left of the sub-matrix is located at (rowSource, colSource);
     * the lower right of the sub-matrix is located at 
     * (lastRowSource,lastColSource).  The sub-matrix is copied into the
     * the target matrix starting at (rowDest, colDest).
     * @param rowSource the top-most row of the sub-matrix
     * @param colSource the left-most column of the sub-matrix
     * @param numRow the number of rows in the sub-matrix
     * @param numCol the number of columns in the sub-matrix
     * @param rowDest the top-most row of the position of the copied sub-matrix
     *                  within the target matrix
     * @param colDest the left-most column of the position of the copied
     *                  sub-matrix within the target matrix
     * @param m the matrix into which the sub-matrix will be copied
     */
    void copySubMatrix(int rowSrc, int colSrc, int lenRowSrc, int lenColSrc,
                       int rowDst, int colDst,
                       GMatrixT<T> &m) const {
        assert(rowSrc >= 0);
        assert(colSrc >= 0);
        assert(rowDst >= 0);
        assert(colDst >= 0);
        assert(numRow >= numRow + rowSrc);
        assert(numCol >= numCol + colSrc);
        assert(m.numRow >= numRow + rowDst);
        assert(m.numCol >= numCol + colDst);
        
        T *src = (*this)[rowSrc] + colSrc;
        T *dst = m[rowDst] + colDst;
        const int len = sizeof(T) * lenColSrc;
        for(int i = 0; i < lenRowSrc; ++i, src += numCol, dst += m.numCol) {
            memcpy(dst, src, len);
        }
    }
  
    /**
     * Returns true if the L-infinite distance between this matrix and 
     * matrix m is less than or equal to the epsilon parameter,
     * otherwise returns false. The L-infinite distance is equal to 
     * MAX[i=0,1,2, . . .n ; j=0,1,2, . . .n ; abs(this.m(i,j) - m.m(i,j)] . 
     * @deprecated The double version of this method should be used.
     * @param m The matrix to be compared to this matrix 
     * @param epsilon the threshold value
     */
    bool epsilonEquals(const GMatrixT<T> &m, T eps) const {
        if (m.numRow != numRow || m.numCol != numCol) {
            return false;
        }
        T *m1 = m.elem;
        T *m2 = elem;
        for (int i = 0; i < numRowCol; ++i, ++m1, ++m2) {
            if (abs(*m1 - *m2) > eps) {
                return false;
            }
        }
        return true;
    }
    
    /**
     * Returns true if the L-infinite distance between this matrix and 
     * matrix m is less than or equal to the epsilon parameter,
     * otherwise returns false. The L-infinite distance is equal to 
     * MAX[i=0,1,2, . . .n ; j=0,1,2, . . .n ; abs(this.m(i,j) - m.m(i,j)] . 
     * @param m The matrix to be compared to this matrix 
     * @param epsilon the threshold value
     */
    bool equals(const GMatrixT<T> &m) const {
        if (m.numRow != numRow || m.numCol != numCol) {
            return false;
        }
        T *m1 = m.elem;
        T *m2 = elem;
        for (int i = 0; i < numRowCol; ++i, ++m1, ++m2) {
            if (*m1 != *m2) {
                return false;
            }
        }
        return true;
    }

    /**
     * Places the values in the this matrix into the matrix m; m
     * should be at least as large as this GMatrix.
     * @param m The matrix that will hold the new values 
     */
    void get(GMatrixT<T> &m) const {
        assert(m.numRow >= numRow);
        assert(m.numCol >= numCol);
        copySubMatrix(0, 0, numRow, numCol, 0, 0, m);
    }

    /**
     * Places the values in the upper 3X3 of this GMatrix into the matrix m.
     * @param m The matrix that will hold the new values 
     */
    void get(Matrix3d &m) const {
        assert(numRow >= 3);
        assert(numCol >= 3);
        T *src = elem;
        m.m00 = double(*src); ++src;
        m.m01 = double(*src); ++src; 
        m.m02 = double(*src);
        src = (*this)[1];
        m.m10 = double(*src); ++src;
        m.m11 = double(*src); ++src;
        m.m12 = double(*src);
        src = (*this)[2];
        m.m20 = double(*src); ++src;
        m.m21 = double(*src); ++src;
        m.m22 = double(*src);
    }

    /**
     * Places the values in the upper 3X3 of this GMatrix into the matrix m.
     * @param m The matrix that will hold the new values 
     */
    void get(Matrix3f &m) const {
        assert(numRow >= 3);
        assert(numCol >= 3);
        T *src = elem;
        m.m00 = float(*src); ++src;
        m.m01 = float(*src); ++src; 
        m.m02 = float(*src);
        src = (*this)[1];
        m.m10 = float(*src); ++src;
        m.m11 = float(*src); ++src;
        m.m12 = float(*src);
        src = (*this)[2];
        m.m20 = float(*src); ++src;
        m.m21 = float(*src); ++src;
        m.m22 = float(*src);
    }

    /**
     * Places the values in the upper 4X4 of this GMatrix into the matrix m.
     * @param m The matrix that will hold the new values 
     */
    void get(Matrix4d &m) const {
        assert(numRow >= 4);
        assert(numCol >= 4);
        T *src = elem;
        m.m00 = double(*src); ++src;
        m.m01 = double(*src); ++src; 
        m.m02 = double(*src); ++src;
        m.m03 = double(*src);
        src = (*this)[1];
        m.m10 = double(*src); ++src;
        m.m11 = double(*src); ++src;
        m.m12 = double(*src); ++src;
        m.m13 = double(*src);
        src = (*this)[2];
        m.m20 = double(*src); ++src;
        m.m21 = double(*src); ++src;
        m.m22 = double(*src); ++src;
        m.m23 = double(*src);
        src = (*this)[3];
        m.m30 = double(*src); ++src;
        m.m31 = double(*src); ++src;
        m.m32 = double(*src); ++src;
        m.m33 = double(*src);
    }

    /**
     * Places the values in the upper 4X4 of this GMatrix into the matrix m.
     * @param m The matrix that will hold the new values 
     */
    void get(Matrix4f &m) const {
        assert(numRow >= 4);
        assert(numCol >= 4);
        T *src = elem;
        m.m00 = float(*src); ++src;
        m.m01 = float(*src); ++src; 
        m.m02 = float(*src); ++src;
        m.m03 = float(*src);
        src = (*this)[1];
        m.m10 = float(*src); ++src;
        m.m11 = float(*src); ++src;
        m.m12 = float(*src); ++src;
        m.m13 = float(*src);
        src = (*this)[2];
        m.m20 = float(*src); ++src;
        m.m21 = float(*src); ++src;
        m.m22 = float(*src); ++src;
        m.m23 = float(*src);
        src = (*this)[3];
        m.m30 = float(*src); ++src;
        m.m31 = float(*src); ++src;
        m.m32 = float(*src); ++src;
        m.m33 = float(*src);
    }

    /**
     * Places the values of the specified column into the array parameter. 
     * @param col the target column number
     * @param dst the array into which the column values will be placed
     */
    void getColumn(int col, T *dst) const {
        assert(numCol > col);
        assert(col >= 0);
        T *src = elem + col;
        for (int i = 0; i < numRow; ++i, ++dst, src += numCol) {
            *dst = *src;
        }
    }

    /**
     * Places the values of the specified column into the vector parameter.
     * @param col the target column number
     * @v vector the vector into which the column values will be placed
     */
    inline void getColumn(int col, GVectorT<T> &v) const;

    T getElement(int row, int col) const {
        assert(numRow > row);
        assert(row >= 0);
        assert(numCol > col);
        assert(col >= 0);
        return (*this)[row][col];
    }
    
    /**
     * Returns the number of colmuns in this matrix.
     * @return number of columns in this matrix
     */
    int getNumCol(void) const {return numCol;}
    
    /**
     * Returns the number of rows in this matrix.
     * @return number of rows in this matrix
     */
    int getNumRow(void) const {return numRow;}

    
    /**
     * Places the values of the specified row into the array parameter. 
     * @row row the target row number
     * @dst array the array into which the row values will be placed
     */
    void getRow(int row, T *dst) const {
        assert(numRow > row);
        assert(row >= 0);
        T *src = (*this)[row];
        memcpy(dst, src, sizeof(T) * numCol);
    }

    /**
     * Places the values of the specified row into the vector parameter.
     * @row row the target row number
     * @v vector the vector into which the row values will be placed
     */
    inline void getRow(int row, GVectorT<T> &v) const;

    /**
     * Returns a hash number based on the data values in this
     * object.  Two different GMatrix objects with identical data values
     * (ie, returns true for equals(GMatrix) ) will return the same hash
     * number.  Two objects with different data members may return the
     * same hash value, although this is not likely.
     * @return the integer hash value
     */
    int hashCode(void) const {
        return VmUtil<T>::hashCode(sizeof *this, this);
    }

    /**
     * Subtracts this matrix from the identity matrix and puts the values
     * back into this (this = I - this).
     */
    void identityMinus(void) {
        negate();
        const int min = numRow < numCol ? numRow : numCol;
        T *dst = elem;
        for (int i = 0; i < min; ++i) {
            *dst += 1.0;
            dst += numCol + 1;
        }
    }

    /**
     * Inverts this matrix in place. 
     */
    inline void invert(void);
    
    /**
     * Inverts matrix m and places the new values into this matrix.  Matrix
     * m is not modified. 
     * @param m the matrix to be inverted
     */
    void invert(const GMatrixT<T> &m) {
        set(m);
        invert();
    }

protected:
    void swapRows(int r1, int r2) {
        T *src1 = (*this)[r1];
        T *src2 = (*this)[r2];
        for (int k = 0; k < numCol; ++k, ++src1, ++src2) {
            T tmp = *src1;
            *src1 = *src2;
            *src2 = tmp;
        }
    }
    
public:

    /**
     * LU Decomposition; this matrix must be a square matrix; the LU GMatrix 
     * parameter must be the same size as this matrix. The matrix LU will be 
     * overwritten as the combination of a lower diagonal and upper diagonal
     * matrix decompostion of this matrix; the diagonal elements of L (unity)
     * are not stored. The GVector parameter records the row permutation 
     * effected by the partial pivoting, and is used as a parameter to the
     * GVector method LUDBackSolve to solve sets of linear equations. This 
     * method returns +/- 1 depending on whether the number of row interchanges
     * was even or odd, respectively. 
     * @param permutation The row permutation effected by the 
     * partial pivoting 
     * @return +-1 depending on whether the number of row interchanges 
     * was even or odd respectively 
     */
    inline int LUD(GMatrixT<T> &LU, GVectorT<T> &permutation) const;
    
    /**
     * Sets the value of this matrix to the result of multiplying itself
     * with matrix m (this = this * m. 
     * @param m the other matrix
     */
    void mul(const GMatrixT<T> &m) {
        mul(*this, m);
    }

    /**
     * Sets the value of this matrix to the result of multiplying
     * the two argument matrices together (this = m1 * m2).
     * @param m1 the first matrix
     * @param m2 the second matrix
     */
    void mul(const GMatrixT<T> &m1, const GMatrixT<T> &m2) {
        assert(numRow == m1.numRow);
        assert(numCol == m2.numCol);
        assert(m1.numCol == m2.numRow);
        
        bool aliasing = (&m1 == this || &m2 == this);
        GMatrixT<T> *m;
        T *dst;
        if (aliasing) {
            m = new GMatrixT<T>(numRow, numCol);
            dst = m->elem;
        } else {
            dst = elem;
        }
        const T *const bsrc2 = m2.elem;
        for (int i = 0; i < numRow; ++i) {
            const T *const bsrc1 = m1[i];
            for (int j = 0; j < numCol; ++j, ++dst) {
                T sum = 0.0;
                const T *src1 = bsrc1;
                const T *src2 = bsrc2 + j;
                for (int k = 0; k < m1.numCol; ++k,++src1,src2 += m2.numCol) {
                    sum += (*src1) * (*src2);
                }
                *dst = sum;
            }
        }
        if (aliasing) {
            capture(*m);
            delete m;
        }
    }
    
    /**
     * Computes the outer product of the two vectors; multiplies the 
     * the first vector by the transpose of the second vector
     * and places the matrix result into this matrix. This matrix must
     * be as big or bigger than getSize(v1)xgetSize(v2).
     * @param v1 the first vector, treated as a row vector 
     * @param v2 the second vector, treated as a column vector
     */
    inline void mul(const GVectorT<T> &v1, const GVectorT<T> &v2);
    
    /**
     * Multiplies the transpose of matrix m1 times the transpose of
     * matrix m2, and places the result into this.
     * @param m1 The matrix on the left hand side of the multiplication
     * @param m2 The matrix on the right hand side of the multiplication
     */
    void mulTransposeBoth(const GMatrixT<T> &m1, const GMatrixT<T> &m2) {
        assert(numRow == m1.numCol);
        assert(numCol == m2.numRow);
        assert(m1.numRow == m2.numCol);
        
        bool aliasing = (&m1 == this || &m2 == this);
        GMatrixT<T> *m;
        T *dst;
        if (aliasing) {
            m = new GMatrixT<T>(numRow, numCol);
            dst = m->elem;
        } else {
            dst = elem;
        }
        for (int i = 0; i < numRow; ++i) {
            const T *const bsrc1 = m1.elem + i;
            for (int j = 0; j < numCol; ++j, ++dst) {
                T sum = 0.0;
                const T *src1 = bsrc1;
                const T *src2 = m2[j];
                for (int k = 0; k < m1.numCol; ++k, src1 += m1.numCol, ++src2){
                    sum += (*src1) * (*src2);
                }
                *dst = sum;
            }
        }
        if (aliasing) {
            capture(*m);
            delete m;
        }
    }
    

    /**
     * Multiplies the transpose of matrix m1 times the matrix m2, and
     * places the result into this.
     * @param m1 The matrix on the left hand side of the multiplication
     * @param m2 The matrix on the right hand side of the multiplication
     */
    void mulTransposeLeft(const GMatrixT<T> &m1, const GMatrixT<T> &m2) {
        assert(numRow == m1.numCol);
        assert(numCol == m2.numCol);
        assert(m1.numRow == m2.numRow);
         
        bool aliasing = (&m1 == this || &m2 == this);
        GMatrixT<T> *m;
        T *dst;
        if (aliasing) {
            m = new GMatrixT<T>(numRow, numCol);
            dst = m->elem;
        } else {
            dst = elem;
        }
        const T *const bsrc2 = m2.elem;
        for (int i = 0; i < numRow; ++i) {
            const T *const bsrc1 = m1.elem + i;
            for (int j = 0; j < numCol; ++j, ++dst) {
                T sum = 0.0;
                const T *src1 = bsrc1;
                const T *src2 = bsrc2 + j;
                for (int k = 0; k < m1.numCol;
                     ++k, src1 += m1.numCol, src2 += m2.numCol) {
                    sum += (*src1) * (*src2);
                }
                *dst = sum;
            }
        }
        if (aliasing) {
            capture(*m);
            delete m;
        }
    }
    
    /**
     * Multiplies matrix m1 times the transpose of matrix m2, and places the
     * result into this. 
     */
    void mulTransposeRight(const GMatrixT<T> &m1, const GMatrixT<T> &m2) {
        assert(numRow == m1.numRow);
        assert(numCol == m2.numRow);
        assert(m1.numCol == m2.numCol);
        bool aliasing = (&m1 == this || &m2 == this);
        GMatrixT<T> *m;
        T *dst;
        if (aliasing) {
            m = new GMatrixT<T>(numRow, numCol);
            dst = m->elem;
        } else {
            dst = elem;
        }
        for (int i = 0; i < numRow; ++i) {
            const T *const bsrc1 = m1[i];
            for (int j = 0; j < numCol; ++j, ++dst) {
                T sum = 0.0;
                const T *src1 = bsrc1;
                const T *src2 = m2[j];
                for (int k = 0; k < numCol; ++k, ++src1, ++src2) {
                    sum += (*src1) * (*src2);
                }
                *dst = sum;
            }
        }
        if (aliasing) {
            capture(*m);
            delete m;
        }
    }
    
    /**
     * Negates the value of this matrix: this = -this.
     */
    void negate(void) {
        T *dst = elem;
        for (int i = 0; i < numRowCol; ++i, ++dst) {
            *dst = -(*dst);
        }
    }
    
    /**
     * Sets the value of this matrix to the negation of the GMatrix parameter.
     * @param m The source matrix
     */
    void negate(const GMatrixT<T> &m) {
        const T *src = m.elem;
        T *dst = elem;
        for (int i = 0; i < numRowCol; ++i, ++dst, ++src) {
            *dst = -(*src);
        }
    }

    /**
     * Sets the value of this matrix to the values found in the array
     * parameter.The values are copied in one row at a time, in row major 
     * fashion.  The array should be at least equal in length to
     * the number of matrix rows times the number of matrix columns
     * in this matrix.
     * @param m the row major source array
     */
    void set(const T *const m) {
        T *dst = elem;
        memcpy(dst, m, numRowCol * sizeof(T));
    }
    
    /**
     * Sets the value of this matrix to the values found in matrix m.
     * @param m the source matrix
     */
    void set(const GMatrixT<T> &m) {
        assert(m.numRow <= numRow);
        assert(m.numCol <= numCol);
        const T *src = m.elem;
        T *dst = elem;
        if (numCol == m.numCol) {
            if (numRow == m.numRow) {
                memcpy(dst, src, numRowCol * sizeof(T));
            } else {
                memcpy(dst, src, numCol * m.numRow * sizeof(T));
            }
        } else {
            const int len = m.numCol * sizeof(T);
            for (int i = 0; i < m.numRow; ++i, src += m.numCol, dst += numCol){
                memcpy(dst, src, len);
            }
        }
    }

    /**
     * Sets the value of this matrix to that of the Matrix3d provided. 
     * @param m the source matrix
     */
    void set(const Matrix3d &m) {
        assert(numRow >= 3);
        assert(numCol >= 3);
        T *dst = elem;
        *dst = T(m.m00); ++dst;
        *dst = T(m.m01); ++dst; 
        *dst = T(m.m02);
        dst = (*this)[1];
        *dst = T(m.m10); ++dst;
        *dst = T(m.m11); ++dst;
        *dst = T(m.m12);
        dst = (*this)[2];
        *dst = T(m.m20); ++dst;
        *dst = T(m.m21); ++dst;
        *dst = T(m.m22);
    }

    /**
     * Sets the value of this matrix to that of the Matrix3f provided. 
     * @param m the source matrix
     */
    void set(const Matrix3f &m) {
        assert(numRow >= 3);
        assert(numCol >= 3);
        T *dst = elem;
        *dst = T(m.m00); ++dst;
        *dst = T(m.m01); ++dst; 
        *dst = T(m.m02);
        dst = (*this)[1];
        *dst = T(m.m10); ++dst;
        *dst = T(m.m11); ++dst;
        *dst = T(m.m12);
        dst = (*this)[2];
        *dst = T(m.m20); ++dst;
        *dst = T(m.m21); ++dst;
        *dst = T(m.m22);
    }

    /**
     * Sets the value of this matrix to that of the Matrix4d provided. 
     * @param m the source matrix
     */
    void set(const Matrix4d &m) {
        assert(numRow >= 4);
        assert(numCol >= 4);
        T *dst = elem;
        *dst = T(m.m00); ++dst;
        *dst = T(m.m01); ++dst; 
        *dst = T(m.m02); ++dst;
        *dst = T(m.m03);
        dst = (*this)[1];
        *dst = T(m.m10); ++dst;
        *dst = T(m.m11); ++dst;
        *dst = T(m.m12); ++dst;
        *dst = T(m.m13);
        dst = (*this)[2];
        *dst = T(m.m20); ++dst;
        *dst = T(m.m21); ++dst;
        *dst = T(m.m22); ++dst;
        *dst = T(m.m23);
        dst = (*this)[3];
        *dst = T(m.m30); ++dst;
        *dst = T(m.m31); ++dst;
        *dst = T(m.m32); ++dst;
        *dst = T(m.m33);
    }

    /**
     * Sets the value of this matrix to that of the Matrix4f provided. 
     * @param m the source matrix
     */
    void set(const Matrix4f &m) {
        assert(numRow >= 4);
        assert(numCol >= 4);
        T *dst = elem;
        *dst = T(m.m00); ++dst;
        *dst = T(m.m01); ++dst; 
        *dst = T(m.m02); ++dst;
        *dst = T(m.m03);
        dst = (*this)[1];
        *dst = T(m.m10); ++dst;
        *dst = T(m.m11); ++dst;
        *dst = T(m.m12); ++dst;
        *dst = T(m.m13);
        dst = (*this)[2];
        *dst = T(m.m20); ++dst;
        *dst = T(m.m21); ++dst;
        *dst = T(m.m22); ++dst;
        *dst = T(m.m23);
        dst = (*this)[3];
        *dst = T(m.m30); ++dst;
        *dst = T(m.m31); ++dst;
        *dst = T(m.m32); ++dst;
        *dst = T(m.m33);
    }
    
    /**
     * Copy the values from the array into the specified column of this
     * matrix. 
     * @param col the column of this matrix into which the array values
     * will be copied.
     * @param src the source array
     */
    void setColumn(int col, const T *src) {
        assert(numCol > col);
        assert(col >= 0);
        T *dst = elem + col;
        for (int i = 0; i < numRow; ++i, ++src, dst += numCol) {
            *dst = *src;
        }
    }

    /**
     * Copy the values from the array into the specified column of this
     * matrix.
     *
     * @param col the column of this matrix into which the vector values
     * will be copied.
     * @param v the source vector
     */
    inline void setColumn(int col, const GVectorT<T> &v);
    
    /**
     * Modifies the value at the specified row and column of this matrix.
     * @param row the row number to be modified (zero indexed)
     * @param col the column number to be modified (zero indexed)
     * @param value the new matrix element value
     */
    void setElement(int row, int col, T value) {
        assert(numRow > row);
        assert(row >= 0);
        assert(numCol > col);
        assert(col >= 0);
        (*this)[row][col] = value;
    }
    
    /**
     * Sets this GMatrix to the identity matrix.
     */
    void setIdentity(void) {
        setZero();
        const int min = numRow < numCol ? numRow : numCol;
        const int step = numCol + 1;
        T *dst = elem;
        for (int i = 0; i < min; ++i, dst += step) {
            *dst = 1.0;
        }
    }
  
    /**
     * Copy the values from the array into the specified row of this
     * matrix.
     * @param row the row of this matrix into which the array values
     *             will be copied.
     * @param src the source array
     */
    void setRow(int row, T *const src) {
        assert(numRow > row);
        assert(row >= 0);
        T *dst = (*this)[row];
        memcpy(dst, src, sizeof(T) * numCol);
    }
    
    /**
     * Copy the values from the array into the specified row of this
     * matrix. 
     * @param row the row of this matrix into which the vector values
     *             will be copied.
     * @param v the source vector
     */
    inline void setRow(int row, const GVectorT<T> &v);

    /**
     * Sets this matrix to a uniform scale matrix; all of the values are reset.
     * @param scale The new scale value 
     */
    void setScale(T scale) {
        setZero();
        const int min = numRow < numCol ? numRow : numCol;
        const int step = numCol + 1;
        T *dst = elem;
        for (int i = 0; i < min; ++i, dst += step) {
            *dst = scale;
        }
    }

    /**
     * Changes the size of this matrix dynamically.  If the size is increased
     * no data values will be lost.  If the size is decreased, only those data
     * values whose matrix positions were eliminated will be lost.
     *
     * @param row number of desired rows in this matrix
     * @param col number of desired columns in this matrix
     */
    void setSize(int row, int col) {
        assert(row >= 0);
        assert(col >= 0);
        if (numRow == row && numCol == col) {return;}
        const int oldRow = numRow;;
        const int oldCol = numCol;
        resize(row, col);
        int i;
        for (i = 0; i < oldRow; ++i) {
            for (int j = oldCol; j < numCol; ++j) {
                (*this)[i][j] = 0.0;
            }
        }
        for (; i < numRow; ++i) {
            for (int j = 0; j < numCol; ++j) {
                (*this)[i][j] = 0.0;
            }
        }
    }

    /**
     * Sets all the values in this matrix to zero.
     */
    void setZero(void) {
        T *dst = elem;
        for (int i = 0; i < numRowCol; ++i, ++dst) {
            *dst = 0.0;
        }
    }

    /**
     * Sets the value of this matrix to the matrix difference of itself
     * and matrix m (this = this - m).
     * @param m the other matrix
     */
    void sub(const GMatrixT<T> &m) {
        assert(numRow == m.numRow);
        assert(numCol == m.numCol);
        const T *src = m.elem;
        T *dst = elem;
        for (int i = 0; i < numRowCol; ++i, ++dst, ++src) {
            *dst -= *src;
        }
    }

    /**
     * Sets the value of this matrix to the matrix difference
     * of matrices m1 and m2 (this = m1 - m2).
     * @param m1 the first matrix
     * @param m2 the second matrix
     */
    void sub(const GMatrixT<T> &m1, const GMatrixT<T> &m2) {
        assert(numRow == m1.numRow);
        assert(numCol == m1.numCol);
        assert(numRow == m2.numRow);
        assert(numCol == m2.numCol);
        const T *src1 = m1.elem;
        const T *src2 = m2.elem;
        T *dst = elem;
        for (int i = 0; i < numRowCol; ++i, ++dst, ++src1, ++src2) {
            *dst = *src1 - *src2;
        }
    }

protected:
    void setDiag(int i, T value) {
        (*this)[i][i] = value;
    }
    T getDiag(int i) const {
        return (*this)[i][i];
    }

    T dpythag(T a, T b) const {
        const T absa = abs(a);
        const T absb = abs(b);
        if (absa > absb) {
            if (absa == 0.0)
                return 0.0;
            T term = absb/absa;
            if (abs(term) <= DBL_MIN)
                return absa;
            return (absa*sqrt(1.0+term*term));
        } else {
            if (absb == 0.0)
                return 0.0;
            T term = absa/absb;
            if (abs(term) <= DBL_MIN)
                return absb;
            return (absb*sqrt(1.0+term*term));
        }
    }
    
public:    
    /**
     * Finds the singular value decomposition (SVD) of this matrix
     * such that this = U*W*transpose(V); and returns the rank of this
     * matrix; the values of U,W,V are all overwritten. Note that the
     * matrix V is output as V, and not transpose(V). If this matrix
     * is mxn, then U is mxm, W is a diagonal matrix that is mxn, and
     * V is nxn. Using the notation W = diag(w), then the inverse of
     * this matrix is: inverse(this) = V*diag(1/w)*tranpose(U), where
     * diag(1/w) is the same matrix as W except that the reciprocal of
     * each of the diagonal components is used.
     * @param U The computed U matrix in the equation this = U*W*transpose(V) 
     * @param W The computed W matrix in the equation this = U*W*transpose(V) 
     * @param V The computed V matrix in the equation this = U*W*transpose(V) 
     * @return The rank of this matrix. 
     */
    int SVD(GMatrixT<T> &u, GMatrixT<T> &w, GMatrixT<T> &v) const {
        assert(u.numRow == numRow);
        assert(u.numCol == numRow);
        assert(v.numRow == numCol);
        assert(v.numCol == numCol);
        assert(w.numCol == numCol);
        assert(w.numRow == numRow);
        
        const int m = numRow;
        const int n = numCol;
        const int imax = m > n ? m : n;
        GMatrixT<T> &A = u;
        GMatrixT<T> &V = v;
        int i, its, j, jj, k;
        int l = 0, nm = 0;
        T anorm, c, f, g, h, s, scale, x, y, z;
        T *rv1 = new T[n];
        
        // copy this to [u]
        this->get(u);
        // pad 0.0 to the other elements.
        for (i = m; i < imax; i++) {
            for (j = 0; j < imax; j++) {
                A[i][j] = 0.0;
            }
        }
        for (j = n; j < imax; j++) {
            for (i = 0; i < imax; i++) {
                A[i][j] = 0.0;
            }
        }
        // pad 0.0 to w
        w.setZero();

        g = scale = anorm = 0.0;
        for (i = 0; i < n; ++i) {
            l = i + 1;
            rv1[i] = scale * g;
            g = s = scale = 0.0;
            if (i < m) {
                for (k = i; k < m; ++k) {
                    scale += abs(A[k][i]);
                }
                if (scale != 0.0) {
                    for (k = i; k < m; ++k) {
                        A[k][i] /= scale;
                        s += A[k][i] * A[k][i];
                    }
                    f = A[i][i];

                    g = (f < 0.0) ? sqrt(s) : -sqrt(s);
                    h = f * g - s;
                    A[i][i] = f - g;
                    for (j = l; j < n; ++j) {
                        for (s = 0.0, k = i; k < m; ++k) {
                            s += A[k][i] * A[k][j];
                        }
                        f = s / h;
                        for (k = i; k < m; ++k) {
                            A[k][j] += f * A[k][i];
                        }
                    }
                    for (k = i; k < m; ++k) {
                        A[k][i] *= scale;
                    }
                }
            }
            w.setDiag(i, scale * g);
            g = s = scale = 0.0;
            if (i < m && i != n-1) {
                for (k = l; k < n; ++k) {
                    scale += abs(A[i][k]);
                }
                if (scale != 0.0) {
                    for (k = l; k < n; ++k) {
                        A[i][k] /= scale;
                        s += A[i][k] * A[i][k];
                    }
                    f = A[i][l];
                    g = (f < 0.0) ? sqrt(s) : -sqrt(s);
                    h = f * g - s;
                    A[i][l] = f - g;
                    for (k = l; k < n; ++k) {
                        rv1[k] = A[i][k] / h;
                    }
                    for (j = l; j < m; ++j) {
                        for (s = 0.0, k = l; k < n; ++k) {
                            s += A[j][k] * A[i][k];
                        }
                        for (k = l; k < n; ++k) {
                            A[j][k] += s * rv1[k];
                        }
                    }
                    for (k = l; k < n; ++k) {
                        A[i][k] *= scale;
                    }
                }
            }
            T a1 = abs(w.getDiag(i)) + abs(rv1[i]);
            if (a1 > anorm) {
                anorm = a1;
            }
        }
        for (i = n-1; i >= 0; --i) {
            if (i < n-1) {
                if (g != 0.0) {
                    for (j = l; j < n; ++j) {
                        V[j][i]=(A[i][j]/A[i][l])/g;
                    }
                    for (j = l; j < n; ++j) {
                        for (s = 0.0, k = l; k < n; ++k) {
                            s += A[i][k] * V[k][j];
                        }
                        for (k = l; k < n; ++k) {
                            V[k][j] += s * V[k][i];
                        }
                    }
                }
                for (j = l; j < n; ++j) {
                    V[i][j] = V[j][i] = 0.0;
                }
            }
            V[i][i] = 1.0;
            g = rv1[i];
            l = i;
        }

        const int imin = (m < n) ? m : n;
        for (i = imin - 1; i >= 0; --i) {
            l = i + 1;
            g = w.getDiag(i);
            for (j = l; j < n; ++j) {
                A[i][j] = 0.0;
            }
            if (g != 0.0) {
                g = 1.0 / g;
                for (j = l; j < n; ++j) {
                    for (s = 0.0, k = l; k < m; ++k) {
                        s += A[k][i] * A[k][j];
                    }
                    f = (s / A[i][i]) * g;
                    for (k = i; k < m; ++k) {
                        A[k][j] += f * A[k][i];
                    }
                }
                for (j = i; j < m; ++j) {
                    A[j][i] *= g;
                }
            } else {
                for (j = i; j < m; ++j) {
                    A[j][i]=0.0;
                }
            }
            ++A[i][i];
        }
        for (k = n - 1; k >= 0; --k) {
            for (its = 1; its <= 30; ++its) {
                bool flag = true;
                for (l = k; l >= 0; --l) {
                    nm = l - 1;
                    if (T(abs(rv1[l]) + anorm) == anorm) {
                        flag = false;
                        break;
                    }
                    if (T(abs(w.getDiag(nm)) + anorm) == anorm) {
                        break;
                    }
                }
                if (flag) {
                    c = 0.0;
                    s = 1.0;
                    for (i = l; i <= k; ++i) {
                        f = s * rv1[i];
                        rv1[i] = c * rv1[i];
                        if (T(abs(f) + anorm) == anorm) {
                            break;
                        }
                        g = w.getDiag(i);
                        h = dpythag(f, g);
                        w.setDiag(i, h);
                        h = 1.0 / h;
                        c = g * h;
                        s = -f * h;
                        for (j = 0; j < m; ++j) {
                            y = A[j][nm];
                            z = A[j][i];
                            A[j][nm] = y * c + z * s;
                            A[j][i] = z * c - y * s;
                        }
                    }
                }
                z = w.getDiag(k);
                if (l == k) {
                    if (z < 0.0) {
                        w.setDiag(k, -z);
                        for (j = 0; j < n; ++j) {
                            V[j][k] = -V[j][k];
                        }
                    }
                    break;      // normal exit
                }
                if (its == 30) {
                    delete[] rv1;
                    return 0;   // not solved
                }
                x = w.getDiag(l);
                nm = k - 1;
                y = w.getDiag(nm);
                g = rv1[nm];
                h = rv1[k];
                f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                g = dpythag(f, 1.0);

                f = ((x - z) * (x + z)
                     + h * ((y / (f + ((f >= 0.0) ? abs(g) : -abs(g)))) - h))
                    / x;
                c = s = 1.0;
                for (j = l; j <= nm; ++j) {
                    i = j + 1;
                    g = rv1[i];
                    y = w.getDiag(i);
                    h = s * g;
                    g = c * g;
                    z = dpythag(f, h);
                    rv1[j] = z;
                    c = f / z;
                    s = h / z;
                    f = x * c + g * s;
                    g = g * c - x * s;
                    h = y * s;
                    y *= c;
                    for (jj = 0; jj < n; ++jj) {
                        x = V[jj][j];
                        z = V[jj][i];
                        V[jj][j] = x * c + z * s;
                        V[jj][i] = z * c - x * s;
                    }
                    z = dpythag(f, h);
                    w.setDiag(j, z);
                    if (z != 0.0) {
                        z = 1.0 / z;
                        c = f * z;
                        s = h * z;
                    }
                    f = c * g + s * y;
                    x = c * y - s * g;
                    for (jj = 0; jj < m; ++jj) {
                        y = A[jj][j];
                        z = A[jj][i];
                        A[jj][j] = y * c + z * s;
                        A[jj][i] = z * c - y * s;
                    }
                }
                rv1[l] = 0.0;
                rv1[k] = f;
                w.setDiag(k, x);
            }
        }

        // find the number of non-zero w which is the rank of this matrix
        int rank = 0;
        for (i = 0; i < n; i++) {
            if (w.getDiag(i) > 0.0) {
                ++rank;
            }
        }

        delete[] rv1;
        return rank;
    }
    
    /**
     * Returns the trace of this matrix.
     * @return the trace of this matrix.
     */
    T trace(void) const {
        const int min = numRow < numCol ? numRow : numCol;
        const int step = numCol + 1;
        const T *src = elem;
        T t = 0.0;
        for (int i = 0; i < min; ++i, src += step) {
            t += *src;
        }
        return t;
    }

    /**
     * Transposes this matrix in place.
     */
    void transpose(void) {
        assert(numRow == numCol);
        for(int i = 0; i < numRow; ++i) {
            for(int j = i + 1; j < numCol; ++j) {
                T &dst1 = (*this)[i][j];
                T &dst2 = (*this)[j][i];
                T tmp = dst1;
                dst1 = dst2;
                dst2 = tmp;
            }
        }
    }

    /**
     * Places the matrix values of the transpose of matrix m into this matrix.
     * @param m1 the matrix to be transposed (but not modified)
     */
    void transpose(const GMatrixT<T> &m) {
        assert(numRow == numCol);
        T *dst = elem;
        for(int i = 0; i < numRow; ++i) {
            for(int j = 0; j < numCol; ++j, ++dst) {
                *dst = m[j][i];
            }
        }
    }
    
#ifdef VM_INCLUDE_TOSTRING
    VM_STRING_STD::string toString() const;
#endif
};

VM_END_NS

#endif // GMATRIX__H

// Local Variables:
// mode: c++
// End:

