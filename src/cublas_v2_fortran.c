/*
    This file is part of Fortrix (version 0.1)

    Fortrix is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Fortrix is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.!

    You should have received a copy of the GNU General Public License
    along with Fortrix.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#if defined(__GNUC__)
#include <stdint.h>
#endif /* __GNUC__ */
#include "cublas_v2_fortran.h"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// CUBLAS HELPER FUNCTIONS:
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cublasStatus_t cublasCreate_f(void **handle)
{
  cublasHandle_t *handle_ptr;
  handle_ptr = malloc(sizeof(cublasHandle_t));
  cublasStatus_t stat = cublasCreate(handle_ptr);
  *handle = handle_ptr;
  return stat;
}
//--------------------------------------------------------------------------------------------------
cublasStatus_t cublasDestroy_f(void *handle)
{
  return cublasDestroy(*(cublasHandle_t*)handle);
}
//--------------------------------------------------------------------------------------------------
cublasStatus_t cublasGetVersion_f(void *handle, int *version)
{
  return cublasGetVersion(*(cublasHandle_t*)handle, version);
}
//--------------------------------------------------------------------------------------------------
cublasStatus_t cublasSetPointerMode_f(void *handle, int mode)
{
  return cublasSetPointerMode(*(cublasHandle_t*)handle,(cublasPointerMode_t)mode);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// CUBLAS LEVEL-1 FUNCTIONS:
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cublasStatus_t cublasScopy_f(void *handle, int n, void *x, int incx, void *y, int incy)
{
  return cublasScopy(*(cublasHandle_t*)handle, n, (float*)x, incx, (float*)y, incy);
}
//--------------------------------------------------------------------------------------------------
cublasStatus_t cublasDcopy_f(void *handle, int n, void *x, int incx, void *y, int incy)
{
  return cublasDcopy(*(cublasHandle_t*)handle, n, (double*)x, incx, (double*)y, incy);
}
//--------------------------------------------------------------------------------------------------
cublasStatus_t cublasSaxpy_f(void *handle, int n, void *alpha, void *x, int incx, void *y, int incy)
{
  return cublasSaxpy(*(cublasHandle_t*)handle, n, (float*)alpha, (float*)x, incx, (float*)y, incy);
}
//--------------------------------------------------------------------------------------------------
cublasStatus_t cublasDaxpy_f(void *handle, int n, void *alpha, void *x, int incx, void *y, int incy)
{
  return cublasDaxpy(*(cublasHandle_t*)handle, n, (double*)alpha,
                     (double*)x, incx, (double*)y, incy);
}
//--------------------------------------------------------------------------------------------------
cublasStatus_t cublasSscal_f(void *handle, int n, void *alpha, void *x, int incx)
{
  return cublasSscal(*(cublasHandle_t*)handle, n, (float*)alpha, (float*)x, incx);
}
//--------------------------------------------------------------------------------------------------
cublasStatus_t cublasDscal_f(void *handle, int n, void *alpha, void *x, int incx)
{
  return cublasDscal(*(cublasHandle_t*)handle, n, (double*)alpha, (double*)x, incx);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// CUBLAS LEVEL-2 FUNCTIONS:
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cublasStatus_t cublasSgbmv_f(void *handle, int trans, int m, int n, int kl, int ku, void *alpha,
                             void *A, int lda, void *x, int incx, void *beta, void *y, int incy)
{
  return cublasSgbmv(*(cublasHandle_t*)handle, (cublasOperation_t)trans, m, n, kl, ku,
                   (float*)alpha, (float*)A, lda, (float*)x, incx, (float*)beta, (float*)y, incy);
}
//--------------------------------------------------------------------------------------------------
cublasStatus_t cublasDgbmv_f(void *handle, int trans, int m, int n, int kl, int ku, void *alpha,
                             void *A, int lda, void *x, int incx, void *beta, void *y, int incy)
{
  return cublasDgbmv(*(cublasHandle_t*)handle, (cublasOperation_t)trans, m, n, kl, ku,
               (double*)alpha, (double*)A, lda, (double*)x, incx, (double*)beta, (double*)y, incy);
}
//--------------------------------------------------------------------------------------------------
cublasStatus_t cublasSgemv_f( void *handle, int trans, int m, int n, void *alpha,
                              void *A, int lda, void *x, int incx, void *beta, void *y, int incy)
{
  return cublasSgemv(*(cublasHandle_t*)handle, (cublasOperation_t)trans, m, n, (float*)alpha,
                      (float*)A, lda, (float*)x, incx, (float*)beta, (float*)y, incy);
}
//--------------------------------------------------------------------------------------------------
cublasStatus_t cublasDgemv_f( void *handle, int trans, int m, int n, void *alpha,
                              void *A, int lda, void *x, int incx, void *beta, void *y, int incy)
{
  return cublasDgemv(*(cublasHandle_t*)handle, (cublasOperation_t)trans, m, n, (double*)alpha,
                      (double*)A, lda, (double*)x, incx, (double*)beta, (double*)y, incy);
}
//--------------------------------------------------------------------------------------------------
cublasStatus_t cublasSger_f(void *handle, int m, int n, void *alpha, void *x, int incx,
                            void *y, int incy, void *A, int lda)
{
  return cublasSger(*(cublasHandle_t*)handle, m, n, (float*)alpha, (float*)x, incx,
                     (float*)y, incy, (float*)A, lda);
}
//--------------------------------------------------------------------------------------------------
cublasStatus_t cublasDger_f(void *handle, int m, int n, void *alpha, void *x, int incx,
                            void *y, int incy, void *A, int lda)
{
  return cublasDger(*(cublasHandle_t*)handle, m, n, (double*)alpha, (double*)x, incx,
                     (double*)y, incy, (double*)A, lda);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// CUBLAS LEVEL-3 FUNCTIONS:
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cublasStatus_t cublasSgemm_f(void *handle, int transa, int transb, int m, int n, int k,
                             void *alpha, void *A, int lda, void *B, int ldb,
                             void *beta, void *C, int ldc)
{
  return cublasSgemm(*(cublasHandle_t*)handle, (cublasOperation_t)transa,
                      (cublasOperation_t)transb, m, n, k, (float*)alpha, (float*)A, lda,
                      (float*)B, ldb, (float*)beta, (float*)C, ldc);
}
//--------------------------------------------------------------------------------------------------
cublasStatus_t cublasDgemm_f(void *handle, int transa, int transb, int m, int n, int k,
                             void *alpha, void *A, int lda, void *B, int ldb,
                             void *beta, void *C, int ldc)
{
  return cublasDgemm(*(cublasHandle_t*)handle, (cublasOperation_t)transa,
                      (cublasOperation_t)transb, m, n, k, (double*)alpha, (double*)A, lda,
                      (double*)B, ldb, (double*)beta, (double*)C, ldc);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// BLAS-LIKE EXTENSIONS:
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cublasStatus_t cublasSgeam_f(void *handle, int transa, int transb, int m, int n, void *alpha,
                             void *A, int lda, void *beta, void *B, int ldb, void *C, int ldc)
{
  return cublasSgeam(*(cublasHandle_t*)handle, (cublasOperation_t)transa,
                     (cublasOperation_t)transb, m, n, (float*)alpha, (float*)A,
                     lda, (float*)beta, (float*)B, ldb, (float*)C, ldc);
}
//--------------------------------------------------------------------------------------------------
cublasStatus_t cublasDgeam_f(void *handle, int transa, int transb, int m, int n, void *alpha,
                             void *A, int lda, void *beta, void *B, int ldb, void *C, int ldc)
{
  return cublasDgeam(*(cublasHandle_t*)handle, (cublasOperation_t)transa,
                     (cublasOperation_t)transb, m, n, (double*)alpha, (double*)A,
                     lda, (double*)beta, (double*)B, ldb, (double*)C, ldc);
}
//--------------------------------------------------------------------------------------------------
cublasStatus_t cublasSdgmm_f(void *handle, int mode, int m, int n, void *A, int lda,
                             void *x, int incx, void *C, int ldc)
{
  return cublasSdgmm(*(cublasHandle_t*)handle, (cublasSideMode_t)mode, m, n, (float*)A,
                      lda, (float*)x, incx, (float*)C, ldc);
}
//--------------------------------------------------------------------------------------------------
cublasStatus_t cublasDdgmm_f(void *handle, int mode, int m, int n, void *A, int lda,
                             void *x, int incx, void *C, int ldc)
{
  return cublasDdgmm(*(cublasHandle_t*)handle, (cublasSideMode_t)mode, m, n, (double*)A,
                      lda, (double*)x, incx, (double*)C, ldc);
}
//--------------------------------------------------------------------------------------------------

