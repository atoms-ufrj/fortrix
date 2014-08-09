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

#include "cublas_v2.h"

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

cublasStatus_t cublasCreate_f(void **handle);
cublasStatus_t cublasDestroy_f(void *handle);
cublasStatus_t cublasGetVersion_f(void *handle, int *version);
cublasStatus_t cublasSetPointerMode_f(void *handle, int mode);

cublasStatus_t cublasScopy_f(void *handle, int n, void *x, int incx, void *y, int incy);
cublasStatus_t cublasDcopy_f(void *handle, int n, void *x, int incx, void *y, int incy);
cublasStatus_t cublasSaxpy_f(void *handle, int n, void *alpha, void *x, int incx, void *y, int incy);
cublasStatus_t cublasDaxpy_f(void *handle, int n, void *alpha, void *x, int incx, void *y, int incy);
cublasStatus_t cublasSscal_f(void *handle, int n, void *alpha, void *x, int incx);
cublasStatus_t cublasDscal_f(void *handle, int n, void *alpha, void *x, int incx);


cublasStatus_t cublasSgbmv_f(void *handle, int trans, int m, int n, int kl, int ku, void *alpha,
                             void *A, int lda, void *x, int incx, void *beta, void *y, int incy);
cublasStatus_t cublasDgbmv_f(void *handle, int trans, int m, int n, int kl, int ku, void *alpha,
                             void *A, int lda, void *x, int incx, void *beta, void *y, int incy);
cublasStatus_t cublasSgemv_f(void *handle, int trans, int m, int n, void *alpha,
                             void *A, int lda, void *x, int incx, void *beta, void *y, int incy);
cublasStatus_t cublasDgemv_f(void *handle, int trans, int m, int n, void *alpha,
                             void *A, int lda, void *x, int incx, void *beta, void *y, int incy);
cublasStatus_t cublasSger_f(void *handle, int m, int n, void *alpha, void *x, int incx,
                            void *y, int incy, void *A, int lda);
cublasStatus_t cublasDger_f(void *handle, int m, int n, void *alpha, void *x, int incx,
                            void *y, int incy, void *A, int lda);


cublasStatus_t cublasSgemm_f(void *handle, int transa, int transb, int m, int n, int k,
                             void *alpha, void *A, int lda, void *B, int ldb,
                             void *beta, void *C, int ldc);
cublasStatus_t cublasDgemm_f(void *handle, int transa, int transb, int m, int n, int k,
                             void *alpha, void *A, int lda, void *B, int ldb,
                             void *beta, void *C, int ldc);


cublasStatus_t cublasSgeam_f(void *handle, int transa, int transb, int m, int n, void *alpha,
                             void *A, int lda, void *beta, void *B, int ldb, void *C, int ldc);
cublasStatus_t cublasDgeam_f(void *handle, int transa, int transb, int m, int n, void *alpha,
                             void *A, int lda, void *beta, void *B, int ldb, void *C, int ldc);
cublasStatus_t cublasSdgmm_f(void *handle, int mode, int m, int n, void *A, int lda,
                             void *x, int incx, void *C, int ldc);
cublasStatus_t cublasDdgmm_f(void *handle, int mode, int m, int n, void *A, int lda,
                             void *x, int incx, void *C, int ldc);


#if defined(__cplusplus)
}
#endif /* __cplusplus */

