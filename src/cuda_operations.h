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

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

__host__ void cudaSexp( int n, void *a );
__global__ void exp_single( float *a );
__host__ void cudaDexp( int n, void *a );
__global__ void exp_double( double *a );

__host__ void cudaSlog( int n, void *a );
__global__ void log_single( float *a );
__host__ void cudaDlog( int n, void *a );
__global__ void log_double( double *a );

__host__ void cudaSinv( int n, void *a );
__global__ void inv_single( float *a );
__host__ void cudaDinv( int n, void *a );
__global__ void inv_double( double *a );

__host__ void cudaSpow( int n, void *a, void *expnt );
__global__ void pow_single( float *a, float *expnt );
__host__ void cudaDpow( int n, void *a, void *expnt );
__global__ void pow_double( double *a, double *expnt );

__host__ void cudaSprod( int n, void *a, void *b );
__global__ void prod_single( float *a, float *b );
__host__ void cudaDprod( int n, void *a, void *b );
__global__ void prod_double( double *a, double *b );

__host__ void cudaSreplic( int n, void *a, void *b );
__global__ void replic_single( float *a, float *b );
__host__ void cudaDreplic( int n, void *a, void *b );
__global__ void replic_double( double *a, double *b );

#if defined(__cplusplus)
}
#endif /* __cplusplus */

