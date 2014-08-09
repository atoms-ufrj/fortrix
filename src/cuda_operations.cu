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

#include "cuda_operations.h"
#include "math.h"

//------------------------------------------------------------------------------

__host__ void cudaSexp( int n, void *a )
{
  exp_single<<<1,n>>>( (float*)a );
}

__global__ void exp_single( float *a )
{
  int i = threadIdx.x;
  a[i] = expf(a[i]);
}

//------------------------------------------------------------------------------

__host__ void cudaDexp( int n, void *a )
{
  exp_double<<<1,n>>>( (double*)a );
}

__global__ void exp_double( double *a )
{
  int i = threadIdx.x;
  a[i] = exp(a[i]);
}

//------------------------------------------------------------------------------

__host__ void cudaSlog( int n, void *a )
{
  log_single<<<1,n>>>( (float*)a );
}

__global__ void log_single( float *a )
{
  int i = threadIdx.x;
  a[i] = logf(a[i]);
}

//------------------------------------------------------------------------------

__host__ void cudaDlog( int n, void *a )
{
  log_double<<<1,n>>>( (double*)a );
}

__global__ void log_double( double *a )
{
  int i = threadIdx.x;
  a[i] = log(a[i]);
}

//------------------------------------------------------------------------------

__host__ void cudaSinv( int n, void *a )
{
  inv_single<<<1,n>>>( (float*)a );
}

__global__ void inv_single( float *a )
{
  int i = threadIdx.x;
  a[i] = 1.0f/a[i];
}

//------------------------------------------------------------------------------

__host__ void cudaDinv( int n, void *a )
{
  inv_double<<<1,n>>>( (double*)a );
}

__global__ void inv_double( double *a )
{
  int i = threadIdx.x;
  a[i] = 1.0/a[i];
}

//------------------------------------------------------------------------------

__host__ void cudaSpow( int n, void *a, void *expnt )
{
  pow_single<<<1,n>>>( (float*)a, (float*)expnt );
}

__global__ void pow_single( float *a, float *expnt )
{
  int i = threadIdx.x;
  a[i] = powf(a[i],*expnt);
}

//------------------------------------------------------------------------------

__host__ void cudaDpow( int n, void *a, void *expnt )
{
  pow_double<<<1,n>>>( (double*)a, (double*)expnt );
}

__global__ void pow_double( double *a, double *expnt )
{
  int i = threadIdx.x;
  a[i] = pow(a[i],*expnt);
}

//------------------------------------------------------------------------------

__host__ void cudaSprod( int n, void *a, void *b )
{
  prod_single<<<1,n>>>( (float*)a, (float*)b );
}

__global__ void prod_single( float *a, float *b )
{
  int i = threadIdx.x;
  b[i] = a[i]*b[i];
}

//------------------------------------------------------------------------------

__host__ void cudaDprod( int n, void *a, void *b )
{
  prod_double<<<1,n>>>( (double*)a, (double*)b );
}

__global__ void prod_double( double *a, double *b )
{
  int i = threadIdx.x;
  b[i] = a[i]*b[i];
}

//------------------------------------------------------------------------------

__host__ void cudaSreplic( int n, void *a, void *b )
{
  replic_single<<<1,n>>>( (float*)a, (float*)b );
}

__global__ void replic_single( float *a, float *b )
{
  int i = threadIdx.x;
  b[i] = a[0];
}

//------------------------------------------------------------------------------

__host__ void cudaDreplic( int n, void *a, void *b )
{
  replic_double<<<1,n>>>( (double*)a, (double*)b );
}

__global__ void replic_double( double *a, double *b )
{
  int i = threadIdx.x;
  b[i] = a[0];
}

