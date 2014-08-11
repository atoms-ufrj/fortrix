#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#if defined(__GNUC__)
#include <stdint.h>
#endif /* __GNUC__ */

#include "blas_arrays.h"

//------------------------------------------------------------------------------
void allocate_blas_array_single(void **array, int n)
{
  *array = malloc(n*sizeof(float));
}
//------------------------------------------------------------------------------
void allocate_blas_array_double(void **array, int n)
{
  *array = malloc(n*sizeof(double));
}
//------------------------------------------------------------------------------
void deallocate_blas_array(void *array)
{
  free(array);
}
//------------------------------------------------------------------------------
void transfer_blas_array_single(void *dest, void *source, int n)
{
  memcpy( dest, source, n*sizeof(float) );
}
//------------------------------------------------------------------------------
void transfer_blas_array_double(void *dest, void *source, int n)
{
  memcpy( dest, source, n*sizeof(double) );
}

