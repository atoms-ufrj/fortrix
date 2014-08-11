#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

void allocate_blas_array_single(void **array, int n);
void allocate_blas_array_double(void **array, int n);
void deallocate_blas_array(void *array);
void transfer_blas_array_single(void *dest, void *source, int n);
void transfer_blas_array_double(void *dest, void *source, int n);

#if defined(__cplusplus)
}
#endif /* __cplusplus */
