/*
 * Nonfunctional stubs to avoid loading errors when Sequential SuperLU
 * is not available.
 */

#include <stdlib.h>

int
c_fortran_dgssv_(int *iopt, int *n, int *nnz, int *nrhs, int *values,
		 int *rowind, int *colptr, int *b, int *ldb,
		 long int *factors, int *info)

{
    *info=-999999;
}

int
c_fortran_zgssv_(int *iopt, int *n, int *nnz, int *nrhs, int *values,
                 int *rowind, int *colptr, int *b, int *ldb,
                 long int *factors, int *info)
 
{
    *info=-999999;
}
