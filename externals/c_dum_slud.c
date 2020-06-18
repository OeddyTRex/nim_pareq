/*
 * Nonfunctional stubs to avoid loading errors when SuperLU_DIST
 * is not available.
 */

#include <stdlib.h>

void
c_fortran_slugrid_(int *iopt, int *slu_comm, int *slu_comm_rank,
                   int *nprow, int *npcol, long int *grid_handle)
{ }

int
c_fortran_pdglob_(int *iopt, int *n, int *nnz, int *nrhs, int *values,
		  int *rowind, int *colptr, int *b, int *ldb,
		  long int *factors, long int *grid_handle,
		  int *info)
{
    *info=-999999;
}

int
c_fortran_pzglob_(int *iopt, int *n, int *nnz, int *nrhs,
		  int *values, int *rowind, int *colptr,
		  int *b, int *ldb,
		  long int *factors, long int *grid_handle,
		  int *info)
{
    *info=-999999;
}

int
c_fortran_pdloc_(int *iopt, int *n, int *mloc, int *nnzloc,
                 int *fstrow, int *nrhs,
                 int *values, int *colind, int *rowptr,
                 int *b, int *ldb,
                 long int *factors, long int *grid_handle,
                 int *info)
{
    *info=-999999;
}

int
c_fortran_pzloc_(int *iopt, int *n, int *mloc, int *nnzloc,
                 int *fstrow, int *nrhs,
                 int *values, int *colind, int *rowptr,
                 int *b, int *ldb,
                 long int *factors, long int *grid_handle,
                 int *info)
{
    *info=-999999;
}
