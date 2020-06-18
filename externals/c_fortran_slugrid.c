#include <stdlib.h>
#include <stdio.h>

#include "config.h"
#include "int_size.h"
#include "superlu_ddefs.h"
#include "util_dist.h"
#include "Cnames.h"
#include "Cnames2.h"

//void c_fortran_slugrid_(int *iopt, MPI_Comm *slu_comm, int *slu_comm_rank,
void c_fortran_slugrid_(int *iopt, MPI_Fint *f_slu_comm, int *slu_comm_rank,
                   int *nprow, int *npcol, long int *grid_handle)
/*
 * This routine provides a fortran call for initializing and 
 * freeing the SuperLU_DIST processor grid.  The pointer for the grid
 * structure is returned in grid_handle.
 *
 * The input option, iopt, controls the functionality:
 *   iopt=1:  allocate and define a new process grid
 *   iopt=2:  free an existing process grid
 *
 * slu_comm is the base communication handle
 * slu_comm_rank is this processor's rank in slu_comm
 * nprow is the number of processors per process grid row
 * npcol is the number of processors per process grid column
 */

{
    MPI_Comm slu_comm;
    gridinfo_t *grid;
    int_t nprow_t=*nprow;
    int_t npcol_t=*npcol;

    slu_comm = MPI_Comm_f2c( *f_slu_comm);

    if ( *iopt == 1 ) {

    /* Allocate the grid structure. */

    grid = (gridinfo_t *) SUPERLU_MALLOC(sizeof(gridinfo_t));

    /* Initialize the process grid. */
   

    superlu_gridinit(slu_comm, nprow_t, npcol_t, grid);

    /* Set the handle passed from fortran, so that the
     * process grid can be reused. */

    *grid_handle = (long) grid;
    }

    else if ( *iopt == 2 ) {

    /* Locate and free the process grid. */

    grid = (gridinfo_t *) *grid_handle;
    superlu_gridexit(grid);
    }
}
