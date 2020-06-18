/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */

/*
 * This is a modified version of c_fortran_dgssv
 * for calling the SuperLU_DIST routine pdgssvx with distributed.
 * matrix storage.
 */

#include <stdlib.h>
#include <stdio.h>

#include "config.h"
#include "int_size.h"
#ifdef HAVE_SUPERLU_DIST5
#include "superlu_defs.h"
#define superlu_options_t superlu_dist_options_t
#endif
#ifdef HAVE_SUPERLU_DIST6
#include "superlu_defs.h"
#define superlu_options_t superlu_dist_options_t
#endif
#include "superlu_ddefs.h"
#include "Cnames2.h"
#include "enum_defs.h"

typedef struct {
    ScalePermstruct_t *ScalePermstruct;
    LUstruct_t *LUstruct;
    SOLVEstruct_t *SOLVEstruct;
    int_t solveinit_flag;
#ifdef OBJ_MEM_PROF
    int lu_mem_id;
    int tot_mem_id;
#endif
} factors_dist_t;

#ifdef OBJ_MEM_PROF
    extern void* memlog_update(int *mem_id, char *obj_type,
                   char *obj_name, int *mem_size, int *c_resize);
#endif

int_t
c_fortran_pdloc_(int *iopt, int_t *n, int_t *mloc, int_t *nnzloc,
                 int_t *fstrow, int *nrhs,
                 double *values, int_t *colind, int_t *rowptr,
                 double *b, int *ldb,
                 long int *factors, long int *grid_handle,
                 int* iparm, double* dparm, int *info)

{
/* 
 * This routine can be called from Fortran.
 *
 * iopt (input) int
 *      Specifies the operation:
 *      = 1, performs LU decomposition for the first time
 *      = 2, performs a subsequent LU decomposition for a new matrix
 *           with the same sparsity pattern
 *      = 3, performs triangular solve
 *      = 4, free all the storage in the end
 *
 * factors (input/output) is a long int
 *      If iopt == 1, it is an output and contains the pointer pointing to
 *                    the structure of the factored matrices.
 *      Otherwise, it it an input.
 *
 * grid_handle is a pointer to the process grid structure, which is
 *	created and freed separately.
 */
 
    superlu_options_t options;
    SuperLUStat_t stat;
    SuperMatrix A;
    ScalePermstruct_t *ScalePermstruct;
    LUstruct_t *LUstruct;
    SOLVEstruct_t *SOLVEstruct;
    double   *berr;
    int_t    nprow, npcol;
    int_t    iam;
    int_t    report;
    int_t    i;
    gridinfo_t *grid;
    factors_dist_t *LUfactors;

/* TMP definitions for temp writes */
   int_t nglob;
   int_t mglob;
   int_t fst_row;
   int_t nnz_loc;
   int_t m_loc;
 
   int_t *rowptr2;
   NRformat_loc *Aptr;
   double *nzval;
 
/* end of definitions for temp writes */

    /*
     * Set option for writing factoring statistics.
     * report = 0: no reporting
     * report = 1: reporting
     */    	
    report = iparm[0];

    /* Locate the process grid. */

    grid = (gridinfo_t *) *grid_handle;
    iam = (*grid).iam;
    nprow = (int_t) (*grid).nprow;
    npcol = (int_t) (*grid).npcol;

    /* Set options. */

    set_default_options_dist(&options);
    options.Equil = NO;
    if (iparm[1]) options.Equil = YES;
    switch (iparm[2]) {
#ifdef HAVE_SUPERLU_DIST6
      case 1: options.RowPerm = LargeDiag_MC64;
      case 2: options.RowPerm = LargeDiag_AWPM;
#else
      case 1: options.RowPerm = LargeDiag;
#endif
      default: options.RowPerm = NO;
    }
    switch (iparm[3]) {
      case 0: options.ColPerm = NATURAL;
      case 1: options.ColPerm = MMD_ATA;
      case 2: options.ColPerm = MMD_AT_PLUS_A;
      case 3: options.ColPerm = COLAMD;
#ifdef HAVE_PARMETIS
      case 4: options.ColPerm = METIS_AT_PLUS_A;
#endif
      default: options.ColPerm = MMD_AT_PLUS_A;
    }
    options.ReplaceTinyPivot = NO;
    options.IterRefine = NO;
    options.Trans = NOTRANS;
    options.SolveInitialized = NO;
    options.RefineInitialized = NO;
    options.PrintStat = NO;
    if (iparm[4]) options.PrintStat = YES;

    if ( *iopt == 1 ) { /* LU decomposition */

        /* if ( iam >= nprow * npcol ) return; */

	PStatInit(&stat);


	dCreate_CompRowLoc_Matrix_dist(&A, *n, *n, *nnzloc, *mloc, *fstrow,
                                       values, colind, rowptr,
			               SLU_NR_loc, SLU_D, SLU_GE);

        /* Set options. */

        options.Fact = DOFACT;

	/* Initialize ScalePermstruct and LUstruct. */

        mglob = (int_t) A.nrow;
        nglob = (int_t) A.ncol;

        ScalePermstruct =
            (ScalePermstruct_t *) SUPERLU_MALLOC(sizeof(ScalePermstruct_t));
        ScalePermstructInit(mglob, nglob, ScalePermstruct);

        LUstruct = (LUstruct_t *) SUPERLU_MALLOC(sizeof(LUstruct_t));
#ifdef HAVE_SUPERLU_DIST6
        LUstructInit(nglob, LUstruct);
#else
        LUstructInit(mglob, nglob, LUstruct);
#endif

        SOLVEstruct = (SOLVEstruct_t *) SUPERLU_MALLOC(sizeof(SOLVEstruct_t));

	/* Call the routine with nrhs=0 to perform the factorization. */

	pdgssvx(&options, &A, ScalePermstruct, NULL, *ldb, 0, 
	        grid, LUstruct, SOLVEstruct, berr, &stat, info);

	if ( *info == 0 ) {
          if ( report == 1 ) PStatPrint(&options, &stat, grid);
	} else {
	    printf("pdgssvx() error returns INFO= %d\n", *info);
	}
	
	/* Save the LU factors in the factors handle */

	LUfactors = (factors_dist_t*) SUPERLU_MALLOC(sizeof(factors_dist_t));
	LUfactors->ScalePermstruct = ScalePermstruct;
	LUfactors->LUstruct = LUstruct;
	LUfactors->SOLVEstruct = SOLVEstruct;
        LUfactors->solveinit_flag = 0;
	*factors = (long) LUfactors;

	/* Free un-wanted storage */

	Destroy_SuperMatrix_Store_dist(&A);
        PStatFree(&stat);

    } else if ( *iopt == 2 ) { /* Factor a modified matrix with the same
	sparsity pattern using existing permutations and L U storage */

	/* Extract the LU factors in the factors handle */

	LUfactors = (factors_dist_t*) *factors;
	ScalePermstruct = LUfactors->ScalePermstruct;
	LUstruct = LUfactors->LUstruct;
	SOLVEstruct = LUfactors->SOLVEstruct;

	PStatInit(&stat);

	/* Reset SuperMatrix pointers. */

	dCreate_CompRowLoc_Matrix_dist(&A, *n, *n, *nnzloc, *mloc, *fstrow,
                                       values, colind, rowptr,
			               SLU_NR_loc, SLU_D, SLU_GE);

	/* Set options. */

        options.Fact = SamePattern_SameRowPerm;

	/* Call the routine with nrhs=0 to perform the factorization. */

	pdgssvx(&options, &A, ScalePermstruct, NULL, *ldb, 0, 
	        grid, LUstruct, SOLVEstruct, berr, &stat, info);

	if ( *info == 0 ) {
          if ( report == 1 ) PStatPrint(&options, &stat, grid);
	} else {
	    printf("pdgssvx() error returns INFO= %d\n", *info);
	}
	
	/* Free un-wanted storage */

	Destroy_SuperMatrix_Store_dist(&A);
        PStatFree(&stat);


    } else if ( *iopt == 3 ) { /* Triangular solve */

	/* Extract the LU factors in the factors handle */

	LUfactors = (factors_dist_t*) *factors;
	ScalePermstruct = LUfactors->ScalePermstruct;
	LUstruct = LUfactors->LUstruct;
	SOLVEstruct = LUfactors->SOLVEstruct;

	PStatInit(&stat);

	/* Reset SuperMatrix pointers. */

	dCreate_CompRowLoc_Matrix_dist(&A, *n, *n, *nnzloc, *mloc, *fstrow,
                                       values, colind, rowptr,
			               SLU_NR_loc, SLU_D, SLU_GE);

        /* Allocate error array. */

        if ( !(berr = doubleMalloc_dist(*nrhs)) )
          ABORT("doubleMalloc_dist fails for berr[].");

	/* Set options. */

        options.Fact = FACTORED;
        if (LUfactors->solveinit_flag == 0) {
          options.SolveInitialized = NO;
          LUfactors->solveinit_flag = 1; 
        } else {
          options.SolveInitialized = YES;
        }

        /* Solve the system A*X=B, overwriting B with X. */

	pdgssvx(&options, &A, ScalePermstruct, b, *ldb, *nrhs, 
	        grid, LUstruct, SOLVEstruct, berr, &stat, info);

	Destroy_SuperMatrix_Store_dist(&A);
        PStatFree(&stat);
        SUPERLU_FREE(berr);

    } else if ( *iopt == 4 ) { /* Free storage */

	/* Free the LU factors in the factors handle */

	LUfactors = (factors_dist_t*) *factors;
	Destroy_LU(*n, grid, LUfactors->LUstruct);
        LUstructFree(LUfactors->LUstruct);
        options.SolveInitialized = YES;
        options.RefineInitialized = NO;
        dSolveFinalize(&options, LUfactors->SOLVEstruct);
	ScalePermstructFree(LUfactors->ScalePermstruct);
        SUPERLU_FREE(LUfactors);

    } else {

	fprintf(stderr, "Invalid iopt=%d passed to c_fortran_pdloc()\n",*iopt);
	exit(-1);

    }

#ifdef OBJ_MEM_PROF
#ifdef HAVE_SUPERLU_DIST6
    superlu_dist_mem_usage_t mem_usage;
#else
    mem_usage_t mem_usage;
#endif
    int for_lu,total,creport=1;
    if ( *iopt == 1 ) {
      LUfactors->lu_mem_id=0;
      LUfactors->tot_mem_id=0;
    }
    if ( *iopt != 4 ) {
#ifdef HAVE_SUPERLU_DIST6
      dQuerySpace_dist(*n,LUfactors->LUstruct,grid,&stat,&mem_usage);
#else
      dQuerySpace_dist(*n,LUfactors->LUstruct,grid,&mem_usage);
#endif
      for_lu=(int)mem_usage.for_lu;
      total=(int)mem_usage.total;
    } else {
      for_lu=0; total=0;
    }
    memlog_update(&LUfactors->lu_mem_id,"SLUd_LU","unknown",&for_lu,&creport);
    memlog_update(&LUfactors->tot_mem_id,"SLUd_total","unknown",&total,
                  &creport);
#endif
    return 0;
}


