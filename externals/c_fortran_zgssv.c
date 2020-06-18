
/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */

#include "slu_zdefs.h"

#define HANDLE_SIZE  8
/* kind of integer to hold a pointer.  Use int.
   This might need to be changed on 64-bit systems. */
typedef long long fptr;  /* Changed to 64-bit */

typedef struct {
    SuperMatrix *L;
    SuperMatrix *U;
    int *perm_c;
    int *perm_r;
    int *etree;
} factors_t;

void
c_fortran_zgssv_(int *iopt, int *n, int *nnz, int *nrhs, 
                 doublecomplex *values, int *rowind, int *colptr,
                 doublecomplex *b, int *ldb,
		 long int *f_factors, /* a handle containing the address
				     pointing to the factored matrices */
		 int *info)

{
/* 
 * This routine can be called from Fortran.
 *
 * iopt (input) int
 *      Specifies the operation:
 *      = 1, performs LU decomposition for the first time
 *      = 2, performs LU decomposition using the same 
 *           sparsity pattern as the first decomposition
 *      = 3, performs triangular solve
 *      = 4, free all the storage in the end
 *
 * f_factors (input/output) fptr* 
 *      If iopt == 1, it is an output and contains the pointer pointing to
 *                    the structure of the factored matrices.
 *      Otherwise, it it an input.
 *
 */
 
    SuperMatrix A, AC, B;
    SuperMatrix *L, *U;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int *etree;  /* column elimination tree */
    SCformat *Lstore;
    NCformat *Ustore;
    int      i, panel_size, permc_spec, relax, report;
    trans_t  trans;
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    factors_t *LUfactors;

    trans = NOTRANS;

    /* Set the input options. */
    set_default_options(&options);
    options.Trans = trans;
    options.DiagPivotThresh = 0.0;
    options.PrintStat = YES;
    // use minimum degree ordering on structure of A'+A (NIM SuperLU 2.0)
    options.ColPerm = MMD_AT_PLUS_A;

    /*
     * Set option for writing factoring statistics.
     * report = 0: no reporting
     * report = 1: reporting
     */    	
    report = 0;

    if ( *iopt == 1 ) { /* LU decomposition */

        /* Initialize the statistics variables. */
        StatInit(&stat);
        
        zCreate_CompCol_Matrix(&A, *n, *n, *nnz, values, rowind, colptr,
        		       SLU_NR, SLU_Z, SLU_GE);
        L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
        U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
        if ( !(perm_r = intMalloc(*n)) ) ABORT("Malloc fails for perm_r[].");
        if ( !(perm_c = intMalloc(*n)) ) ABORT("Malloc fails for perm_c[].");
        if ( !(etree = intMalloc(*n)) ) ABORT("Malloc fails for etree[].");
        
        /*
         * Get column permutation vector perm_c[], according to permc_spec:
         *   permc_spec = 0: natural ordering 
         *   permc_spec = 1: minimum degree on structure of A'*A
         *   permc_spec = 2: minimum degree on structure of A'+A
         *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
         */    	
        permc_spec = options.ColPerm;
        get_perm_c(permc_spec, &A, perm_c);
        
        sp_preorder(&options, &A, perm_c, etree, &AC);
        
        panel_size = sp_ienv(1);
        relax = sp_ienv(2);
        
        zgstrf(&options, &AC, relax, panel_size, etree,
                NULL, 0, perm_c, perm_r, L, U, &stat, info);
        
        if ( *info == 0 ) {
          if (report == 1) {
            Lstore = (SCformat *) L->Store;
            Ustore = (NCformat *) U->Store;
            printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
            printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
            printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);
            zQuerySpace(L, U, &mem_usage);
            printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
        	   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6); }
        } else {
            printf("zgstrf() error returns INFO= %d\n", *info);
            if ( *info <= *n ) { /* factorization completes */
        	zQuerySpace(L, U, &mem_usage);
        	printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
        	       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
            }
        }
        
        /* Save the LU factors in the factors handle */
        LUfactors = (factors_t*) SUPERLU_MALLOC(sizeof(factors_t));
        LUfactors->L = L;
        LUfactors->U = U;
        LUfactors->perm_c = perm_c;
        LUfactors->perm_r = perm_r;
        LUfactors->etree = etree;
        *f_factors = (fptr) LUfactors;
        
        /* Free un-wanted storage */
        Destroy_SuperMatrix_Store(&A);
        Destroy_CompCol_Permuted(&AC);
        StatFree(&stat);

    } else if ( *iopt == 2 ) { /* Factor a modified matrix with the same
                  sparsity pattern using existing perm_c and L U storage */

        /* Extract the LU factors in the factors handle */
        LUfactors = (factors_t*) *f_factors;
        L = LUfactors->L;
        U = LUfactors->U;
        perm_r = LUfactors->perm_r;
        perm_c = LUfactors->perm_c;
        etree = LUfactors->etree;
      
        /* Set the input option for factorization */
        //options.Fact = SamePattern;
        options.Fact = SamePattern_SameRowPerm;
      
        panel_size = sp_ienv(1);
        relax = sp_ienv(2);
      
        StatInit(&stat);
      
        /* New space is needed for the matrix itself. */
      
        zCreate_CompCol_Matrix(&A, *n, *n, *nnz, values, rowind, colptr,
                   SLU_NR, SLU_Z, SLU_GE);
      
      	sp_preorder(&options, &A, perm_c, etree, &AC);
      
        zgstrf(&options, &AC, relax, panel_size, 
               etree, NULL, 0, perm_r, perm_c, L, U, &stat, info);
      
        if ( *info == 0 ) {
                if (report == 1) {
            Lstore = (SCformat *) L->Store;
            Ustore = (NCformat *) U->Store;
            printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
            printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
            printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);
            //dQuerySpace(L, U, panel_size, &mem_usage);
            zQuerySpace(L, U, &mem_usage);
      	    printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
      		   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6); }
        } else {
            printf("dgstrf() error returns INFO= %d\n", *info);
            if ( *info <= *n ) { /* factorization completes */
      		zQuerySpace(L, U, &mem_usage);
      		printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
      		       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
            }
        }
      
        /* Free un-wanted storage */
        Destroy_SuperMatrix_Store(&A);
        Destroy_CompCol_Permuted(&AC);
              StatFree(&stat);

    } else if ( *iopt == 3 ) { /* Triangular solve */
        /* Initialize the statistics variables. */
        StatInit(&stat);
        
        /* Extract the LU factors in the factors handle */
        LUfactors = (factors_t*) *f_factors;
        L = LUfactors->L;
        U = LUfactors->U;
        perm_c = LUfactors->perm_c;
        perm_r = LUfactors->perm_r;
        
        zCreate_Dense_Matrix(&B, *n, *nrhs, b, *ldb, SLU_DN, SLU_Z, SLU_GE);
        
        /* Solve the system A*X=B, overwriting B with X. */
        zgstrs (trans, L, U, perm_c, perm_r, &B, &stat, info);
        
        Destroy_SuperMatrix_Store(&B);
        StatFree(&stat);

    } else if ( *iopt == 4 ) { /* Free storage */
        /* Free the LU factors in the factors handle */
        LUfactors = (factors_t*) *f_factors;
        SUPERLU_FREE (LUfactors->perm_r);
        SUPERLU_FREE (LUfactors->perm_c);
        SUPERLU_FREE (LUfactors->etree);
        Destroy_SuperNode_Matrix(LUfactors->L);
        Destroy_CompCol_Matrix(LUfactors->U);
        SUPERLU_FREE (LUfactors->L);
        SUPERLU_FREE (LUfactors->U);
        SUPERLU_FREE (LUfactors);
    } else {
        fprintf(stderr,"Invalid iopt=%d passed to c_fortran_zgssv()\n",*iopt);
        exit(-1);
    }
}
