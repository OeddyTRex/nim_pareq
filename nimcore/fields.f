c-----------------------------------------------------------------------
c     file fields.f
c     module containing block and block-connection structures; all grid
c     and spatially dependent variables are within.
c     note:  nbl_total = total # of blocks in problem (on all
c     processors), and nbl = # of blocks stored on this processor.
c     However, for nimset, only nbl is used--nimset does not run in
c     parallel.
c     nrbl (nrbl_total) is the number of rectangular blocks numbered 
c     from 1 to nrbl, and the block of unstructured triangles are
c     indexed from nrbl+1 to nbl (nrbl_total+1 to nbl_total).
c-----------------------------------------------------------------------
      MODULE fields
      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE vector_type_mod
      IMPLICIT NONE

      INTEGER(i4) :: nbl_total,nrbl_total
      INTEGER(i4) :: nbl,nrbl
      TYPE(rblock_type), DIMENSION(:), POINTER :: rb
      TYPE(tblock_type), DIMENSION(:), POINTER :: tb
c-----------------------------------------------------------------------
c     the following are pointers to nodal data.  the form is the same
c     for all blocks.
c-----------------------------------------------------------------------
      TYPE(cvector_type), DIMENSION(:), POINTER ::
     $  be,ja,ve,pres,prese,nd,tele,tion,conc,work1,work2,work3,work4,
     $  work5,work6,auxb,auxv,w6v1,w6v2
      TYPE(cvector_type), DIMENSION(:), POINTER :: be_old,ja_old,ve_old,
     $  p_old,pe_old,nd_old,te_old,ti_old
      TYPE(vector_type), DIMENSION(:), POINTER ::
     $  be_eq,ja_eq,ve_eq,pres_eq,prese_eq,nd_eq,diff_shape,
     $  be_n0,pres_n0,nd_n0,rwork1,rwork2,rwork3,rwork4,e_applied,
     $  q_applied,si_nl_pres,tele_eq,tion_eq,nd_eq2,ve_n0,elecd_n0,
     $  kappli_n0,kapple_n0,kaprpi_n0,kaprpe_n0,pflux,fllen
c-----------------------------------------------------------------------
c     close module
c-----------------------------------------------------------------------
      END MODULE fields
