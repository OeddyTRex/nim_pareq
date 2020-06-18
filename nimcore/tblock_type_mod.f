c-----------------------------------------------------------------------
c     file tblock_type_mod.f
c     contains a module that defines tblock types.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE tblock_type_mod
      USE local
      USE tri_linear
      IMPLICIT NONE

c-----------------------------------------------------------------------
c     types for saving data at quadrature points
c-----------------------------------------------------------------------
      TYPE :: tb_real_qp_type
        REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: qpf,qpfr,qpfz
      END TYPE tb_real_qp_type

      TYPE :: tb_comp_qp_type
        COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: qpf,qpfr,qpfz
      END TYPE tb_comp_qp_type
c-----------------------------------------------------------------------
c     data saved in tblock types
c-----------------------------------------------------------------------
      TYPE :: tblock_type
        CHARACTER(64) :: name
        INTEGER(i4) :: id
        INTEGER(i4) :: mvert,mcell
        INTEGER(i4) :: ng
        TYPE(tri_linear_geom_type) :: tgeom

        TYPE(tri_linear_2D_type) :: be_eq,ja_eq,ve_eq,pres_eq,prese_eq,
     $            nd_eq,diff_shape,pflux,fllen,be_n0,pres_n0,nd_n0,
     $            e_applied,q_applied,si_nl_pres,tele_eq,
     $            tion_eq,nd_eq2,rwork1,rwork2,rwork3,rwork4,ve_n0
        TYPE(tri_linear_type) :: be,ja,ve,pres,prese,nd,tele,tion,
     $            conc,work1,work2,work3,work4,work5,work6,w6v1,w6v2
        TYPE(tri_linear_type) :: auxv,auxb,mwork1,mwork2,mwork3,mwork4

        TYPE(tb_real_qp_type) :: qbe_eq,qja_eq,qve_eq,qpres_eq,
     $            qprese_eq,qnd_eq,qbe_n0,qpres_n0,qnd_n0,
     $            qe_applied,qq_applied,qsi_nl_pres,
     $            qbb,qtele_eq,qtion_eq,qrwork1,qrwork4,qbe_tot,
     $            qnd_tot,qelecd_phi,qelecd_n0,qelecd_eq,qdiff_shape,
     $            qkappli_phi,qkapple_phi,qkappli_n0,qkapple_n0,
     $            qkaprpi_phi,qkaprpe_phi,qkaprpi_n0,qkaprpe_n0,
     $            qte_b2,qti_b2,qbcrgte,qbcrgti,
     $            qve_n0,qve_tot,qgrdv,qja_tot,qdvv_eq,qeq_force,
     $            qdart,qupw_phi,qupw_n0,qvv,qupti_phi,qupti_n0,
     $            qupte_phi,qupte_n0,qti_n0,qte_n0,qgrdveq,qpi_veq,
     $            qpi_pareq,qpi_gyreq,qgrdb,qpr_tot,qgrdp,qti_tot,
     $            qndiff_phi,qndiff_n0
        TYPE(tb_comp_qp_type) :: qbe,qja,qve,qpres,qprese,qnd,qtele,
     $            qtion,qconc,qvisc,qbe_off,qja_off,qve_off,qnd_off,
     $            qwork1,qwork2,qwork3,qwork4,qkapple_off,qkaprpe_off,
     $            qkappli_off,qkaprpi_off,qhyph,qndiff,qndiffa

        REAL(r8), DIMENSION(:), ALLOCATABLE :: wg
        REAL(r8), DIMENSION(:), ALLOCATABLE :: cell_vol,r_cent
      END TYPE tblock_type
c-----------------------------------------------------------------------
c     type used for defining the mass matrices.
c-----------------------------------------------------------------------
      TYPE :: matrix_element_type1
        REAL(r8), DIMENSION(:), ALLOCATABLE :: element
      END TYPE matrix_element_type1
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE tblock_type_mod
