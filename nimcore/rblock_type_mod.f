c-----------------------------------------------------------------------
c     file rblock_type_mod.f
c     module containing rblock_type definition.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE rblock_type_mod
      USE local
      USE bicube
      USE lagr_quad_mod
      USE lagr_disc_mod
      USE modal_disc_mod
      IMPLICIT NONE

c-----------------------------------------------------------------------
c     types for saving data at quadrature points
c-----------------------------------------------------------------------
      TYPE :: rb_real_qp_type
        REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: qpf,qpfr,qpfz
      END TYPE rb_real_qp_type

      TYPE :: rb_comp_qp_type
        COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: qpf,qpfr,qpfz
      END TYPE rb_comp_qp_type
c-----------------------------------------------------------------------
c     types for holding basis-function values and their derivatives
c     at quadrature points.
c-----------------------------------------------------------------------
      TYPE :: rb_basis_type
        INTEGER(i4) :: poly_deg_basis
        INTEGER(i4) :: poly_degmin_basis
        INTEGER(i4) :: poly_degmax_basis
        REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: alpha,dalpdr,dalpdz,
     $            dalpdrc,alpham,dalpmdr,dalpmdz,dalpmdrc
      END TYPE rb_basis_type
c-----------------------------------------------------------------------
c     data saved in rb types
c-----------------------------------------------------------------------
      TYPE :: rblock_type
        CHARACTER(64) :: name
        INTEGER(i4) :: id
        INTEGER(i4) :: mx,my
        INTEGER(i4) :: ng
        LOGICAL :: degenerate,r0block,periodicy

        TYPE(lagr_quad_2D_type) :: rz,be_eq,ja_eq,ve_eq,pres_eq,
     $            prese_eq,nd_eq,diff_shape,pflux,fllen
        TYPE(lagr_quad_2D_type) :: be_n0,pres_n0,nd_n0,
     $            e_applied,q_applied,si_nl_pres,tele_eq,
     $            tion_eq,nd_eq2,rwork1,rwork2,rwork3,rwork4,ve_n0
        TYPE(lagr_quad_type) :: be,ja,ve,pres,prese,nd,tele,tion,
     $            conc,work1,work2,work3,work4,work5,work6,w6v1,w6v2
        TYPE(modal_quad_type) :: auxv,auxb,mwork1,mwork2,mwork3,mwork4

        TYPE(rb_real_qp_type) :: qbe_eq,qja_eq,qve_eq,qpres_eq,
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
        TYPE(rb_comp_qp_type) :: qbe,qja,qve,qpres,qprese,qnd,qtele,
     $            qtion,qconc,qvisc,qbe_off,qja_off,qve_off,qnd_off,
     $            qwork1,qwork2,qwork3,qwork4,qkapple_off,qkaprpe_off,
     $            qkappli_off,qkaprpi_off,qhyph,qndiff,qndiffa
        TYPE(rb_basis_type), DIMENSION (:), POINTER :: base_pd,
     $            base_disc,base_modal

        REAL(r8), DIMENSION(:), ALLOCATABLE :: xg,yg,wg
        REAL(r8), DIMENSION(:,:), ALLOCATABLE :: dxdr,dxdz,dydr,dydz,
     $            bigr,wjac,jac2d
        REAL(r8), DIMENSION(:,:), ALLOCATABLE :: cell_vol,r_cent
      END TYPE rblock_type
c-----------------------------------------------------------------------
c     close module
c-----------------------------------------------------------------------
      END MODULE rblock_type_mod
