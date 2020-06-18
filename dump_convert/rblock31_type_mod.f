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
      IMPLICIT NONE

c-----------------------------------------------------------------------
c     types for saving data at quadrature points
c-----------------------------------------------------------------------
      TYPE :: rb_real_qp_type
        REAL(r8), DIMENSION(:,:,:,:), POINTER :: qpf,qpfr,qpfz
      END TYPE rb_real_qp_type

      TYPE :: rb_comp_qp_type
        COMPLEX(r8), DIMENSION(:,:,:,:,:), POINTER :: qpf,qpfr,qpfz
      END TYPE rb_comp_qp_type
c-----------------------------------------------------------------------
c     data saved in rb types
c-----------------------------------------------------------------------
      TYPE :: rblock_type
        CHARACTER(64) :: name
        INTEGER(i4) :: id
        INTEGER(i4) :: mx,my
        INTEGER(i4) :: ng
        LOGICAL :: degenerate

        TYPE(bicube_type) :: rz,be_eq,ja_eq,ve_eq,pres_eq,prese_eq,
     $            nd_eq,diff_shape
        TYPE(lagr_quad_2D_type) :: be_n0,pres_n0,nd_n0,neo_mp,neo_p_eq,
     $            e_applied,q_applied,si_nl_pres,neo_p_eqe,tele_eq,
     $            tion_eq,nd_eq2,rwork1
        TYPE(lagr_quad_type) :: be,ja,ve,pres,prese,nd,tele,tion,
     $            conc,neo_p,work1,work2,work3,work4,diff_sh2,neo_pe 

        TYPE(rb_real_qp_type) :: qbe_eq,qja_eq,qve_eq,qpres_eq,
     $            qprese_eq,qnd_eq,qbe_n0,qpres_n0,qnd_n0,qneo_mp,
     $            qneo_p_eq,qe_applied,qq_applied,qsi_nl_pres,
     $            qneo_p_eqe,qbb,qtele_eq,qtion_eq,qrwork1,qbe_tot,
     $            qnd_tot,qelecd_phi,qelecd_n0,qelecd_eq
        TYPE(rb_comp_qp_type) :: qbe,qja,qve,qpres,qprese,qnd,qtele,
     $            qtion,qconc,qneo_p,qdiff_sh2,qneo_pe,qvisc,
     $            qwork1,qwork2,qwork3,qwork4

        REAL(r8), DIMENSION(:), POINTER :: xg,yg,wg
        REAL(r8), DIMENSION(:,:,:), POINTER :: dxdr,dxdz,dydr,dydz,
     $            bigr,wjac
        REAL(r8), DIMENSION(:,:,:,:), POINTER :: alpha,dalpdr,dalpdz
        REAL(r8), DIMENSION(:,:), POINTER :: cell_vol,r_cent
      END TYPE rblock_type
c-----------------------------------------------------------------------
c     close module
c-----------------------------------------------------------------------
      END MODULE rblock_type_mod
