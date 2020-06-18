c-----------------------------------------------------------------------
c     file global.f
c     contains run-time information that is not directly associated
c     with input.
c-----------------------------------------------------------------------
      MODULE global
      USE local
      IMPLICIT NONE

      REAL(r8) :: t,dt=0,dt_old=0,delta2,kdivb_2,smallnum
      REAL(r8) :: fl_cfl=0,lin_cfl=0,nl_cfl=0
      REAL(r8) :: cross_section,total_volume,cross_s_overr
      REAL(r8) :: i_eq,i_n0,flux_eq,flux_n0,ff,theta
      REAL(r8) :: volt=0,volt_old=0,i_n0_old
      REAL(r8) :: e_vert=0,e_vert_old=0
      REAL(r8) :: coefjvi,coefjve,coefhll,coefgpe,coefme1,coefme2
      REAL(r8), DIMENSION(:), POINTER :: keff,keff_total,k2ef
      INTEGER(i4) :: istep,mhdits,hallits,jaits,vmhdits,viscits,
     $               teits,tiits,dbits,ndits,nphi,nstop,neoits,
     $               jmode,nmodes,nmodes_total,mpseudo,ipseust,ipseuen,
     $               nsym_bpre_band,nsym_tpre_band,ncontb,ndiscb
      INTEGER(i4) :: n_dt_increase=0
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: mps_block,nindex,
     $             nindex_total,ipolst_block,ipolen_block,mpsq_block,
     $             ipqst_block,ipqen_block
      CHARACTER(32) :: integrand_flag
      LOGICAL :: b0_changed=.false.,p0_changed=.false.,
     $           nl_changed=.false.,n0_changed=.false.,
     $           eta_changed=.false.,kpll_changed=.false.,
     $           v0_changed=.false.
      LOGICAL :: impladv=.false.,threedeta=.false.,disc_nd=.false.,
     $           closure_n0_only=.true.

      CONTAINS
c-----------------------------------------------------------------------
c     this routine is used to set constant coefficients for two-fluid
c     contributions in different equations.  defining them in global
c     ensures consistency.
c-----------------------------------------------------------------------
      SUBROUTINE set_2fl_coefs(ohmslaw,advect,sep_pe,me,meomi,elemq)

      CHARACTER(*), INTENT(IN) :: ohmslaw,advect
      LOGICAL, INTENT(IN) :: sep_pe
      REAL(r8), INTENT(IN) :: me,meomi,elemq
c-----------------------------------------------------------------------
c     defaults represent MHD modeling.
c-----------------------------------------------------------------------
      coefhll=0._r8
      coefgpe=0._r8
      coefme1=0._r8
      coefme2=0._r8
      coefjvi=0._r8
      coefjve=0._r8
c-----------------------------------------------------------------------
c     coefficients for the Hall, grad-pe, and partial(J)/partial(t)
c     terms in Ohm's law.
c-----------------------------------------------------------------------
      SELECT CASE(ohmslaw)
      CASE('hall','mhd&hall')
        coefhll=(1._r8-meomi)/(elemq*(1._r8+meomi))
        coefgpe=1._r8/(elemq*(1._r8+meomi))
      CASE('2fl')
        coefhll=(1._r8-meomi)/(elemq*(1._r8+meomi))
        coefgpe=1._r8/(elemq*(1._r8+meomi))
        coefme1=me/(elemq**2*(1._r8+meomi))
      END SELECT
c-----------------------------------------------------------------------
c     coefficient for the advective part of electron inertia.
c-----------------------------------------------------------------------
      IF (ohmslaw=='2fl'.AND.advect=='all')
     $  coefme2=me/(elemq**2*(1._r8+meomi))
c-----------------------------------------------------------------------
c     coefficients for finding separate species flow velocities from
c     COM-V and J, which depend on Ohm's law to achieve energy
c     conservation.
c-----------------------------------------------------------------------
      IF (sep_pe) THEN
        SELECT CASE(ohmslaw)
        CASE('hall','mhd&hall')
          coefjve=-1._r8/(elemq*(1._r8+meomi))
        CASE('2fl')
          coefjve=-1._r8/(elemq*(1._r8+meomi))
          coefjvi=meomi/(elemq*(1._r8+meomi))
        END SELECT
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE set_2fl_coefs
c-----------------------------------------------------------------------
c     close module
c-----------------------------------------------------------------------
      END MODULE global
