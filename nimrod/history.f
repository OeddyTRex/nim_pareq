c-----------------------------------------------------------------------
c     file history.f:  module that produces time histories of data
c     at a point on the grid and integrated data across the grid.
c     the binary output is intended for Alan Glasser's xdraw graphics
c     package.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. hist_mod.
c     1. probe_hist.
c     2. probe_hist_init.
c     3. energies.
c     4. discharge_hist.
c     5. divb_check.
c     6. current_flux.
c-----------------------------------------------------------------------
c     module declaration.
c-----------------------------------------------------------------------
      MODULE hist_mod
      USE local
      USE fields
      USE input
      USE global
      USE mpi_nim
      USE pardata
      USE time
      IMPLICIT NONE

c saved variables for history routines
c initially set by history_init, used by history
c histbl = local block # (1-nrbl) that contains history point
c histx,histy = coords within histbl of history point
c history_node = processor id of owner of history point (-1 = no owner)

      INTEGER(i4) :: histbl
      INTEGER(i4) :: histx,histy
      INTEGER(i4) :: history_node = -1_i4
      REAL(r8) :: total_energy,internal_energy,internal_ien,internal_een
      REAL(r8), DIMENSION(:), ALLOCATABLE :: kinetic_energy,
     $          magnetic_energy

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. probe_hist.
c     processor who owns probe history grid point writes its data vs.
c     time to the history file.
c-----------------------------------------------------------------------
      SUBROUTINE probe_hist(close_nh)

      LOGICAL, INTENT(IN) :: close_nh

      INTEGER(i4) :: ierror,jhbl,jhx,jhy,im,imt
      CHARACTER(6), SAVE :: hist_pos='REWIND'
      REAL(r8), DIMENSION(6,nmodes_total) :: be,ja,ve,tmp
      REAL(r8), DIMENSION(2,nmodes_total) :: pr,pe,nd,co,tmps,te,ti
      LOGICAL, SAVE :: nh_opened=.false.
c-----------------------------------------------------------------------
c     open and close the file according to input.
c-----------------------------------------------------------------------
c     CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(timestart)

      node_if: IF (node == history_node) THEN

        IF (ilayer==0.AND..NOT.nh_opened) THEN
          IF (hist_binary) THEN
            CALL open_bin(hst_unit,'nimhist.bin','UNKNOWN',
     $        hist_pos,32_i4)
          ELSE
            OPEN(UNIT=hst_unit,FILE='nimhist.txt',STATUS='UNKNOWN',
     $        POSITION=hist_pos)
          ENDIF
          hist_pos='APPEND'
          nh_opened=.true.
        ENDIF

        jhbl=histbl
        jhx=histx
        jhy=histy
c-----------------------------------------------------------------------
c       communicate accross layers.
c-----------------------------------------------------------------------

        be=0
        ja=0
        ve=0
        pr=0
        pe=0
        nd=0
        co=0
        te=0
        ti=0
        DO im=1,nmodes
          imt=im+mode_lo-1
          be(1:3,imt)=rb(jhbl)%be%fs(:,jhx,jhy,im)
          be(4:6,imt)=-(0,1)*rb(jhbl)%be%fs(:,jhx,jhy,im)
          ja(1:3,imt)=rb(jhbl)%ja%fs(:,jhx,jhy,im)
          ja(4:6,imt)=-(0,1)*rb(jhbl)%ja%fs(:,jhx,jhy,im)
          ve(1:3,imt)=rb(jhbl)%ve%fs(:,jhx,jhy,im)
          ve(4:6,imt)=-(0,1)*rb(jhbl)%ve%fs(:,jhx,jhy,im)
          pr(1,imt)=rb(jhbl)%pres%fs(1,jhx,jhy,im)
          pr(2,imt)=-(0,1)*rb(jhbl)%pres%fs(1,jhx,jhy,im)
          pe(1,imt)=rb(jhbl)%prese%fs(1,jhx,jhy,im)
          pe(2,imt)=-(0,1)*rb(jhbl)%prese%fs(1,jhx,jhy,im)
          nd(1,imt)=rb(jhbl)%nd%fs(1,jhx,jhy,im)
          nd(2,imt)=-(0,1)*rb(jhbl)%nd%fs(1,jhx,jhy,im)
          co(1,imt)=rb(jhbl)%conc%fs(1,jhx,jhy,im)
          co(2,imt)=-(0,1)*rb(jhbl)%conc%fs(1,jhx,jhy,im)
          te(1,imt)=rb(jhbl)%tele%fs(1,jhx,jhy,im)
          te(2,imt)=-(0,1)*rb(jhbl)%tele%fs(1,jhx,jhy,im)
          ti(1,imt)=rb(jhbl)%tion%fs(1,jhx,jhy,im)
          ti(2,imt)=-(0,1)*rb(jhbl)%tion%fs(1,jhx,jhy,im)
        ENDDO

        IF (nlayers>1) THEN
          CALL mpi_allreduce(be,tmp,6*nmodes_total,
     $         mpi_nim_real,mpi_sum,comm_mode,ierror)
          be=tmp
          CALL mpi_allreduce(ja,tmp,6*nmodes_total,
     $         mpi_nim_real,mpi_sum,comm_mode,ierror)
          ja=tmp
          CALL mpi_allreduce(ve,tmp,6*nmodes_total,
     $         mpi_nim_real,mpi_sum,comm_mode,ierror)
          ve=tmp
          CALL mpi_allreduce(pr,tmps,2*nmodes_total,
     $         mpi_nim_real,mpi_sum,comm_mode,ierror)
          pr=tmps
          CALL mpi_allreduce(pe,tmps,2*nmodes_total,
     $         mpi_nim_real,mpi_sum,comm_mode,ierror)
          pe=tmps
          CALL mpi_allreduce(nd,tmps,2*nmodes_total,
     $         mpi_nim_real,mpi_sum,comm_mode,ierror)
          nd=tmps
          CALL mpi_allreduce(co,tmps,2*nmodes_total,
     $         mpi_nim_real,mpi_sum,comm_mode,ierror)
          co=tmps
          CALL mpi_allreduce(te,tmps,2*nmodes_total,
     $         mpi_nim_real,mpi_sum,comm_mode,ierror)
          te=tmps
          CALL mpi_allreduce(ti,tmps,2*nmodes_total,
     $         mpi_nim_real,mpi_sum,comm_mode,ierror)
          ti=tmps
        ENDIF
       
c-----------------------------------------------------------------------
c       array constructors are needed to avoid a write error on the c90
c       and j90.
c-----------------------------------------------------------------------
        write_if: IF (ilayer==0) THEN

        IF (hist_binary) THEN

          DO im=1,nmodes_total
            WRITE(hst_unit) (/
     $        REAL(istep,4),REAL(t,4),REAL(im,4),REAL(keff_total(im),4),
     $        REAL(be(:,im),4),REAL(ja(:,im),4),REAL(ve(:,im),4),
     $        REAL(pr(:,im),4),REAL(pe(:,im),4),REAL(nd(:,im),4),
     $        REAL(co(:,im),4),REAL(te(:,im),4),REAL(ti(:,im),4) /)
          ENDDO
          WRITE(hst_unit)
          IF (close_nh) THEN
            CALL close_bin(hst_unit,'nimhist.bin')
            nh_opened=.false.
          ENDIF

        ELSE

          DO im=1,nmodes_total
            WRITE(hst_unit,*)
     $        REAL(istep,4),REAL(t,4),REAL(im,4),REAL(keff_total(im),4),
     $        REAL(be(:,im),4),REAL(ja(:,im),4),REAL(ve(:,im),4),
     $        REAL(pr(:,im),4),REAL(pe(:,im),4),REAL(nd(:,im),4),
     $        REAL(co(:,im),4),REAL(te(:,im),4),REAL(ti(:,im),4)
          ENDDO
          WRITE(hst_unit,*) ' '
          IF (close_nh) THEN
            CLOSE(UNIT=hst_unit)
            nh_opened=.false.
          ENDIF

        ENDIF

        ENDIF write_if

      ENDIF node_if


c     CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(timeend)
      time_io = time_io + timeend-timestart
      
      RETURN
      END SUBROUTINE probe_hist
c-----------------------------------------------------------------------
c     subprogram 2. probe_hist_init.
c     determine which single grid point is the history grid point
c     and what processor owns it assumes a rectangular, regular grid
c     blocking
c-----------------------------------------------------------------------
      SUBROUTINE probe_hist_init

      INTEGER(i4) :: ibx,iby,ibl,ierror,nsearch_limit
      LOGICAL :: found

c     CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(timestart)

c-----------------------------------------------------------------------
c     search for point horizontally and vertically across global grid 
c     rblocks using block_sizes data structure
c-----------------------------------------------------------------------

      nsearch_limit=nrbl_total

      found = .false.
      histx = ihist
      DO ibx = 1,nsearch_limit,nybl
        IF (histx <= block_sizes(2,ibx)) THEN
          found = .true.
          EXIT
        ENDIF
        histx = histx - block_sizes(2,ibx)
      ENDDO

      IF (.NOT.found) THEN
        CALL nim_stop('History:  Point is not on the grid.')
      ENDIF

      found = .false.
      histy = jhist
      DO iby = 1,nybl
        IF (histy <= block_sizes(3,iby)) THEN
          found = .true.
          EXIT
        ENDIF
        histy = histy - block_sizes(3,iby)
      ENDDO

      IF (.NOT.found) THEN
        CALL nim_stop('History:  Point is not on the grid.')
      ENDIF
c-----------------------------------------------------------------------
c     convert global grid block to local one on a particular processor
c-----------------------------------------------------------------------
      ibl = ibx + iby - 1
      histbl = global2local(ibl)
      history_node = block2proc(ibl)

c     CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(timeend)
      time_io = time_io + timeend-timestart

      RETURN
      END SUBROUTINE probe_hist_init
c-----------------------------------------------------------------------
c     subprogram 3. energies.
c     compute the energies associated with each mode.
c-----------------------------------------------------------------------
      SUBROUTINE energies(close_en)
      USE rblock
      USE tblock
      USE diagnostic_ints
      USE computation_pointers

      LOGICAL, INTENT(IN) :: close_en

      INTEGER(i4) :: ibl,ierror,imode,iq
      REAL(r8), DIMENSION(nmodes_total) :: tmp
      REAL(r8) :: tmps
      CHARACTER(6), SAVE :: hist_pos='REWIND'
      LOGICAL, SAVE :: en_opened=.false.
c-----------------------------------------------------------------------
c     allocate arrays if this is the first call to energies.
c-----------------------------------------------------------------------
      IF (.NOT.ALLOCATED(magnetic_energy)) THEN
        ALLOCATE(magnetic_energy(nmodes_total))
        ALLOCATE(kinetic_energy(nmodes_total))
      ENDIF
      magnetic_energy=0
      internal_energy=0
      internal_ien=0
      internal_een=0
      kinetic_energy=0
c-----------------------------------------------------------------------
c     accumulate integral[b**2/2mu0] and integral[ro*v**2/2] for each
c     mode and integral[p/(gamma-1)] for n=0.  for nonlinear runs, this
c     includes the equilibrium fields, but for linear runs it doesn't.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        mpseudo=mpsq_block(ibl)
        ipseust=ipqst_block(ibl)
        ipseuen=ipqen_block(ibl)
        IF (ibl<=nrbl) THEN
          CALL rblock_get_rhs(rb(ibl),cell_crhs(ibl),energy_density,
     $                        4_i4,nmodes)
        ELSE
          CALL tblock_get_rhs(tb(ibl),cell_crhs(ibl),energy_density,
     $                        4_i4,nmodes)
        ENDIF
        DO imode=mode_lo,mode_hi
          magnetic_energy(imode)=magnetic_energy(imode)
     $       +SUM(cell_crhs(ibl)%arri(1,1,:,:,imode-mode_lo+1))
          kinetic_energy(imode)=kinetic_energy(imode)
     $       +SUM(cell_crhs(ibl)%arri(2,1,:,:,imode-mode_lo+1))
          internal_een=internal_een
     $       +SUM(cell_crhs(ibl)%arri(3,1,:,:,imode-mode_lo+1))
          internal_ien=internal_ien
     $       +SUM(cell_crhs(ibl)%arri(4,1,:,:,imode-mode_lo+1))
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     accumulate the sums from all processors.
c-----------------------------------------------------------------------
      IF (nprocs>1) THEN
        CALL mpi_allreduce(magnetic_energy,tmp,nmodes_total,
     $       mpi_nim_real,mpi_sum,mpi_comm_world,ierror)
        magnetic_energy=tmp
        CALL mpi_allreduce(kinetic_energy ,tmp,nmodes_total,
     $       mpi_nim_real,mpi_sum,mpi_comm_world,ierror)
        kinetic_energy =tmp
        CALL mpi_allreduce(internal_een,tmps,1,mpi_nim_real,
     $       mpi_sum,mpi_comm_world,ierror)
        internal_een=tmps
        CALL mpi_allreduce(internal_ien,tmps,1,mpi_nim_real,
     $       mpi_sum,mpi_comm_world,ierror)
        internal_ien=tmps
      ENDIF
c-----------------------------------------------------------------------
c     when a nonlinear computation has continuity=full, the complete
c     3D computation of kinetic energy is stored in the n=0 index,
c     and other contributions must be subtracted.  since the index
c     is over all components, the n=0 is always first.
c-----------------------------------------------------------------------
      IF (nonlinear.AND.continuity=='full') 
     $  kinetic_energy(1)=kinetic_energy(1)-SUM(kinetic_energy(2:))
c-----------------------------------------------------------------------
c     total energy:
c-----------------------------------------------------------------------
      internal_energy=internal_een+internal_ien
      total_energy=SUM(magnetic_energy+kinetic_energy)+internal_energy
c-----------------------------------------------------------------------
c     geometric factor:
c-----------------------------------------------------------------------
      IF (geom=='lin') THEN
        magnetic_energy=magnetic_energy*per_length
        kinetic_energy=kinetic_energy*per_length
        internal_energy=internal_energy*per_length
        internal_een=internal_een*per_length
        internal_ien=internal_ien*per_length
        total_energy=total_energy*per_length
      ELSE
        magnetic_energy=magnetic_energy*twopi
        kinetic_energy=kinetic_energy*twopi
        internal_energy=internal_energy*twopi
        internal_een=internal_een*twopi
        internal_ien=internal_ien*twopi
        total_energy=total_energy*twopi
      ENDIF
c-----------------------------------------------------------------------
c     output:  open and close the file according to input.
c-----------------------------------------------------------------------
c     CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(timestart)

      IF (node == 0) THEN
        IF (.NOT.en_opened) THEN
          IF (hist_binary) THEN
            CALL open_bin(en_unit,'energy.bin','UNKNOWN',hist_pos,32_i4)
          ELSE
            OPEN(UNIT=en_unit,FILE='energy.txt',STATUS='UNKNOWN',
     $           POSITION=hist_pos)
          ENDIF
          hist_pos='APPEND'
          en_opened=.true.
        ENDIF

c-----------------------------------------------------------------------
c       array constructors are needed to avoid a write error on the c90
c       and j90.
c-----------------------------------------------------------------------
        IF (hist_binary) THEN
          DO iq=1,nmodes_total
            WRITE(en_unit) (/ REAL(istep,4),REAL(t,4),REAL(iq,4),
     $                        REAL(keff_total(iq),4),
     $                        REAL(magnetic_energy(iq),4),
     $                        REAL(kinetic_energy (iq),4) /)
          ENDDO
          WRITE(en_unit)
          IF (close_en) THEN
            CALL close_bin(en_unit,'energy.bin')
            en_opened=.false.
          ENDIF
        ELSE
          DO iq=1,nmodes_total
            WRITE(en_unit,*) REAL(istep,4),REAL(t,4),REAL(iq,4),
     $                       REAL(keff_total(iq),4),
     $                       REAL(magnetic_energy(iq),4),
     $                       REAL(kinetic_energy (iq),4)
          ENDDO
          WRITE(en_unit,*) ' '
          IF (close_en) THEN
            CLOSE(UNIT=en_unit)
            en_opened=.false.
          ENDIF
        ENDIF
      ENDIF

c     CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(timeend)
      time_io = time_io + timeend-timestart
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE energies
c-----------------------------------------------------------------------
c     subprogram 4. discharge_hist.
c     compute and write global discharge parameters to an xdraw file.
c     this should be called after energies, divb_check and current_flux.
c-----------------------------------------------------------------------
      SUBROUTINE discharge_hist(close_dis)
      USE physdat
      USE global

      LOGICAL, INTENT(IN) :: close_dis

      INTEGER(i4) :: ierror
      CHARACTER(6), SAVE :: hist_pos='REWIND'
      REAL(r8), SAVE :: ek_old=-1,t_old
      REAL(r8) :: log_etot,log_ei,log_ek,growth=0,ran
      LOGICAL, SAVE :: dis_opened=.false.
c-----------------------------------------------------------------------
c     output:  open and close the file according to input.
c-----------------------------------------------------------------------
c     CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(timestart)

      IF (node == 0) THEN
        IF (.NOT.dis_opened) THEN
          IF (hist_binary) THEN
            CALL open_bin(dis_unit,'discharge.bin','UNKNOWN',
     $                    hist_pos,32_i4)
          ELSE
            OPEN(UNIT=dis_unit,FILE='discharge.txt',STATUS='UNKNOWN',
     $           POSITION=hist_pos)
          ENDIF
          hist_pos='APPEND'
          dis_opened=.true.
        ENDIF
c-----------------------------------------------------------------------
c       array constructors are needed to avoid a write error on the c90
c       and j90.
c       random number keeps xdraw from tick-crashing.
c-----------------------------------------------------------------------
        CALL RANDOM_NUMBER(ran)
        ran=1+ran/1000
        log_etot=LOG(MAX(total_energy,ran*TINY(log_etot)))
        log_ei  =LOG(MAX(internal_energy,ran*TINY(log_ei)))
        log_ek  =LOG(MAX(kinetic_energy(1),ran*TINY(log_ek)))
        IF (ek_old>0.AND.t>t_old)
     $    growth=(log_ek-LOG(ek_old))/(2*(t-t_old))
        ek_old=MAX(kinetic_energy(1),ran*TINY(log_ek))
        t_old=t
        IF (hist_binary) THEN
          WRITE(dis_unit) (/
     $      REAL(istep,4),REAL(t,4),
     $      REAL(kdivb_2,4),
     $      REAL(total_energy,4),REAL(internal_energy,4),
     $      REAL(internal_een,4),REAL(internal_ien,4),
     $      REAL(log_etot,4),REAL(log_ei,4),REAL(growth,4),
     $      REAL(i_eq+i_n0,4),REAL(i_n0,4),REAL(volt,4),
     $      REAL(flux_eq+flux_n0,4),REAL(flux_n0,4),
     $      REAL(ff,4),REAL(theta,4),
     $      REAL(lin_cfl,4),REAL(nl_cfl,4),REAL(fl_cfl,4) /)
          IF (close_dis) THEN
            CALL close_bin(dis_unit ,'discharge.bin' )
            dis_opened=.false.
          ENDIF
        ELSE
          WRITE(dis_unit,*)
     $      REAL(istep,4),REAL(t,4),
     $      REAL(kdivb_2,4),
     $      REAL(total_energy,4),REAL(internal_energy,4),
     $      REAL(internal_een,4),REAL(internal_ien,4),
     $      REAL(log_etot,4),REAL(log_ei,4),REAL(growth,4),
     $      REAL(i_eq+i_n0,4),REAL(i_n0,4),REAL(volt,4),
     $      REAL(flux_eq+flux_n0,4),REAL(flux_n0,4),
     $      REAL(ff,4),REAL(theta,4),
     $      REAL(lin_cfl,4),REAL(nl_cfl,4),REAL(fl_cfl,4)
          IF (close_dis) THEN
            CLOSE(UNIT=dis_unit )
            dis_opened=.false.
          ENDIF
        ENDIF
      ENDIF

c     CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(timeend)
      time_io = time_io + timeend-timestart
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE discharge_hist
c-----------------------------------------------------------------------
c     subprogram 5. divb_check.
c     determine the wavenumber associated with divergence(B), and
c     gracefully stop the simulation if necessary.
c-----------------------------------------------------------------------
      SUBROUTINE divb_check
      USE rblock
      USE tblock
      USE diagnostic_ints
      USE computation_pointers

      INTEGER(i4) :: ibl,ierror
      REAL(r8) :: bsq,tmp
c-----------------------------------------------------------------------
c     find integral[div(b)**2] and integral[b**2]
c-----------------------------------------------------------------------
      bsq=0
      kdivb_2=0
      DO ibl=1,nbl
        IF (ibl<=nrbl) THEN
          CALL rblock_get_rhs(rb(ibl),cell_rhs(ibl),div_b,2_i4)
        ELSE
          CALL tblock_get_rhs(tb(ibl),cell_rhs(ibl),div_b,2_i4)
        ENDIF
        kdivb_2=kdivb_2+SUM(cell_rhs(ibl)%arri(1,1,:,:))
        bsq    =    bsq+SUM(cell_rhs(ibl)%arri(2,1,:,:))
      ENDDO
c-----------------------------------------------------------------------
c     accumulate the sums from all processors.
c-----------------------------------------------------------------------
      IF (nprocs>1) THEN
        CALL mpi_allreduce(kdivb_2,tmp,1,mpi_nim_real,mpi_sum,
     $       mpi_comm_world,ierror)
        kdivb_2=tmp
        CALL mpi_allreduce(bsq,tmp,1,mpi_nim_real,mpi_sum,
     $       mpi_comm_world,ierror)
        bsq=tmp
      ENDIF
c-----------------------------------------------------------------------
c     find the ratio for an average k**2.
c-----------------------------------------------------------------------
      IF (bsq>1.e-30_r8) THEN
        kdivb_2=kdivb_2/bsq
      ELSE
        kdivb_2=0
      ENDIF
c-----------------------------------------------------------------------
c     compare with the specified limit.
c-----------------------------------------------------------------------
      IF (kdivb_2>kdivb_2_limit) THEN
        CALL nim_output(.false.)
        CALL nim_stop('Divb_check:  kdivb_2 > kdivb_2_limit.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE divb_check
c-----------------------------------------------------------------------
c     subprogram 6. current_flux.
c     compute current and magnetic flux.
c-----------------------------------------------------------------------
      SUBROUTINE current_flux
      USE rblock
      USE tblock
      USE diagnostic_ints
      USE computation_pointers
      USE seam_storage_mod

      INTEGER(i4) :: ibl,ierror,mxb,myb,mc,ibe,ix,iy,ixp,iyp,iv,nv,ivp,
     $               iseg
      REAL(r8) :: rbphi_av,circumf,dx,dy,dl,x,y
      REAL(r8), DIMENSION(2) :: tmpi,tmpo
      LOGICAL, SAVE :: eq_computed=.false.
c-----------------------------------------------------------------------
c     find the current and flux associated with the equilibrium fields.
c-----------------------------------------------------------------------
      IF (.NOT.eq_computed) THEN
        i_eq=0
        flux_eq=0
        DO ibl=1,nbl
          IF (ibl<=nrbl) THEN
            CALL rblock_get_rhs(rb(ibl),cell_rhs(ibl),eq_i_phi,2_i4)
          ELSE
            CALL tblock_get_rhs(tb(ibl),cell_rhs(ibl),eq_i_phi,2_i4)
          ENDIF
          i_eq=i_eq+SUM(cell_rhs(ibl)%arri(1,1,:,:))
          flux_eq=flux_eq+SUM(cell_rhs(ibl)%arri(2,1,:,:))
        ENDDO
c-----------------------------------------------------------------------
c       accumulate the sums from all processors.
c-----------------------------------------------------------------------
        IF (nprocs>1) THEN
          tmpi=(/i_eq,flux_eq/)
          CALL mpi_allreduce(tmpi,tmpo,2,mpi_nim_real,
     $         mpi_sum,comm_layer,ierror)
          i_eq=tmpo(1)
          flux_eq=tmpo(2)
        ENDIF
        eq_computed=.true.
      ENDIF
c-----------------------------------------------------------------------
c     find the current and flux for the n=0 Fourier component for 
c     nonlinear runs only.  n=0 components are always in layer 0
c     for nonlinear cases.
c-----------------------------------------------------------------------
      i_n0_old=i_n0
      i_n0=0
      flux_n0=0
      IF (nonlinear) THEN
        IF (ilayer==0) THEN
          DO ibl=1,nbl
            IF (ibl<=nrbl) THEN
              CALL rblock_get_rhs(rb(ibl),cell_rhs(ibl),n0_i_phi,2_i4)
            ELSE
              CALL tblock_get_rhs(tb(ibl),cell_rhs(ibl),n0_i_phi,2_i4)
            ENDIF
            i_n0=i_n0+SUM(cell_rhs(ibl)%arri(1,1,:,:))
            flux_n0=flux_n0+SUM(cell_rhs(ibl)%arri(2,1,:,:))
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       accumulate the sums from all processors.
c-----------------------------------------------------------------------
        IF (nprocs>1) THEN
          tmpi=(/i_n0,flux_n0/)
          CALL mpi_allreduce(tmpi,tmpo,2,mpi_nim_real,
     $         mpi_sum,mpi_comm_world,ierror)
          i_n0=tmpo(1)
          flux_n0=tmpo(2)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     compute ff and theta, the reversal and pinch parameters,
c     respectively.  for an arbitrary cross section in toroidal
c     geometry, we define ff as 
c     mu0*(shell current)*integral(dA/r) / (2*pi*flux), which is
c     equivalent to ave(R*B_phi)*integral(dA/r) / flux, where the
c     average is based on a loop integral around the poloidal cross
c     section for the toroidal average B_phi.
c
c-PRE tblock higher order computation?
c-----------------------------------------------------------------------
      rbphi_av=0
      circumf=0
      IF (ilayer==0) THEN
        DO ibl=1,SIZE(exblock_list)
          ibe=exblock_list(ibl)
          nv=seam(ibe)%nvert
          DO iv=1,nv
            ivp=iv-1
            IF (ivp<1) ivp=nv
            IF (seam(ibe)%expoint(iv).AND.seam(ibe)%expoint(ivp)) THEN
              ix=seam(ibe)%vertex(iv)%intxy(1)
              iy=seam(ibe)%vertex(iv)%intxy(2)
              ixp=seam(ibe)%vertex(ivp)%intxy(1)
              iyp=seam(ibe)%vertex(ivp)%intxy(2)
              block_type: IF (ibe<=nrbl) THEN
                dx=REAL(ix-ixp)/poly_degree
                dy=REAL(iy-iyp)/poly_degree
                DO iseg=1,poly_degree
                  x=ixp+(iseg-0.5_r8)*dx
                  y=iyp+(iseg-0.5_r8)*dy
                  CALL lagr_quad_eval(rb(ibe)%rz,x,y,1_i4)
                  dl=SQRT(SUM((rb(ibe)%rz%fx*dx)**2
     $                       +(rb(ibe)%rz%fy*dy)**2))
                  circumf=circumf+dl
                  CALL lagr_quad_eval(rb(ibe)%be_eq,x,y,0_i4)
                  rbphi_av=rbphi_av+dl*rb(ibe)%be_eq%f(3)
                  IF (nonlinear) THEN
                    CALL lagr_quad_eval(rb(ibe)%be,x,y,0_i4)
                    IF (geom=='tor') THEN
                      rbphi_av=rbphi_av
     $                        +dl*rb(ibe)%be%f(3,1)*rb(ibe)%rz%f(1)
                    ELSE
                      rbphi_av=rbphi_av+dl*rb(ibe)%be%f(3,1)
                    ENDIF
                  ENDIF
                ENDDO
              ELSE block_type
                dl=SQRT((tb(ibe)%tgeom%xs(ix)-tb(ibe)%tgeom%xs(ixp))**2
     $                 +(tb(ibe)%tgeom%ys(ix)-tb(ibe)%tgeom%ys(ixp))**2)
                circumf=circumf+dl
                rbphi_av=rbphi_av+dl*0.5*(
     $             tb(ibe)%be_eq%fs(3,ix ,iy )
     $            +tb(ibe)%be_eq%fs(3,ixp,iyp))
                IF (nonlinear) THEN
                  IF (geom=='tor') THEN
                    rbphi_av=rbphi_av+dl*0.5*(
     $                 tb(ibe)%be%fs(3,ix ,iy ,1)*tb(ibe)%tgeom%xs(ix )
     $                +tb(ibe)%be%fs(3,ixp,iyp,1)*tb(ibe)%tgeom%xs(ixp))
                  ELSE
                    rbphi_av=rbphi_av+dl*0.5*(
     $                 tb(ibe)%be%fs(3,ix ,iy ,1)
     $                +tb(ibe)%be%fs(3,ixp,iyp,1))
                  ENDIF
                ENDIF
              ENDIF block_type
            ENDIF
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     sum over processors.
c-----------------------------------------------------------------------
      IF (nprocs>1) THEN
        tmpi=(/rbphi_av,circumf/)
        CALL mpi_allreduce(tmpi,tmpo,2,mpi_nim_real,
     $       mpi_sum,mpi_comm_world,ierror)
        rbphi_av=tmpo(1)
        circumf=tmpo(2)
      ENDIF
      IF (gridshape=='rect'.AND.periodicity=='both'
     $    .OR.flux_n0+flux_eq==0) THEN
        ff=0
        theta=0
      ELSE
        rbphi_av=rbphi_av/circumf
        ff=rbphi_av*cross_s_overr/(flux_n0+flux_eq)
        theta=mu0*(i_n0+i_eq)*cross_section/(circumf*(flux_n0+flux_eq))
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE current_flux
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE hist_mod
