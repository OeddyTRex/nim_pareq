c-----------------------------------------------------------------------
c     file utilities.f:  contains utility subprograms for the physics
c     kernel.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  nim_stat.
c     2.  nim_version.
c     3.  nim_output.
c     4.  nim_stop.
c     5.  new_dt.
c     6.  add_mass.
c     7.  add_mass_comp.
c     8.  iter_out.
c     9.  loop_voltage.
c     10. vertical_efield.
c     11. qp_fft_save.
c     12. qp_fft_noeq_save.
c     13. soln_save.
c     14. zero_last_comp.
c     15. fourier_damp.
c     16. set_nodal_min.
c-----------------------------------------------------------------------
c     subprogram 1. nim_stat.
c     writes time-step status information to standard out and 
c     nimrod.out
c-----------------------------------------------------------------------
      SUBROUTINE nim_stat(close_wr)
      USE local
      USE global
      USE pardata
      USE input
      USE time
      USE mpi_nim
      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: close_wr

      INTEGER(i4) :: ierror
      INTEGER(i4), PARAMETER :: nits=10
      INTEGER(i4), DIMENSION(nits) :: tmp
c-----------------------------------------------------------------------
c     format definition
c-----------------------------------------------------------------------
  100 FORMAT (/,' Cycle = ',i6,' Time = ',es12.5,' dt = ',es12.5,/,
     $        '   BMHD its  =',i4,  ' Hall its  =',i4,' Div_B its =',i4,
     $        /,'   VMHD its  =',i4,' Visc its  =',i4,
     $        ' Te its    =',i4,    ' Ti its    =',i4,
     $        /,'   Ja its    =',i4,' Neo its   =',i4,' Cont its  =',i4,
     $        /,'   Flow CFL  =',es10.3,' Nonlinear CFL =',es10.3,
     $        ' k_divb**2 =',es10.3)
c-----------------------------------------------------------------------
c     accumulate iteration count over multiple processor layers.
c-----------------------------------------------------------------------
      IF (nlayers>1) THEN
        CALL mpi_allreduce
     $    ((/mhdits,hallits,dbits,vmhdits,viscits,teits,tiits,
     $       jaits,neoits,ndits/),
     $     tmp,nits,mpi_nim_int,mpi_sum,comm_mode,ierror)
        mhdits =tmp(1)
        hallits=tmp(2)
        dbits  =tmp(3)
        vmhdits=tmp(4)
        viscits=tmp(5)
        teits  =tmp(6)
        tiits  =tmp(7)
        jaits  =tmp(8)
        neoits =tmp(9)
        ndits  =tmp(10)
      ENDIF
c-----------------------------------------------------------------------
c     write output.  the file is opened and closed as directed by
c     input.
c-----------------------------------------------------------------------
c     CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(timestart)
      IF (node == 0) THEN
        IF (.NOT.out_opened) THEN
          OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $         POSITION='APPEND')
          out_opened=.true.
        ENDIF
        WRITE(out_unit,100) istep,t,dt,mhdits,hallits,dbits,vmhdits,
     $     viscits,teits,tiits,jaits,neoits,ndits,fl_cfl,nl_cfl,kdivb_2
        IF (close_wr) THEN
          CLOSE(UNIT=out_unit)
          out_opened=.false.
        ENDIF
      ENDIF
c     CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(timeend)
      time_io = time_io + timeend-timestart
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE nim_stat
c-----------------------------------------------------------------------
c     subprogram 2. nim_version.
c     reads the nimrod version number from the README file and writes it
c     to nimrod.out.
c-----------------------------------------------------------------------
      SUBROUTINE nim_version
      USE local
      USE time
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      CHARACTER(128) :: version,arch_type='********'
      LOGICAL :: rstat
      INTEGER, DIMENSION(8) :: tval
      INTEGER(i4) :: ierror
c-----------------------------------------------------------------------
c     check for README file.
c-----------------------------------------------------------------------
c     CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(timestart)
      IF (node == 0) THEN
        INQUIRE(FILE='../nimrod/README',exist=rstat)
        IF (rstat) THEN
          OPEN(UNIT=in_unit,FILE='../nimrod/README',
     $         STATUS='OLD',PAD='YES')
          READ(UNIT=in_unit,FMT='(a24)') version
          READ(UNIT=in_unit,FMT='(a31)',ADVANCE='NO') version
          READ(UNIT=in_unit,FMT='(a24)') version
          CLOSE(UNIT=in_unit)
        ELSE
          version='********'
        ENDIF
c-----------------------------------------------------------------------
c       open the output file for the first time.
c-----------------------------------------------------------------------
        OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $       POSITION='REWIND')
c-----------------------------------------------------------------------
c       write version #.
c-----------------------------------------------------------------------
        WRITE(out_unit,'(2a,/)')
     $    'Output from NIMROD version # ',TRIM(version)
c-----------------------------------------------------------------------
c       simulation date and time.
c-----------------------------------------------------------------------
        CALL DATE_AND_TIME(VALUES=tval)
        WRITE(out_unit,'(a,i2,a,i2,a,i4,a,i2,a,i2,a,i2,/)')
     $    'Simulation started on ',tval(2),'/',tval(3),'/',tval(1),
     $    ' at ',tval(5),':',tval(6),':',tval(7)
c-----------------------------------------------------------------------
c       architecture type.
c-----------------------------------------------------------------------
        CALL system('uname > temporary_uname_file')
        INQUIRE(FILE='temporary_uname_file',exist=rstat)
        IF (rstat) THEN
          OPEN(UNIT=temp_unit,FILE='temporary_uname_file')
          READ(temp_unit,'(a)') arch_type
          CLOSE(temp_unit)
          CALL system('rm temporary_uname_file')
        ENDIF
        WRITE(out_unit,'(2a,/)')
     $    'Run with computer architecture type ',TRIM(arch_type)
c-----------------------------------------------------------------------
c       close the output file to force a buffer flush.
c-----------------------------------------------------------------------
        CLOSE(UNIT=out_unit)
      ENDIF
c     CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(timeend)
      time_io = time_io + timeend-timestart
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE nim_version
c-----------------------------------------------------------------------
c     subprogram 3. nim_output.
c     call graphic output and data dump routines.
c-----------------------------------------------------------------------
      SUBROUTINE nim_output(check_step)
      USE local
      USE input
      USE hist_mod
      USE diagnose
      USE pardata
      USE global
      USE dump
      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: check_step
      LOGICAL :: close_files
      LOGICAL, SAVE :: close_next=.false.
      INTEGER(i4) :: ierror
c-----------------------------------------------------------------------
c     always compute currents; they are now used during the advance. 
c-----------------------------------------------------------------------
      CALL current_flux
c-----------------------------------------------------------------------
c     write output if istep is at the end of a plotting cycle.
c-----------------------------------------------------------------------
      IF (check_step) THEN
        IF (ndump>0.AND.MODULO(istep,ndump)==0.AND.istep>nstop-nstep)
     $    close_next=.true.
        IF (nhist > 0 .AND. MODULO(istep,nhist) == 0) THEN
          SELECT CASE(hist_flush)
          CASE("always")
            close_files=.true.
          CASE("at dumps")
            IF (close_next) THEN
              close_files=.true.
              close_next=.false.
            ELSE
              close_files=.false.
            ENDIF
          CASE DEFAULT
            close_files=.false.
          END SELECT
          CALL energies(close_files)
          CALL probe_hist(close_files)
          CALL discharge_hist(close_files)
          CALL nim_stat(close_files)
        ENDIF
c-----------------------------------------------------------------------
c       slices across the grid for single processor runs.  parallel
c       computing requires graphical post-processing from the dump file.
c-----------------------------------------------------------------------
        IF (nprocs==1) THEN
          IF (xy_stride > 0) THEN
            IF (MODULO(istep,xy_stride) == 0) CALL xy_slice(istep,t)
          ENDIF
          IF (xt_stride > 0) THEN
            IF (MODULO(istep,xt_stride) == 0)
     $        CALL xt_slice(y0fac,istep,t)
          ENDIF
          IF (yt_stride > 0) THEN
            IF (MODULO(istep,yt_stride) == 0)
     $        CALL yt_slice(x0fac,istep,t)
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       write restart dump.  nodal pressures are generated for plotting.
c-----------------------------------------------------------------------
        IF (ndump>0.AND.MODULO(istep,ndump)==0.AND.istep>nstop-nstep)
     $    THEN
          CALL p_from_nt('all')
          CALL dump_write(nmodes,nmodes_total,keff_total,t,istep,
     $                    dump_name)
        ENDIF
c-----------------------------------------------------------------------
c     write output regardless of stride--usually done before stopping.
c     time hisories are first.
c-----------------------------------------------------------------------
      ELSE
        IF (nhist>0 .AND. istep>nstop-nstep) THEN
          CALL energies(.true.)
          CALL probe_hist(.true.)
          CALL discharge_hist(.true.)
          CALL nim_stat(.true.)
        ENDIF
c-----------------------------------------------------------------------
c       slices across the grid.
c-----------------------------------------------------------------------
        IF (nprocs==1) THEN
          IF (xy_stride>0) CALL xy_slice(istep,t)
          IF (xt_stride>0) CALL xt_slice(y0fac,istep,t)
          IF (yt_stride>0) CALL yt_slice(x0fac,istep,t)
        ENDIF
c-----------------------------------------------------------------------
c       dump file.  nodal pressures are generated for plotting.
c-----------------------------------------------------------------------
        CALL p_from_nt('all')
        CALL dump_write(nmodes,nmodes_total,keff_total,t,istep,
     $                  dump_name)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE nim_output
c-----------------------------------------------------------------------
c     subprogram 4. nim_stop.
c     closes output files, writes completion message and stops nimrod.
c-----------------------------------------------------------------------
      SUBROUTINE nim_stop(message)
      USE local
      USE input
      USE mpi_nim
      USE pardata
      USE time
      USE global
      USE diagnose
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: message
      INTEGER(i4) :: ierror
      LOGICAL, EXTERNAL :: direct_check
c-----------------------------------------------------------------------
c     for parallel runs, there are lots of cases where nim_stop is
c     called by a subset of the processors.  a less graceful exit is
c     required to stop all processes.
c-----------------------------------------------------------------------
      IF (nprocs>1.AND..NOT.(message=="Normal termination.".OR.
     $                       message=="CPU time limit reached.")) THEN
        IF (.NOT.out_opened) THEN
          OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $         POSITION='APPEND')
        ENDIF
        WRITE(out_unit,'(a,i3,2a)') 'NIM_STOP from node ',node,' => ',
     $    TRIM(message)
        CLOSE(UNIT=out_unit)
        WRITE(nim_wr,'(a,i3,2a)') 'NIM_STOP from node ',node,' => ',
     $    TRIM(message)
        CALL mpi_abort(mpi_comm_world,node,ierror)
      ENDIF
c-----------------------------------------------------------------------
c     complete and close IBM Data Explorer output and other graphics.
c-----------------------------------------------------------------------
      IF (nprocs==1) THEN
        IF (xt_stride>0 .OR. yt_stride>0) CALL time_slice_close
      ENDIF
c-----------------------------------------------------------------------
c     finish timing information.
c-----------------------------------------------------------------------
c     CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(time_total_end)
c-----------------------------------------------------------------------
c     write completion message.
c-----------------------------------------------------------------------
      IF (message=="Normal termination.".OR.
     $    message=="CPU time limit reached.") CALL timer_stats
      IF (node == 0) THEN
        IF (.NOT.out_opened) THEN
          OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $         POSITION='APPEND')
        ENDIF
        WRITE(out_unit,'(2a)') 'NIM_STOP => ', TRIM(message)
        WRITE(nim_wr,'(2a)') 'NIM_STOP => ', TRIM(message)
c-----------------------------------------------------------------------
c     close files.
c-----------------------------------------------------------------------
        CLOSE(UNIT=out_unit)
        IF (hist_binary) THEN
          CALL open_bin(dis_unit,'discharge.bin','UNKNOWN',
     $                  'APPEND',32_i4)
          WRITE(dis_unit)
          CALL close_bin(dis_unit ,'discharge.bin' )
        ENDIF
        IF (itflag) CLOSE(UNIT=it_unit)
      ENDIF
c-----------------------------------------------------------------------
c     release the processor grid used by SuperLU.
c-----------------------------------------------------------------------
      IF (direct_check(solver)) THEN
        CALL c_fortran_slugrid(2_i4,comm_layer,node_layer,slu_nrowp,
     $                         slu_ncolp,slugrid_handle)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL mpi_finalize(ierror)
      STOP
      END SUBROUTINE nim_stop
c-----------------------------------------------------------------------
c     subprogram 5. new_dt.
c     determines time step.
c-----------------------------------------------------------------------
      SUBROUTINE new_dt(converged)
      USE local
      USE physdat
      USE input
      USE global
      USE fields
      USE fft_mod
      USE mpi_nim
      USE pardata
      USE math_tran
      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: converged

      CHARACTER(64) :: msg
      INTEGER(i4) :: ibl,ix,iy,mxb,myb,ic,mc,ip,j,ierror,nph,imode,nq,
     $               iq,jq,rem,cell_per_layer,icst,icen,ig,ng,jg,icq,icg
      INTEGER(i4), DIMENSION(4) :: iv
      INTEGER(i4), SAVE :: n_same_dt=0
      REAL(r8) :: dtfl,cfl,dtnl,tmp,cl2t,avei,pol_fac
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: drz,real_v
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: real_vx,real_t
      REAL(r8), DIMENSION(:,:,:), POINTER :: real_b,real_n
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: b0,v0,v0x
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: cl2,md,kph,bigr,nd_cell
      REAL(r8), DIMENSION(3) :: bt
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: v,vx,b,teff,br,bz
      LOGICAL, SAVE :: first_step=.true.
      REAL(r8), PARAMETER :: dt_reduce=0.5
c-----------------------------------------------------------------------
c     avei is the inverse of the # of velocity datums points
c     averaged around an rblock cell.  pol_fac is the relation between
c     the effective finite element wavenumber and a cell dimensionl.
c     the sqrt(3) is accurate for linear elements, including the effect
c     of the mass matrix, but the coefficient for higher degree
c     polynomials is just approximated.
c
c     the n0 fields are now treated as bilinear.
c-----------------------------------------------------------------------
      avei=1._r8/(poly_degree+1)**2
      IF (poly_degree==1) THEN
        pol_fac=SQRT(3._r8)
      ELSE
        pol_fac=pi*poly_degree
      ENDIF
c-----------------------------------------------------------------------
c     compute the time step for a flow cfl of unity (vertex info is
c     averaged to the cell centers).  also compute cfls based on 
c     the magneto-acoustic wave speed (as if it had an isotropic
c     dispersion relation for a conservative, i.e. high, estimate
c     of the cfl) from linear and nonlinear terms separately.
c     rblocks first:
c-----------------------------------------------------------------------
      IF (.NOT.first_step) dt_old=dt
      lin_cfl=0
      nl_cfl=0
      cfl=TINY(0._r8)
      DO ibl=1,nrbl
        mxb=rb(ibl)%mx
        myb=rb(ibl)%my
        ng=rb(ibl)%ng
        ig=MAX(1_i4,ng/2)
        ALLOCATE(drz(2,2,mxb,myb),kph(mxb,myb),bigr(mxb,myb),
     $           nd_cell(mxb,myb))
        drz(:,1,:,:)=0.5*(rb(ibl)%rz%fs(:,1:mxb  ,0:myb-1)
     $                   +rb(ibl)%rz%fs(:,1:mxb  ,1:myb  )
     $                   -rb(ibl)%rz%fs(:,0:mxb-1,0:myb-1)
     $                   -rb(ibl)%rz%fs(:,0:mxb-1,1:myb  ))
        drz(:,2,:,:)=0.5*(rb(ibl)%rz%fs(:,0:mxb-1,1:myb  )
     $                   +rb(ibl)%rz%fs(:,1:mxb  ,1:myb  )
     $                   -rb(ibl)%rz%fs(:,0:mxb-1,0:myb-1)
     $                   -rb(ibl)%rz%fs(:,1:mxb  ,0:myb-1))
        nd_cell=RESHAPE(rb(ibl)%qnd_eq%qpf(1,ig,:),(/mxb,myb/))
        IF (nonlinear.AND.
     $      (continuity=='n=0 only'.OR.continuity=='full'))
     $    nd_cell=nd_cell+RESHAPE(rb(ibl)%qnd_n0%qpf(1,ig,:),
     $                            (/mxb,myb/))
        IF (geom=='tor') THEN
          bigr=0.25*(rb(ibl)%rz%fs(1,1:mxb  ,0:myb-1)
     $              +rb(ibl)%rz%fs(1,1:mxb  ,1:myb  )
     $              +rb(ibl)%rz%fs(1,0:mxb-1,0:myb-1)
     $              +rb(ibl)%rz%fs(1,0:mxb-1,1:myb  ))
          kph=MAXVAL(ABS(keff_total))/bigr
        ELSE
          bigr=1
          kph=MAXVAL(ABS(keff_total))
        ENDIF
c-----------------------------------------------------------------------
c       flow cfl in rblocks--find the maximum value of v.x/x**2, where
c       x is each of the cell dimensions.  this is now a parallel
c       computation.
c-----------------------------------------------------------------------
        IF (nonlinear.OR.eq_flow/='none') THEN
          IF (nonlinear) THEN
            ALLOCATE(v(3,mxb,myb,nmodes),vx(3,mxb,myb,nmodes))
            DO imode=1,nmodes
              v(:,:,:,imode)=
     $          avei*(ve(ibl)%arr(:,1:mxb  ,0:myb-1,imode)
     $               +ve(ibl)%arr(:,1:mxb  ,1:myb  ,imode)
     $               +ve(ibl)%arr(:,0:mxb-1,0:myb-1,imode)
     $               +ve(ibl)%arr(:,0:mxb-1,1:myb  ,imode))
              IF (poly_degree>1) THEN
                v(:,:,:,imode)=v(:,:,:,imode)+
     $            avei*SUM(ve(ibl)%arrh(:,:,1:mxb,0:myb-1,imode)
     $                    +ve(ibl)%arrh(:,:,1:mxb,1:myb  ,imode)
     $                    +ve(ibl)%arrv(:,:,0:mxb-1,1:myb,imode)
     $                    +ve(ibl)%arrv(:,:,1:mxb  ,1:myb,imode),2)+
     $            avei*SUM(ve(ibl)%arri(:,:,:,:,imode),2)
              ENDIF
              IF (eq_flow/='none'.AND.keff(imode)==0) THEN
                v(:,:,:,imode)=v(:,:,:,imode)
     $            +0.25*(ve_eq(ibl)%arr(:,1:mxb  ,0:myb-1)
     $                  +ve_eq(ibl)%arr(:,1:mxb  ,1:myb  )
     $                  +ve_eq(ibl)%arr(:,0:mxb-1,0:myb-1)
     $                  +ve_eq(ibl)%arr(:,0:mxb-1,1:myb  ))
              ENDIF
              vx(1,:,:,imode)=
     $          pol_fac*SUM(v(1:2,:,:,imode)*drz(:,1,:,:),1)
     $                 /SUM(drz(:,1,:,:)**2,1)
              vx(2,:,:,imode)=
     $          pol_fac*SUM(v(1:2,:,:,imode)*drz(:,2,:,:),1)
     $                 /SUM(drz(:,2,:,:)**2,1)
              vx(3,:,:,imode)=kph*v(3,:,:,imode)
            ENDDO
            DEALLOCATE(v)
            ALLOCATE(real_vx(3,mps_block(ibl),nphi))
            CALL fft_nim('inverse',mxb*myb,mps_block(ibl),lphi,
     $                   3_i4,vx,real_vx,dealiase)
            cfl=MAX(cfl,MAXVAL(ABS(real_vx)))
            DEALLOCATE(vx,real_vx)
          ELSE
            ALLOCATE(v0(3,mxb,myb),v0x(3,mxb,myb))
            v0=0.25*(ve_eq(ibl)%arr(:,1:mxb  ,0:myb-1)
     $              +ve_eq(ibl)%arr(:,1:mxb  ,1:myb  )
     $              +ve_eq(ibl)%arr(:,0:mxb-1,0:myb-1)
     $              +ve_eq(ibl)%arr(:,0:mxb-1,1:myb  ))
            v0x(1,:,:)=pol_fac*ABS(SUM(v0(1:2,:,:)*drz(:,1,:,:),1))
     $                            /SUM(drz(:,1,:,:)**2,1)
            v0x(2,:,:)=pol_fac*ABS(SUM(v0(1:2,:,:)*drz(:,2,:,:),1))
     $                            /SUM(drz(:,2,:,:)**2,1)
            v0x(3,:,:)=kph*ABS(v0(3,:,:))
            cfl=MAX(cfl,MAXVAL(v0x))
            DEALLOCATE(v0,v0x)
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       electron flow cfl in rblocks only if electron pressure is
c       evolved separately.  the first part of this is just finding
c       a cell-centered J then V_electron.
c-----------------------------------------------------------------------
        IF (coefjve/=0._r8) THEN
          IF (nonlinear) THEN
            ALLOCATE(v(3,mxb,myb,nmodes),vx(3,mxb,myb,nmodes))
            DO imode=1,nmodes
              v(:,:,:,imode)=
     $          avei*(ja(ibl)%arr(:,1:mxb  ,0:myb-1,imode)
     $               +ja(ibl)%arr(:,1:mxb  ,1:myb  ,imode)
     $               +ja(ibl)%arr(:,0:mxb-1,0:myb-1,imode)
     $               +ja(ibl)%arr(:,0:mxb-1,1:myb  ,imode))
              IF (poly_degree>1) THEN
                v(:,:,:,imode)=v(:,:,:,imode)+
     $            avei*SUM(ja(ibl)%arrh(:,:,1:mxb,0:myb-1,imode)
     $                    +ja(ibl)%arrh(:,:,1:mxb,1:myb  ,imode)
     $                    +ja(ibl)%arrv(:,:,0:mxb-1,1:myb,imode)
     $                    +ja(ibl)%arrv(:,:,1:mxb  ,1:myb,imode),2)+
     $            avei*SUM(ja(ibl)%arri(:,:,:,:,imode),2)
              ENDIF
              IF (keff(imode)==0) THEN
                v(1:2,:,:,imode)=v(1:2,:,:,imode)
     $            +0.25*(ja_eq(ibl)%arr(1:2,1:mxb  ,0:myb-1)
     $                  +ja_eq(ibl)%arr(1:2,1:mxb  ,1:myb  )
     $                  +ja_eq(ibl)%arr(1:2,0:mxb-1,0:myb-1)
     $                  +ja_eq(ibl)%arr(1:2,0:mxb-1,1:myb  ))
                v(3,:,:,imode)=v(3,:,:,imode)
     $            +0.25*(ja_eq(ibl)%arr(3,1:mxb  ,0:myb-1)
     $                  +ja_eq(ibl)%arr(3,1:mxb  ,1:myb  )
     $                  +ja_eq(ibl)%arr(3,0:mxb-1,0:myb-1)
     $                  +ja_eq(ibl)%arr(3,0:mxb-1,1:myb  ))*bigr
              ENDIF
            ENDDO
            DO imode=1,nmodes
              DO iq=1,3
                v(iq,:,:,imode)=coefjve*v(iq,:,:,imode)/nd_cell
     $            +avei*(ve(ibl)%arr(iq,1:mxb  ,0:myb-1,imode)
     $                  +ve(ibl)%arr(iq,1:mxb  ,1:myb  ,imode)
     $                  +ve(ibl)%arr(iq,0:mxb-1,0:myb-1,imode)
     $                  +ve(ibl)%arr(iq,0:mxb-1,1:myb  ,imode))
              ENDDO
              IF (poly_degree>1) THEN
                v(:,:,:,imode)=v(:,:,:,imode)+
     $            avei*SUM(ve(ibl)%arrh(:,:,1:mxb,0:myb-1,imode)
     $                    +ve(ibl)%arrh(:,:,1:mxb,1:myb  ,imode)
     $                    +ve(ibl)%arrv(:,:,0:mxb-1,1:myb,imode)
     $                    +ve(ibl)%arrv(:,:,1:mxb  ,1:myb,imode),2)+
     $            avei*SUM(ve(ibl)%arri(:,:,:,:,imode),2)
              ENDIF
              IF (eq_flow/='none'.AND.keff(imode)==0) THEN
                v(:,:,:,imode)=v(:,:,:,imode)
     $            +0.25*(ve_eq(ibl)%arr(:,1:mxb  ,0:myb-1)
     $                  +ve_eq(ibl)%arr(:,1:mxb  ,1:myb  )
     $                  +ve_eq(ibl)%arr(:,0:mxb-1,0:myb-1)
     $                  +ve_eq(ibl)%arr(:,0:mxb-1,1:myb  ))
              ENDIF
              vx(1,:,:,imode)=
     $          pol_fac*SUM(v(1:2,:,:,imode)*drz(:,1,:,:),1)
     $                 /SUM(drz(:,1,:,:)**2,1)
              vx(2,:,:,imode)=
     $          pol_fac*SUM(v(1:2,:,:,imode)*drz(:,2,:,:),1)
     $                 /SUM(drz(:,2,:,:)**2,1)
              vx(3,:,:,imode)=kph*v(3,:,:,imode)
            ENDDO
            DEALLOCATE(v)
            ALLOCATE(real_vx(3,mps_block(ibl),nphi))
            CALL fft_nim('inverse',mxb*myb,mps_block(ibl),lphi,
     $                   3_i4,vx,real_vx,dealiase)
            cfl=MAX(cfl,MAXVAL(ABS(real_vx)))
            DEALLOCATE(vx,real_vx)
          ELSE
            ALLOCATE(v0(3,mxb,myb),v0x(3,mxb,myb))
            DO iq=1,2
              v0(iq,:,:)=-0.25*(ja_eq(ibl)%arr(iq,1:mxb  ,0:myb-1)
     $                         +ja_eq(ibl)%arr(iq,1:mxb  ,1:myb  )
     $                         +ja_eq(ibl)%arr(iq,0:mxb-1,0:myb-1)
     $                         +ja_eq(ibl)%arr(iq,0:mxb-1,1:myb  ))
     $                        /(elementary_q*nd_cell)
            ENDDO
            v0(3,:,:)=-0.25*(ja_eq(ibl)%arr(3,1:mxb  ,0:myb-1)
     $                      +ja_eq(ibl)%arr(3,1:mxb  ,1:myb  )
     $                      +ja_eq(ibl)%arr(3,0:mxb-1,0:myb-1)
     $                      +ja_eq(ibl)%arr(3,0:mxb-1,1:myb  ))*bigr
     $                     /(elementary_q*nd_cell)
            IF (eq_flow/='none') THEN
              v0=v0+0.25*(ve_eq(ibl)%arr(:,1:mxb  ,0:myb-1)
     $                   +ve_eq(ibl)%arr(:,1:mxb  ,1:myb  )
     $                   +ve_eq(ibl)%arr(:,0:mxb-1,0:myb-1)
     $                   +ve_eq(ibl)%arr(:,0:mxb-1,1:myb  ))
            ENDIF
            v0x(1,:,:)=pol_fac*ABS(SUM(v0(1:2,:,:)*drz(:,1,:,:),1))
     $                            /SUM(drz(:,1,:,:)**2,1)
            v0x(2,:,:)=pol_fac*ABS(SUM(v0(1:2,:,:)*drz(:,2,:,:),1))
     $                            /SUM(drz(:,2,:,:)**2,1)
            v0x(3,:,:)=kph*ABS(v0(3,:,:))
            cfl=MAX(cfl,MAXVAL(v0x))
            DEALLOCATE(v0,v0x)
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       magneto-acoustic speed computations in rblocks--'linear' part,
c       i.e., equilibrium + saved n=0 component first.
c-----------------------------------------------------------------------
        kph=MAX(kph**2,pol_fac**2/MINVAL(SUM(drz**2,1),1))
        ALLOCATE(cl2(mxb,myb),md(mxb,myb),b0(3,mxb,myb))
        md=mtot*nd_cell
        b0=RESHAPE(rb(ibl)%qbe_eq%qpf(:,ig,:),(/3,mxb,myb/))
        IF (nonlinear)
     $    b0=b0+RESHAPE(rb(ibl)%qbe_n0%qpf(:,ig,:),(/3,mxb,myb/))
        cl2=SUM(b0**2,1)/mu0
        IF (beta>0) THEN
          cl2=cl2+gamma*
     $      RESHAPE(rb(ibl)%qpres_eq%qpf(1,ig,:),(/mxb,myb/))
          IF (nonlinear) cl2=cl2+gamma*
     $      RESHAPE(rb(ibl)%qpres_n0%qpf(1,ig,:),(/mxb,myb/))
        ENDIF
        lin_cfl=MAX(lin_cfl,MAXVAL(cl2*kph/md))
        DEALLOCATE(cl2,nd_cell)
c-----------------------------------------------------------------------
c       magneto-acoustic speed computations in rblocks--'nonlinear'
c       part, i.e., total c**2 - 'linear' c**2 from above.  the
c       layer-parallel computation require selecting the correct portion
c       of the poloidal extent of this block for the configuration-space
c       computations.
c-----------------------------------------------------------------------
        IF (nonlinear) THEN
          icst=ipqst_block(ibl)
          icen=ipqen_block(ibl)
          real_b=>rb(ibl)%qbe_tot%qpf
          IF (beta<=0) THEN
            ic=0
            yloopr1: DO iy=1,myb
              DO ix=1,mxb
                DO jg=1,ng
                  ic=ic+1
                  IF (ic<icst.OR.jg/=ig) CYCLE
                  IF (ic>icen) EXIT yloopr1
                  icg=ic-icst+1
                  DO ip=1,nphi
                    nl_cfl=MAX(nl_cfl,kph(ix,iy)/(mu0*md(ix,iy))*
     $                ABS(SUM(real_b(:,icg,ip)**2)-SUM(b0(:,ix,iy)**2)))
                  ENDDO
                ENDDO
              ENDDO
            ENDDO yloopr1
          ELSE
            real_n=>rb(ibl)%qnd_tot%qpf
            ALLOCATE(teff(1,ng,mxb*myb,nmodes),
     $               real_t(1,mpsq_block(ibl),nphi))
            DO imode=1,nmodes
              teff(:,:,:,imode)=
     $          kboltz*(rb(ibl)%qtele%qpf(:,:,:,imode)+
     $                  rb(ibl)%qtion%qpf(:,:,:,imode)/zeff)
            ENDDO
            IF (ilayer==0) teff(:,:,:,1)=teff(:,:,:,1)+
     $          kboltz*(rb(ibl)%qtele_eq%qpf+
     $                  rb(ibl)%qtion_eq%qpf/zeff)
            CALL fft_nim('inverse',ng*mxb*myb,mpsq_block(ibl),lphi,
     $                   1_i4,teff,real_t,dealiase)
            DEALLOCATE(teff)
            ic=0
            yloopr2: DO iy=1,myb
              DO ix=1,mxb
                DO jg=1,ng
                  ic=ic+1
                  IF (ic<icst.OR.jg/=ig) CYCLE
                  IF (ic>icen) EXIT yloopr2
                  icq=(iy-1)*mxb+ix
                  icg=ic-icst+1
                  DO ip=1,nphi
                    nl_cfl=MAX(nl_cfl,kph(ix,iy)/md(ix,iy)*
     $                 ABS( SUM(real_b(:,icg,ip)**2)/mu0+
     $                      gamma* real_n(1,icg,ip)*real_t(1,icg,ip) -
     $                     (SUM(b0(:,ix,iy)**2)/mu0+
     $                      gamma*(rb(ibl)%qpres_eq%qpf(1,ig,icq)+
     $                             rb(ibl)%qpres_n0%qpf(1,ig,icq))) ) )
                  ENDDO
                ENDDO
              ENDDO
            ENDDO yloopr2
            DEALLOCATE(real_t)
          ENDIF
        ENDIF
        DEALLOCATE(b0,drz,md,kph,bigr)
      ENDDO
c-----------------------------------------------------------------------
c     tblocks:  the dimension vector is normal to each side,
c     intersecting the opposite vertex.
c-----------------------------------------------------------------------
      DO ibl=nrbl+1,nbl
        mc=tb(ibl)%mcell
        ng=tb(ibl)%ng
        ig=MAX(1_i4,ng/2)
        ALLOCATE(drz(2,3,mc,1),kph(mc,1),bigr(mc,1),nd_cell(mc,1))
        DO ic=1,mc
          iv=(/tb(ibl)%tgeom%vertex(ic,1),tb(ibl)%tgeom%vertex(ic,2),
     $         tb(ibl)%tgeom%vertex(ic,3),tb(ibl)%tgeom%vertex(ic,1)/)
          DO j=1,3
            drz(:,j,ic,1)=(/tb(ibl)%tgeom%ys(iv(j+1))
     $                     -tb(ibl)%tgeom%ys(iv(j  )),
     $                     -tb(ibl)%tgeom%xs(iv(j+1))
     $                     +tb(ibl)%tgeom%xs(iv(j  ))/)
            drz(:,j,ic,1)=2*tb(ibl)%tgeom%area(ic)*drz(:,j,ic,1)
     $                          /SUM(drz(:,j,ic,1)*drz(:,j,ic,1))
          ENDDO
          nd_cell(ic,1)=tb(ibl)%qnd_eq%qpf(1,ig,ic)
          IF (nonlinear.AND.
     $        (continuity=='n=0 only'.OR.continuity=='full'))
     $      nd_cell(ic,1)=nd_cell(ic,1)+tb(ibl)%qnd_n0%qpf(1,ig,ic)
          IF (geom=='tor') THEN
            bigr(ic,1)=
     $        (tb(ibl)%tgeom%xs(tb(ibl)%tgeom%vertex(ic,1))
     $        +tb(ibl)%tgeom%xs(tb(ibl)%tgeom%vertex(ic,2))
     $        +tb(ibl)%tgeom%xs(tb(ibl)%tgeom%vertex(ic,3)))/3
            kph(ic,1)=MAXVAL(ABS(keff_total))/bigr(ic,1)
          ELSE
            bigr(ic,1)=1
            kph(ic,1)=MAXVAL(ABS(keff_total))
          ENDIF
        ENDDO
c-----------------------------------------------------------------------
c       flow cfl in tblocks:
c-PRE
c-----------------------------------------------------------------------
        IF (nonlinear.OR.eq_flow/='none') THEN
          IF (nonlinear) THEN
            ALLOCATE(v(3,mc,1,nmodes),vx(4,mc,1,nmodes))
            DO imode=1,nmodes
              DO ic=1,mc
                v(:,ic,1,imode)=
     $           (ve(ibl)%arr(:,tb(ibl)%tgeom%vertex(ic,1),0,imode)
     $           +ve(ibl)%arr(:,tb(ibl)%tgeom%vertex(ic,2),0,imode)
     $           +ve(ibl)%arr(:,tb(ibl)%tgeom%vertex(ic,3),0,imode))/3
              ENDDO
              IF (eq_flow/='none'.AND.keff(imode)==0) THEN
                DO ic=1,mc
                  v(:,ic,1,imode)=v(:,ic,1,imode)
     $             +(ve_eq(ibl)%arr(:,tb(ibl)%tgeom%vertex(ic,1),0)
     $              +ve_eq(ibl)%arr(:,tb(ibl)%tgeom%vertex(ic,2),0)
     $              +ve_eq(ibl)%arr(:,tb(ibl)%tgeom%vertex(ic,3),0))/3
                ENDDO
              ENDIF
              DO j=1,3
                vx(j,:,:,imode)=
     $            pol_fac*SUM(v(1:2,:,:,imode)*drz(:,j,:,:),1)
     $                   /SUM(drz(:,j,:,:)**2,1)
              ENDDO
              vx(4,:,:,imode)=kph*v(3,:,:,imode)
            ENDDO
            DEALLOCATE(v)
            ALLOCATE(real_vx(4,mps_block(ibl),nphi))
            CALL fft_nim('inverse',mc,mps_block(ibl),lphi,
     $                   4_i4,vx,real_vx,dealiase)
            cfl=MAX(cfl,MAXVAL(ABS(real_vx)))
            DEALLOCATE(vx,real_vx)
          ELSE
            ALLOCATE(v0(3,mc,1),v0x(4,mc,1))
            DO ic=1,mc
              v0(:,ic,1)=
     $          +(ve_eq(ibl)%arr(:,tb(ibl)%tgeom%vertex(ic,1),0)
     $           +ve_eq(ibl)%arr(:,tb(ibl)%tgeom%vertex(ic,2),0)
     $           +ve_eq(ibl)%arr(:,tb(ibl)%tgeom%vertex(ic,3),0))/3
            ENDDO
            DO j=1,3
              v0x(j,:,:)=pol_fac*SUM(v0(1:2,:,:)*drz(:,j,:,:),1)
     $                          /SUM(drz(:,j,:,:)**2,1)
            ENDDO
            v0x(4,:,:)=kph*ABS(v0(3,:,:))
            cfl=MAX(cfl,MAXVAL(v0x))
            DEALLOCATE(v0,v0x)
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       electron flow cfl in tblocks:
c-----------------------------------------------------------------------
        IF (coefjve/=0._r8) THEN
          IF (nonlinear) THEN
            ALLOCATE(v(3,mc,1,nmodes),vx(4,mc,1,nmodes))
            DO imode=1,nmodes
              v(:,ic,1,imode)=1/3._r8*
     $          (ja(ibl)%arr(:,tb(ibl)%tgeom%vertex(ic,1),0,imode)
     $          +ja(ibl)%arr(:,tb(ibl)%tgeom%vertex(ic,2),0,imode)
     $          +ja(ibl)%arr(:,tb(ibl)%tgeom%vertex(ic,3),0,imode))
              IF (keff(imode)==0) THEN
                DO ic=1,mc
                  v(1:2,ic,1,imode)=v(1:2,ic,1,imode)+1/3._r8*
     $              (ja_eq(ibl)%arr(1:2,tb(ibl)%tgeom%vertex(ic,1),0)
     $              +ja_eq(ibl)%arr(1:2,tb(ibl)%tgeom%vertex(ic,2),0)
     $              +ja_eq(ibl)%arr(1:2,tb(ibl)%tgeom%vertex(ic,3),0))
                  v(3,ic,1,imode)=v(3,ic,1,imode)+bigr(ic,1)/3._r8*
     $              (ja_eq(ibl)%arr(3,tb(ibl)%tgeom%vertex(ic,1),0)
     $              +ja_eq(ibl)%arr(3,tb(ibl)%tgeom%vertex(ic,2),0)
     $              +ja_eq(ibl)%arr(3,tb(ibl)%tgeom%vertex(ic,3),0))
                ENDDO
              ENDIF
            ENDDO
            DO imode=1,nmodes
              DO ic=1,mc
                v(:,ic,1,imode)=coefjve*v(:,ic,1,imode)/nd_cell(ic,1)
     $          +(tb(ibl)%ve%fs(:,tb(ibl)%tgeom%vertex(ic,1),0,imode)
     $           +tb(ibl)%ve%fs(:,tb(ibl)%tgeom%vertex(ic,2),0,imode)
     $           +tb(ibl)%ve%fs(:,tb(ibl)%tgeom%vertex(ic,3),0,imode))/3
              ENDDO
              IF (eq_flow/='none'.AND.keff(imode)==0) THEN
                DO ic=1,mc
                  v(:,ic,1,imode)=v(:,ic,1,imode)
     $             +(ve_eq(ibl)%arr(:,tb(ibl)%tgeom%vertex(ic,1),0)
     $              +ve_eq(ibl)%arr(:,tb(ibl)%tgeom%vertex(ic,2),0)
     $              +ve_eq(ibl)%arr(:,tb(ibl)%tgeom%vertex(ic,3),0))/3
                ENDDO
              ENDIF
              DO j=1,3
                vx(j,:,:,imode)=
     $            pol_fac*SUM(v(1:2,:,:,imode)*drz(:,j,:,:),1)
     $                   /SUM(drz(:,j,:,:)**2,1)
              ENDDO
              vx(4,:,:,imode)=kph*v(3,:,:,imode)
            ENDDO
            DEALLOCATE(v)
            ALLOCATE(real_vx(4,mps_block(ibl),nphi))
            CALL fft_nim('inverse',mc,mps_block(ibl),lphi,
     $                   4_i4,vx,real_vx,dealiase)
            cfl=MAX(cfl,MAXVAL(ABS(real_vx)))
            DEALLOCATE(vx,real_vx)
          ELSE
            ALLOCATE(v0(3,mc,1),v0x(4,mc,1))
            DO ic=1,mc
              DO iq=1,2
                v0(iq,ic,1)=
     $            -(ja_eq(ibl)%arr(iq,tb(ibl)%tgeom%vertex(ic,1),0)
     $             +ja_eq(ibl)%arr(iq,tb(ibl)%tgeom%vertex(ic,2),0)
     $             +ja_eq(ibl)%arr(iq,tb(ibl)%tgeom%vertex(ic,3),0))
     $            /(3*elementary_q*nd_cell(ic,1))
              ENDDO
              v0(3,ic,1)=
     $          -(ja_eq(ibl)%arr(3,tb(ibl)%tgeom%vertex(ic,1),0)
     $           +ja_eq(ibl)%arr(3,tb(ibl)%tgeom%vertex(ic,2),0)
     $           +ja_eq(ibl)%arr(3,tb(ibl)%tgeom%vertex(ic,3),0))
     $          *bigr(ic,1)/(3*elementary_q*nd_cell(ic,1))
            ENDDO
            IF (eq_flow/='none') THEN
              DO ic=1,mc
                v0(:,ic,1)=v0(:,ic,1)
     $            +(ve_eq(ibl)%arr(:,tb(ibl)%tgeom%vertex(ic,1),0)
     $             +ve_eq(ibl)%arr(:,tb(ibl)%tgeom%vertex(ic,2),0)
     $             +ve_eq(ibl)%arr(:,tb(ibl)%tgeom%vertex(ic,3),0))/3
              ENDDO
            ENDIF
            DO j=1,3
              v0x(j,:,:)=pol_fac*SUM(v0(1:2,:,:)*drz(:,j,:,:),1)
     $                          /SUM(drz(:,j,:,:)**2,1)
            ENDDO
            v0x(4,:,:)=kph*ABS(v0(3,:,:))
            cfl=MAX(cfl,MAXVAL(v0x))
            DEALLOCATE(v0,v0x)
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       'linear' magneto-acoustic speed computations in tblocks:
c-----------------------------------------------------------------------
        kph=MAX(kph**2,pol_fac**2/MINVAL(SUM(drz**2,1),1))
        ALLOCATE(md(mc,1),b0(3,mc,1))
        DO ic=1,mc
          md(ic,1)=mtot*nd_cell(ic,1)
          bt=tb(ibl)%qbe_eq%qpf(:,ig,ic)
          IF (nonlinear) bt=bt+tb(ibl)%qbe_n0%qpf(:,ig,ic)
          b0(:,ic,1)=bt
          cl2t=SUM(bt**2)/mu0
          IF (beta>0) THEN
            cl2t=cl2t+gamma*tb(ibl)%qpres_eq%qpf(1,ig,ic)
            IF (nonlinear)
     $        cl2t=cl2t+gamma*tb(ibl)%qpres_n0%qpf(1,ig,ic)
          ENDIF
          lin_cfl=MAX(lin_cfl,cl2t*kph(ic,1)/md(ic,1))
        ENDDO
c-----------------------------------------------------------------------
c       'nonlinear' magneto-acoustic speed computations in tblocks:
c-----------------------------------------------------------------------
        IF (nonlinear) THEN
          icst=ipqst_block(ibl)
          icen=ipqen_block(ibl)
          real_b=>tb(ibl)%qbe_tot%qpf
          IF (beta<=0) THEN
            ic=0
            yloopt1: DO iy=1,mc
              DO jg=1,ng
                ic=ic+1
                IF (ic<icst.OR.jg/=ig) CYCLE
                IF (ic>icen) EXIT yloopt1
                icg=ic-icst+1
                DO ip=1,nphi
                  nl_cfl=MAX(nl_cfl,kph(iy,1)/(mu0*md(iy,1))*
     $              ABS(SUM(real_b(:,icg,ip)**2)-SUM(b0(:,iy,1)**2)))
                ENDDO
              ENDDO
            ENDDO yloopt1
          ELSE
            real_n=>tb(ibl)%qnd_tot%qpf
            ALLOCATE(teff(1,ng,mc,nmodes),
     $               real_t(1,mpsq_block(ibl),nphi))
            DO imode=1,nmodes
              teff(:,:,:,imode)=
     $          kboltz*(tb(ibl)%qtele%qpf(:,:,:,imode)+
     $                  tb(ibl)%qtion%qpf(:,:,:,imode)/zeff)
            ENDDO
            IF (ilayer==0) teff(:,:,:,1)=teff(:,:,:,1)+
     $          kboltz*(tb(ibl)%qtele_eq%qpf+
     $                  tb(ibl)%qtion_eq%qpf/zeff)
            CALL fft_nim('inverse',ng*mxb*myb,mpsq_block(ibl),lphi,
     $                   1_i4,teff,real_t,dealiase)
            DEALLOCATE(teff)
            ic=0
            yloopt2: DO iy=1,mc
              DO jg=1,ng
                ic=ic+1
                IF (ic<icst.OR.jg/=ig) CYCLE
                IF (ic>icen) EXIT yloopt2
                icg=ic-icst+1
                DO ip=1,nphi
                  nl_cfl=MAX(nl_cfl,kph(iy,1)/md(iy,1)*
     $               ABS( SUM(real_b(:,icg,ip)**2)/mu0+
     $                    gamma* real_n(1,icg,ip)*real_t(1,icg,ip) -
     $                   (SUM(b0(:,iy,1)**2)/mu0+
     $                    gamma*(tb(ibl)%qpres_eq%qpf(1,ig,iy)+
     $                           tb(ibl)%qpres_n0%qpf(1,ig,iy))) ) )
                ENDDO
              ENDDO
            ENDDO yloopt2
            DEALLOCATE(real_t)
          ENDIF
        ENDIF
        DEALLOCATE(b0,drz,kph,md,bigr,nd_cell)
      ENDDO
c-----------------------------------------------------------------------
c     compute the time step allowed by the flow cfl.
c-----------------------------------------------------------------------
      IF (nonlinear.OR.eq_flow/='none'.OR.coefjve/=0._r8) THEN
        IF (nprocs>1) THEN
          CALL mpi_allreduce(cfl,tmp,1,mpi_nim_real,mpi_max,
     $         mpi_comm_world,ierror)
          cfl=tmp
        ENDIF
        dtfl=v_cfl/(cfl+100.0*TINY(cfl))
      ELSE
        dtfl=dtm
      ENDIF
c-----------------------------------------------------------------------
c     time step allowed by nonlinear wave speed cfl.
c-----------------------------------------------------------------------
      lin_cfl=SQRT(lin_cfl)
      nl_cfl =SQRT(nl_cfl)
      IF (nprocs>1) THEN
        CALL mpi_allreduce(lin_cfl,tmp,1,mpi_nim_real,mpi_max,
     $       mpi_comm_world,ierror)
        lin_cfl=tmp
        CALL mpi_allreduce(nl_cfl,tmp,1,mpi_nim_real,mpi_max,
     $       mpi_comm_world,ierror)
        nl_cfl=tmp
      ENDIF
      IF (nonlinear.AND.nl_cfl>0) THEN
        dtnl=nl_cfl_lim/nl_cfl
      ELSE
        dtnl=dtm
      ENDIF
c-----------------------------------------------------------------------
c     apply limits to the time step.  dt may creap up if conditions
c     allow, but it is not permitted to sit at the flow limit to 
c     avoid an oscillating time step.
c-----------------------------------------------------------------------
      IF (.NOT.converged) dt=dt*dt_reduce
      IF (first_step) THEN
        IF (dt_initial/=0) THEN
          dt=MIN(dt_initial,dtm,dtfl,dtnl)
        ELSE
          dt=MIN(dtm,dtfl,dtnl)
        ENDIF
      ELSE
        tmp=MIN(dtm,dtfl,dtnl,dt_incr*dt)
        IF (tmp>dt*(1+dt_change_frac).OR.
     $     (tmp>dt.AND.n_same_dt>=n_dt_release.AND.tmp/=dtfl)) THEN
          dt=tmp
        ELSE IF (tmp<dt) THEN
          dt=MIN(tmp,dt*(1-dt_change_frac))
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     this prevents another dt drop on the next step.
c-----------------------------------------------------------------------
      IF (nonlinear.AND.dt==dtfl) dt=0.98_r8*dt
      IF (dt==dt_old) THEN
        n_same_dt=n_same_dt+1
      ELSE
        n_same_dt=0
      ENDIF
      IF (dt>dt_old) THEN
        n_dt_increase=0
      ELSE
        n_dt_increase=n_dt_increase+1
      ENDIF
      IF (nonlinear.OR.eq_flow/='none'.OR.coefjve/=0._r8) fl_cfl=dt*cfl
      IF (nonlinear) nl_cfl=dt*nl_cfl
      lin_cfl=dt*lin_cfl
c-----------------------------------------------------------------------
c     if this call precedes the time advance loop, do not change the
c     first_step flag.
c-----------------------------------------------------------------------
      IF (istep==nstop-nstep) RETURN
      first_step=.false.
c-----------------------------------------------------------------------
c     stop if the time step is too small.
c-----------------------------------------------------------------------
      IF (dt<=dt_stop*dtm) THEN
        CALL nim_output(.false.)
        WRITE(msg,'(a,i5,a)')
     $    'Newimpdt: dt <= dt_stop*dtm at cycle number ', istep,'.'
        CALL nim_stop(TRIM(msg))
      ENDIF
c-----------------------------------------------------------------------
c     this a good place to zero out iteration counters.
c-----------------------------------------------------------------------
      mhdits=0
      hallits=0
      jaits=0
      vmhdits=0
      viscits=0
      tiits=0
      teits=0
      dbits=0
      neoits=0
      ndits=0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE new_dt
c-----------------------------------------------------------------------
c     subprogram 6. add_mass.
c     add the mass matrix to the lmat operator.  lumping is no longer
c     allowed.
c-----------------------------------------------------------------------
      SUBROUTINE add_mass(mat,nqty)
      USE local
      USE fields
      USE matrix_storage_mod
      IMPLICIT NONE

      TYPE(global_matrix_type), INTENT(INOUT) :: mat
      INTEGER(i4), INTENT(IN) :: nqty

      INTEGER(i4) :: ibl,iq,jq,iqt,jqt,iqm,jqm,iv,mv,ityp,jtyp
      REAL(r8), DIMENSION(:,:,:,:,:,:), POINTER :: rmat,mass
      TYPE(matrix_element_type3), DIMENSION(:), POINTER :: tmat,tmass
c-----------------------------------------------------------------------
c     loop over rblocks first.  any regularity conditions are set
c     after this operation, so mass_mat(1) is always used.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        DO ityp=1,mat%rbl_mat(ibl)%nbtype
          DO jtyp=1,mat%rbl_mat(ibl)%nbtype
            rmat=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
            mass=>mass_mat(1)%rbl_mat(ibl)%mat(jtyp,ityp)%arr
            DO iq=1,nqty
              DO iqm=1,mass_mat(1)%rbl_mat(ibl)%nq_type(ityp)
                iqt=nqty*(iqm-1)+iq
                DO jqm=1,mass_mat(1)%rbl_mat(ibl)%nq_type(jtyp)
                  jqt=nqty*(jqm-1)+iq
                  rmat(jqt,:,:,iqt,:,:)=rmat(jqt,:,:,iqt,:,:)
     $                                 +mass(jqm,:,:,iqm,:,:)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     now tblocks.
c-----------------------------------------------------------------------
      DO ibl=nrbl+1,nbl
        tmat=>mat%tbl_mat(ibl)%lmat
        tmass=>mass_mat(1)%tbl_mat(ibl)%lmat
        mv=tb(ibl)%mvert
        DO iv=0,mv
          DO iq=1,nqty
            tmat(iv)%element(iq,iq,:)
     $          =tmat (iv)%element(iq,iq,:)
     $          +tmass(iv)%element(1 ,1 ,:)
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE add_mass
c-----------------------------------------------------------------------
c     subprogram 7. add_mass_comp.
c     add the mass matrix to the complex lmat operator.
c-----------------------------------------------------------------------
      SUBROUTINE add_mass_comp(mat,nqty)
      USE local
      USE fields
      USE matrix_storage_mod
      IMPLICIT NONE

      TYPE(complex_matrix_type), INTENT(INOUT) :: mat
      INTEGER(i4), INTENT(IN) :: nqty

      INTEGER(i4) :: ibl,iq,jq,iqt,jqt,iqm,jqm,iv,mv,ityp,jtyp
      INTEGER(i4), DIMENSION(nqty) :: mass_type
      COMPLEX(r8), DIMENSION(:,:,:,:,:,:), POINTER :: rmat
      REAL(r8), DIMENSION(:,:,:,:,:,:), POINTER :: mass
      TYPE(comp_matrix_element_type3), DIMENSION(:), POINTER :: tmat
      TYPE(matrix_element_type3), DIMENSION(:), POINTER :: tmass
c-----------------------------------------------------------------------
c     loop over rblocks first.  any regularity conditions are set
c     after this operation, so mass_mat(1) is always used.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        DO ityp=1,SIZE(mat%rbl_mat(ibl)%mat,2)
          DO jtyp=1,SIZE(mat%rbl_mat(ibl)%mat,1)
            rmat=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
            mass=>mass_mat(1)%rbl_mat(ibl)%mat(jtyp,ityp)%arr
            DO iq=1,nqty
              DO iqm=1,mass_mat(1)%rbl_mat(ibl)%nq_type(ityp)
                iqt=nqty*(iqm-1)+iq
                DO jqm=1,mass_mat(1)%rbl_mat(ibl)%nq_type(jtyp)
                  jqt=nqty*(jqm-1)+iq
                  rmat(jqt,:,:,iqt,:,:)=rmat(jqt,:,:,iqt,:,:)
     $                                 +mass(jqm,:,:,iqm,:,:)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     now tblocks.
c-----------------------------------------------------------------------
      DO ibl=nrbl+1,nbl
        tmat=>mat%tbl_mat(ibl)%lmat
        tmass=>mass_mat(1)%tbl_mat(ibl)%lmat
        mv=tb(ibl)%mvert
        DO iv=0,mv
          DO iq=1,nqty
            tmat(iv)%element(iq,iq,:)
     $          =tmat (iv)%element(iq,iq,:)
     $          +tmass(iv)%element(1 ,1 ,:)
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE add_mass_comp
c-----------------------------------------------------------------------
c     subprogram 8. iter_out.
c     write iterative solver information to an ASCII file.
c-----------------------------------------------------------------------
      SUBROUTINE iter_out(eqn,seed,iterations,error)
      USE local
      USE global
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: eqn,seed
      INTEGER(i4), INTENT(IN) :: iterations
      REAL(r8), INTENT(IN) :: error

      WRITE(it_unit,'(a,i5,a,a9,a,a8,a,i3,a,1pe13.5)')
     $  'Cycle: ',istep,
     $  ' Eqn: ',eqn,
     $  ' Seed: ',seed,
     $  ' Iterations: ',iterations,
     $  ' Error: ',error
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_out
c-----------------------------------------------------------------------
c     subprogram 9. loop_voltage
c     find the loop voltage from a specified function of time, or
c     to relax the current to a desired value. 
c-----------------------------------------------------------------------
      SUBROUTINE loop_voltage
      USE local
      USE input
      USE global
      IMPLICIT NONE

      LOGICAL, SAVE :: initialized=.false.
c-----------------------------------------------------------------------
c     set old voltage:
c-----------------------------------------------------------------------
      IF (.NOT.initialized) THEN
        initialized=.true.
        volt_old=loop_volt
        IF (i_desired/=0) THEN
          volt=loop_volt
          RETURN
        ENDIF
      ELSE
        volt_old=volt
      ENDIF
c-----------------------------------------------------------------------
c     two options:
c-----------------------------------------------------------------------
      IF (i_desired/=0) THEN
        volt=volt_old+dt*loop_rate*(i_desired-i_n0)
        IF (i_n0_old/=0.AND.dt_old/=0)
     $    volt=volt-dt*loop_rate2*(i_n0-i_n0_old)/dt_old
      ELSE
        IF (tloopv1>=tloopv0) THEN
          IF (t+dt>=tloopv1) THEN
            volt=loop_volt
          ELSE IF (t+dt<=tloopv0) THEN
            volt=0
          ELSE
            volt=loop_volt*(t+dt-tloopv0)/MAX(tloopv1-tloopv0,smallnum)
          ENDIF
        ELSE
          IF (t+dt>=tloopv0) THEN
            volt=0
          ELSE IF (t+dt<=tloopv1) THEN
            volt=loop_volt
          ELSE
            volt=loop_volt
     $          *(1-(t+dt-tloopv1)/MAX(tloopv0-tloopv1,smallnum))
          ENDIF
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE loop_voltage
c-----------------------------------------------------------------------
c     subprogram 10. vertical_efield
c     set an applied electric field for the present time.
c-----------------------------------------------------------------------
      SUBROUTINE vertical_efield
      USE local
      USE input
      USE global
      IMPLICIT NONE

      LOGICAL, SAVE :: initialized=.false.
c-----------------------------------------------------------------------
c     set old voltage:
c-----------------------------------------------------------------------
      e_vert_old=e_vert
c-----------------------------------------------------------------------
c     set e_vert according to t.
c-----------------------------------------------------------------------
      IF (t_e_vert1>=t_e_vert0) THEN
        IF (t+dt>=t_e_vert1) THEN
          e_vert=e_vertical
        ELSE IF (t+dt<=t_e_vert0) THEN
          e_vert=0
        ELSE
          e_vert=e_vertical*(t+dt-t_e_vert0)
     $          /MAX(t_e_vert1-t_e_vert0,smallnum)
        ENDIF
      ELSE
        IF (t+dt>=t_e_vert0) THEN
          e_vert=0
        ELSE IF (t+dt<=t_e_vert1) THEN
          e_vert=e_vertical
        ELSE
          e_vert=e_vertical
     $          *(1-(t+dt-t_e_vert1)/MAX(t_e_vert0-t_e_vert1,smallnum))
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     set old value to the new if this is called during initialization.
c-----------------------------------------------------------------------
      IF (.NOT.initialized) THEN
        initialized=.true.
        e_vert_old=e_vert
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vertical_efield
c-----------------------------------------------------------------------
c     subprogram 11. qp_fft_save.
c     save real data as a function of toroidal angle at quadrature
c     points to reduce the number of ffts called during 3D matrix
c     iterations.
c-----------------------------------------------------------------------
      SUBROUTINE qp_fft_save(fcmp,fphi,lx,ly,mps,nq,ng,feq)
      USE local
      USE global
      USE input
      USE fft_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: lx,ly,mps,nq,ng
      REAL(r8), DIMENSION(nq,mps,nphi), INTENT(INOUT) :: fphi
      REAL(r8), DIMENSION(nq,ng,lx*ly), INTENT(IN) :: feq
      COMPLEX(r8), DIMENSION(nq,ng,lx*ly,nmodes), INTENT(IN) :: fcmp

      INTEGER(i4) :: ig,im,iq,ix,iy,nqm,ipol,iphi
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: fcmp2
c-----------------------------------------------------------------------
c     reordering the quadrature-point index eliminates the data shuffle
c     that was here, and mps includes a combination of poloidal indices
c     and quadrature-point indices.
c-----------------------------------------------------------------------
      ALLOCATE(fcmp2(nq,ng,lx*ly,nmodes))
      DO im=1,nmodes
        IF (keff(im)==0) THEN
          fcmp2(:,:,:,im)=fcmp(:,:,:,im)+feq
        ELSE
          fcmp2(:,:,:,im)=fcmp(:,:,:,im)
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     perform ffts of the packed data.
c-----------------------------------------------------------------------
      CALL fft_nim('inverse',ng*lx*ly,mps,lphi,nq,fcmp2,fphi,dealiase)
      DEALLOCATE(fcmp2)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE qp_fft_save
c-----------------------------------------------------------------------
c     subprogram 12. qp_fft_noeq_save.
c     save real data as a function of toroidal angle at quadrature
c     points to reduce the number of ffts called during 3D matrix
c     iterations.  this version does not have an 'equilibrium'
c     component to be added.
c-----------------------------------------------------------------------
      SUBROUTINE qp_fft_noeq_save(fcmp,fphi,lx,ly,mps,nq,ng)
      USE local
      USE global
      USE input
      USE fft_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: lx,ly,mps,nq,ng
      REAL(r8), DIMENSION(nq,mps,nphi), INTENT(INOUT) :: fphi
      COMPLEX(r8), DIMENSION(nq,ng,lx*ly,nmodes), INTENT(INOUT) :: fcmp
c-----------------------------------------------------------------------
c     reordering the quadrature-point index eliminates the data shuffle
c     that was here, and mps includes a combination of poloidal indices
c     and quadrature-point indices.  this routine now calls the fft
c     routine directly.
c-----------------------------------------------------------------------
      CALL fft_nim('inverse',ng*lx*ly,mps,lphi,nq,fcmp,fphi,dealiase)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE qp_fft_noeq_save
c-----------------------------------------------------------------------
c     subprogram 13. soln_save.
c     save and restore old solutions in case the iterative solver does
c     not converge.
c-----------------------------------------------------------------------
      SUBROUTINE soln_save(converged)
      USE local
      USE input
      USE fields
      USE global
      USE vector_type_mod
      USE extrap_mod
      USE pardata
      USE computation_pointers
      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: converged

      LOGICAL :: first_call=.true.
      INTEGER(i4) :: ibl,mxb,myb
c-----------------------------------------------------------------------
c     create the space at the first call.
c-----------------------------------------------------------------------
      IF (first_call) THEN
        ALLOCATE(be_old(nbl),ja_old(nbl),ve_old(nbl),p_old(nbl),
     $           pe_old(nbl),nd_old(nbl),te_old(nbl),ti_old(nbl))
        DO ibl=1,nbl
          IF (ibl<=nrbl) THEN
            mxb=rb(ibl)%mx
            myb=rb(ibl)%my
          ELSE
            mxb=tb(ibl)%mvert
            myb=0
          ENDIF
          CALL vector_type_alloc(be_old(ibl),poly_degree,mxb,myb,
     $                           3_i4,nmodes)
          CALL vector_type_alloc(ja_old(ibl),poly_degree,mxb,myb,
     $                           3_i4,nmodes)
          CALL vector_type_alloc(ve_old(ibl),poly_degree,mxb,myb,
     $                           3_i4,nmodes)
          CALL vector_type_alloc(p_old(ibl),poly_degree,mxb,myb,
     $                           1_i4,nmodes)
          CALL vector_type_alloc(pe_old(ibl),poly_degree,mxb,myb,
     $                           1_i4,nmodes)
          CALL vector_type_alloc(nd_old(ibl),poly_degree,mxb,myb,
     $                           1_i4,nmodes)
          CALL vector_type_alloc(te_old(ibl),poly_degree,mxb,myb,
     $                           1_i4,nmodes)
          CALL vector_type_alloc(ti_old(ibl),poly_degree,mxb,myb,
     $                           1_i4,nmodes)
        ENDDO
        first_call=.false.
      ENDIF
c-----------------------------------------------------------------------
c     copy solution if solves have converged.  if not, restore the
c     old solution and reduce dt.
c-----------------------------------------------------------------------
      IF (converged) THEN
        DO ibl=1,nbl
          be_old(ibl)=be(ibl)
          ja_old(ibl)=ja(ibl)
          ve_old(ibl)=ve(ibl)
          p_old(ibl)=pres(ibl)
          pe_old(ibl)=prese(ibl)
          nd_old(ibl)=nd(ibl)
          te_old(ibl)=tele(ibl)
          ti_old(ibl)=tion(ibl)
        ENDDO
c-----------------------------------------------------------------------
c     if not, restore the old solution and reduce dt.
c     also, deallocate any temporary space for algebraic operations
c     that is left from the management routine where the solver
c     failed.
c-----------------------------------------------------------------------
      ELSE
        IF (node == 0) THEN
          WRITE(nim_wr,'(/a,i5)')
     $     ' Soln_save is restoring previous data due to lack of'//
     $     ' convergence at cycle ',istep
          IF (.NOT.out_opened) THEN
            OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $           POSITION='APPEND')
            out_opened=.true.
          ENDIF
          WRITE(out_unit,'(/a,i5)')
     $     ' Soln_save is restoring previous data due to lack of'//
     $     ' convergence at cycle ',istep
        ENDIF
        n_dt_increase=-1
        DO ibl=1,nbl
          CALL vector_type_dealloc(crhs(ibl))
          CALL vector_type_dealloc(cvecn(ibl))
          CALL vector_type_dealloc(sln(ibl))
          CALL vector_type_dealloc(csln(ibl))
          CALL vector_type_dealloc(vectr(ibl))
          CALL vector_type_dealloc(cvectr(ibl))
          be(ibl)=be_old(ibl)
          ja(ibl)=ja_old(ibl)
          ve(ibl)=ve_old(ibl)
          pres(ibl)=p_old(ibl)
          prese(ibl)=pe_old(ibl)
          nd(ibl)=nd_old(ibl)
          tele(ibl)=te_old(ibl)
          tion(ibl)=ti_old(ibl)
        ENDDO
        CALL extrap_init_array
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE soln_save
c-----------------------------------------------------------------------
c     subprogram 14. zero_last_comp.
c     this subroutine is used to zero-out the last Fourier component in
c     nonlinear computations without standard dealiasing.  these
c     components are necessarily incomplete, anyway, in that they must
c     be real for the ffts.
c-----------------------------------------------------------------------
      SUBROUTINE zero_last_comp(nb,compv)
      USE local
      USE vector_type_mod
      USE global
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nb
      TYPE(cvector_type), INTENT(INOUT), DIMENSION(nb) :: compv

      INTEGER(i4) :: ib,im
c-----------------------------------------------------------------------
c     loop over blocks.
c-----------------------------------------------------------------------
      DO im=1,nmodes
        IF (nindex(im)==nphi/2) THEN
          DO ib=1,nb
            IF (ASSOCIATED(compv(ib)%arr)) compv(ib)%arr(:,:,:,im)=0
            IF (ASSOCIATED(compv(ib)%arrh)) compv(ib)%arrh(:,:,:,:,im)=0
            IF (ASSOCIATED(compv(ib)%arrv)) compv(ib)%arrv(:,:,:,:,im)=0
            IF (ASSOCIATED(compv(ib)%arri)) compv(ib)%arri(:,:,:,:,im)=0
          ENDDO
          EXIT
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE zero_last_comp
c-----------------------------------------------------------------------
c     subprogram 15. fourier_damp.
c     this subroutine is used to damp the largest coefficients of
c     a Fourier expansion.
c-----------------------------------------------------------------------
      SUBROUTINE fourier_damp(nb,dtstep,compv)
      USE local
      USE vector_type_mod
      USE global
      USE input
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nb
      REAL(r8), INTENT(IN) :: dtstep
      TYPE(cvector_type), INTENT(INOUT), DIMENSION(nb) :: compv

      INTEGER(i4) :: ib,im,nmax
      REAL(r8) :: dfac
c-----------------------------------------------------------------------
c     loop over blocks.
c-----------------------------------------------------------------------
      nmax=nindex_total(nmodes_total)
      DO im=1,nmodes
        IF (nindex(im)<=nmax-nfdamp) CYCLE
        dfac=MAX(0._r8,1._r8-fmx_drate*dtstep*
     $                 REAL(nindex(im)-nmax+nfdamp,r8)/nfdamp)
        DO ib=1,nb
          IF (ASSOCIATED(compv(ib)%arr))
     $      compv(ib)%arr(:,:,:,im)=dfac*compv(ib)%arr(:,:,:,im)
          IF (ASSOCIATED(compv(ib)%arrh))
     $      compv(ib)%arrh(:,:,:,:,im)=dfac*compv(ib)%arrh(:,:,:,:,im)
          IF (ASSOCIATED(compv(ib)%arrv))
     $      compv(ib)%arrv(:,:,:,:,im)=dfac*compv(ib)%arrv(:,:,:,:,im)
          IF (ASSOCIATED(compv(ib)%arri))
     $      compv(ib)%arri(:,:,:,:,im)=dfac*compv(ib)%arri(:,:,:,:,im)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fourier_damp
c-----------------------------------------------------------------------
c     subprogram 16. set_nodal_min.
c     this subroutine is used to enforce a minimum value on a field at
c     the nodes of the spectral-element/Fourier expansion.  the
c     input for this routine is the data structures for the solution
c     field and for the corresponding equilibrium field, the number of
c     blocks, and the desired minimum value.  only the solution field
c     is altered.
c-----------------------------------------------------------------------
      SUBROUTINE set_nodal_min(nb,minv,compv,realv)
      USE local
      USE vector_type_mod
      USE global
      USE input
      USE pardata
      USE fft_mod
      USE mpi_nim
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nb
      REAL(r8), INTENT(IN) :: minv
      TYPE(cvector_type), DIMENSION(:), POINTER :: compv
      TYPE(vector_type), DIMENSION(:), POINTER :: realv

      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: carr
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: rarr

      INTEGER(i4) :: ib,ipol,npol,npolr,rep,im,ix,iy,is,ns,nq,mxb,myb,
     $               ierror
      LOGICAL, PARAMETER :: check_min=.false.
      REAL(r8) :: global_min,tmp,global_min0
c-----------------------------------------------------------------------
c     determine the amount of data to transform.
c-----------------------------------------------------------------------
      global_min=HUGE(global_min)
      global_min0=HUGE(global_min0)
      DO ib=1,nb
        npol=0
        IF (ASSOCIATED(compv(ib)%arr)) THEN
          nq=SIZE(compv(ib)%arr,1)
          mxb=SIZE(compv(ib)%arr,2)-1
          myb=SIZE(compv(ib)%arr,3)-1
          npol=npol+(mxb+1)*(myb+1)
        ENDIF
        IF (ASSOCIATED(compv(ib)%arrh)) THEN
          nq=SIZE(compv(ib)%arrh,1)
          mxb=SIZE(compv(ib)%arrh,3)
          myb=SIZE(compv(ib)%arrh,4)-1
          npol=npol+SIZE(compv(ib)%arrh,2)*mxb*(myb+1)
        ENDIF
        IF (ASSOCIATED(compv(ib)%arrv)) THEN
          nq=SIZE(compv(ib)%arrv,1)
          mxb=SIZE(compv(ib)%arrv,3)-1
          myb=SIZE(compv(ib)%arrv,4)
          npol=npol+SIZE(compv(ib)%arrv,2)*(mxb+1)*myb
        ENDIF
        IF (ASSOCIATED(compv(ib)%arri)) THEN
          nq=SIZE(compv(ib)%arri,1)
          mxb=SIZE(compv(ib)%arri,3)
          myb=SIZE(compv(ib)%arri,4)
          npol=npol+SIZE(compv(ib)%arri,2)*mxb*myb
        ENDIF
c-----------------------------------------------------------------------
c       allocate and fill arrays for the transform.
c-----------------------------------------------------------------------
        npolr=npol/nlayers
        rep=MODULO(npol,nlayers)
        npolr=npolr+MIN(rep,ilayer+1)-MIN(rep,ilayer)
        ALLOCATE(carr(nq,npol,nmodes))
        ALLOCATE(rarr(nq,npolr,nphi))
        DO im=1,nmodes
          ipol=1
          IF (ASSOCIATED(compv(ib)%arr)) THEN
            IF (keff(im)==0._r8)
     $        compv(ib)%arr(:,:,:,im)=compv(ib)%arr(:,:,:,im)+
     $                                realv(ib)%arr
            DO iy=0,myb
              DO ix=0,mxb
                carr(1:nq,ipol,im)=compv(ib)%arr(1:nq,ix,iy,im)
                ipol=ipol+1
              ENDDO
            ENDDO
          ENDIF
          IF (ASSOCIATED(compv(ib)%arrh)) THEN
            ns=SIZE(compv(ib)%arrh,2)
            IF (keff(im)==0._r8)
     $        compv(ib)%arrh(:,:,:,:,im)=compv(ib)%arrh(:,:,:,:,im)+
     $                                   realv(ib)%arrh
            DO iy=0,myb
              DO ix=1,mxb
                DO is=1,ns
                  carr(1:nq,ipol,im)=compv(ib)%arrh(1:nq,is,ix,iy,im)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          IF (ASSOCIATED(compv(ib)%arrv)) THEN
            ns=SIZE(compv(ib)%arrv,2)
            IF (keff(im)==0._r8)
     $        compv(ib)%arrv(:,:,:,:,im)=compv(ib)%arrv(:,:,:,:,im)+
     $                                   realv(ib)%arrv
            DO iy=1,myb
              DO ix=0,mxb
                DO is=1,ns
                  carr(1:nq,ipol,im)=compv(ib)%arrv(1:nq,is,ix,iy,im)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          IF (ASSOCIATED(compv(ib)%arri)) THEN
            ns=SIZE(compv(ib)%arri,2)
            IF (keff(im)==0._r8)
     $        compv(ib)%arri(:,:,:,:,im)=compv(ib)%arri(:,:,:,:,im)+
     $                                   realv(ib)%arri
            DO iy=1,myb
              DO ix=1,mxb
                DO is=1,ns
                  carr(1:nq,ipol,im)=compv(ib)%arri(1:nq,is,ix,iy,im)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
c-----------------------------------------------------------------------
c       transform to the nodes of the Fourier expansion, apply the
c       lower bound, then transform back to Fourier coefficients.
c-----------------------------------------------------------------------
        CALL fft_nim('inverse',npol,npolr,lphi,nq,carr,rarr,dealiase)
        IF (check_min) THEN
          CALL fft_nim('inverse',npol,npolr,lphi,nq,carr,rarr,dealiase)
          global_min0=MIN(global_min0,MINVAL(rarr))
        ENDIF
        rarr=MAX(rarr,minv)
        CALL fft_nim('forward',npol,npolr,lphi,nq,carr,rarr,dealiase)
c-----------------------------------------------------------------------
c       with dealiasing, values set at nodes are altered.  this
c       check reports the minimum value that results after the
c       forward transform.  it should only be used for test cases.
c-----------------------------------------------------------------------
        IF (check_min) THEN
          CALL fft_nim('inverse',npol,npolr,lphi,nq,carr,rarr,dealiase)
          global_min=MIN(global_min,MINVAL(rarr))
        ENDIF
c-----------------------------------------------------------------------
c       put the adjusted data into the solution arrays, then subtract
c       the steady-state part from the n=0 component.
c-----------------------------------------------------------------------
        DO im=1,nmodes
          ipol=1
          IF (ASSOCIATED(compv(ib)%arr)) THEN
            DO iy=0,myb
              DO ix=0,mxb
                compv(ib)%arr(1:nq,ix,iy,im)=carr(1:nq,ipol,im)
                ipol=ipol+1
              ENDDO
            ENDDO
            IF (keff(im)==0._r8)
     $        compv(ib)%arr(:,:,:,im)=compv(ib)%arr(:,:,:,im)-
     $                                realv(ib)%arr
          ENDIF
          IF (ASSOCIATED(compv(ib)%arrh)) THEN
            ns=SIZE(compv(ib)%arrh,2)
            DO iy=0,myb
              DO ix=1,mxb
                DO is=1,ns
                  compv(ib)%arrh(1:nq,is,ix,iy,im)=carr(1:nq,ipol,im)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
            IF (keff(im)==0._r8)
     $        compv(ib)%arrh(:,:,:,:,im)=compv(ib)%arrh(:,:,:,:,im)-
     $                                   realv(ib)%arrh
          ENDIF
          IF (ASSOCIATED(compv(ib)%arrv)) THEN
            ns=SIZE(compv(ib)%arrv,2)
            DO iy=1,myb
              DO ix=0,mxb
                DO is=1,ns
                  compv(ib)%arrv(1:nq,is,ix,iy,im)=carr(1:nq,ipol,im)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
            IF (keff(im)==0._r8)
     $        compv(ib)%arrv(:,:,:,:,im)=compv(ib)%arrv(:,:,:,:,im)-
     $                                   realv(ib)%arrv
          ENDIF
          IF (ASSOCIATED(compv(ib)%arri)) THEN
            ns=SIZE(compv(ib)%arri,2)
            DO iy=1,myb
              DO ix=1,mxb
                DO is=1,ns
                  compv(ib)%arri(1:nq,is,ix,iy,im)=carr(1:nq,ipol,im)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
            IF (keff(im)==0._r8)
     $        compv(ib)%arri(:,:,:,:,im)=compv(ib)%arri(:,:,:,:,im)-
     $                                   realv(ib)%arri
          ENDIF
        ENDDO
        DEALLOCATE(carr,rarr)
      ENDDO
c-----------------------------------------------------------------------
c     complete the global_min computation and output.
c-----------------------------------------------------------------------
      IF (check_min) THEN
        IF (nprocs>1) THEN
          CALL mpi_allreduce(global_min,tmp,1,mpi_nim_real,mpi_min,
     $         mpi_comm_world,ierror)
          global_min=tmp
          CALL mpi_allreduce(global_min0,tmp,1,mpi_nim_real,mpi_min,
     $         mpi_comm_world,ierror)
          global_min0=tmp
        ENDIF
        IF (node==0) THEN
          WRITE(nim_wr,'(a,es12.4)') " Min before: ",global_min0
          WRITE(nim_wr,'(a,es12.4,a,es12.4)') " Min specified: ",minv,
     $      " Min after transform: ",global_min
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE set_nodal_min
