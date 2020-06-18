c routines for timing statistics
c storage for accumulating incremental times

      MODULE time
      USE local
      IMPLICIT NONE

c cummulative times and counters

      INTEGER(i4) :: seamcount = 0
      INTEGER(i4) :: segcount = 0
      INTEGER(i4) :: seamcount_loop,segcount_loop

      REAL(r8) :: time_io = 0.0
      REAL(r8) :: time_seam = 0.0
      REAL(r8) :: time_seg = 0.0
      REAL(r8) :: time_iter = 0.0
      REAL(r8) :: time_fac = 0.0
      REAL(r8) :: time_fft = 0.0
      REAL(r8) :: time_mat = 0.0
      REAL(r8) :: time_rhs = 0.0
      REAL(r8) :: time_line = 0.0
      REAL(r8) :: time_stcon = 0.0
      REAL(r8) :: time_exsolv = 0.0
      REAL(r8) :: time_exfact = 0.0
      REAL(r8) :: time_matvec = 0.0
      REAL(r8) :: time_qpfld = 0.0
      REAL(r8) :: time_io_loop,time_seam_loop,time_seg_loop
      REAL(r8) :: time_iter_loop,time_fft_loop,time_mat_loop
      REAL(r8) :: time_fac_loop,time_line_loop,time_stcon_loop
      REAL(r8) :: time_rhs_loop,time_exsolv_loop,time_exfact_loop
      REAL(r8) :: time_matvec_loop,time_qpfld_loop

      REAL(r8) :: time_total_start,time_total_end
      REAL(r8) :: time_loop_start,time_loop_end
      REAL(r8) :: timestart,timeend

      CONTAINS

c initialize cummulative values

      SUBROUTINE timer_init
      USE local
      IMPLICIT NONE

      seamcount = 0
      segcount = 0
      time_io = 0.0
      time_seam = 0.0
      time_seg = 0.0
      time_iter = 0.0
      time_fac = 0.0
      time_fft = 0.0
      time_mat = 0.0
      time_rhs = 0.0
      time_line = 0.0
      time_stcon = 0.0
      time_exsolv = 0.0
      time_exfact = 0.0
      time_matvec = 0.0
      time_qpfld = 0.0

      RETURN
      END SUBROUTINE timer_init

c store io & seam times/counters from within main timestepping loop

      SUBROUTINE timer_close
      USE local
      IMPLICIT NONE

      seamcount_loop = seamcount
      segcount_loop = segcount
      time_io_loop = time_io
      time_seam_loop = time_seam
      time_seg_loop = time_seg
      time_iter_loop = time_iter
      time_fac_loop = time_fac
      time_fft_loop = time_fft
      time_mat_loop = time_mat
      time_rhs_loop = time_rhs
      time_line_loop = time_line
      time_stcon_loop = time_stcon
      time_exsolv_loop = time_exsolv
      time_exfact_loop = time_exfact
      time_matvec_loop = time_matvec
      time_matvec_loop = time_matvec
      time_qpfld_loop = time_qpfld

      RETURN
      END SUBROUTINE timer_close

c print-out stats, averaged across procs

      SUBROUTINE timer_stats
      USE local
      USE mpi_nim
      USE pardata
      IMPLICIT NONE
      REAL(r8) :: tmp,time_setup
      REAL(r8) :: time_total,time_loop
      INTEGER(i4) :: ierror,uwrite,iwr

c average cummulative timers across all procs

      IF (nprocs>1) THEN
        CALL mpi_allreduce(time_io_loop,tmp,1,mpi_nim_real,mpi_sum,
     $       mpi_comm_world,ierror)
        time_io_loop = tmp/nprocs
        CALL mpi_allreduce(time_seam_loop,tmp,1,mpi_nim_real,mpi_sum,
     $       mpi_comm_world,ierror)
        time_seam_loop = tmp/nprocs
        CALL mpi_allreduce(time_seg_loop,tmp,1,mpi_nim_real,mpi_sum,
     $       mpi_comm_world,ierror)
        time_seg_loop = tmp/nprocs
        CALL mpi_allreduce(time_iter_loop,tmp,1,mpi_nim_real,mpi_sum,
     $       mpi_comm_world,ierror)
        time_iter_loop = tmp/nprocs
        CALL mpi_allreduce(time_fac_loop,tmp,1,mpi_nim_real,mpi_sum,
     $       mpi_comm_world,ierror)
        time_fac_loop = tmp/nprocs
        CALL mpi_allreduce(time_fft_loop,tmp,1,mpi_nim_real,mpi_sum,
     $       mpi_comm_world,ierror)
        time_fft_loop = tmp/nprocs
        CALL mpi_allreduce(time_mat_loop,tmp,1,mpi_nim_real,mpi_sum,
     $       mpi_comm_world,ierror)
        time_mat_loop = tmp/nprocs
        CALL mpi_allreduce(time_rhs_loop,tmp,1,mpi_nim_real,mpi_sum,
     $       mpi_comm_world,ierror)
        time_rhs_loop = tmp/nprocs
        CALL mpi_allreduce(time_line_loop,tmp,1,mpi_nim_real,mpi_sum,
     $       mpi_comm_world,ierror)
        time_line_loop = tmp/nprocs
        CALL mpi_allreduce(time_stcon_loop,tmp,1,mpi_nim_real,mpi_sum,
     $       mpi_comm_world,ierror)
        time_stcon_loop = tmp/nprocs
        CALL mpi_allreduce(time_exsolv_loop,tmp,1,mpi_nim_real,mpi_sum,
     $       mpi_comm_world,ierror)
        time_exsolv_loop = tmp/nprocs
        CALL mpi_allreduce(time_exfact_loop,tmp,1,mpi_nim_real,mpi_sum,
     $       mpi_comm_world,ierror)
        time_exfact_loop = tmp/nprocs
        CALL mpi_allreduce(time_matvec_loop,tmp,1,mpi_nim_real,mpi_sum,
     $       mpi_comm_world,ierror)
        time_matvec_loop = tmp/nprocs
        CALL mpi_allreduce(time_qpfld_loop,tmp,1,mpi_nim_real,mpi_sum,
     $       mpi_comm_world,ierror)
        time_qpfld_loop = tmp/nprocs
      ENDIF

c compute derivative times

      time_total = MAX(time_total_end - time_total_start,
     $             SQRT(TINY(time_total)))
      time_loop = MAX(time_loop_end - time_loop_start,
     $             SQRT(TINY(time_total)))
      time_setup = time_total - time_loop

c print-out

   10 format(a,2es12.5)

      if (node == 0) then

        IF (.NOT.out_opened) THEN
          OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $         POSITION='APPEND')
          out_opened=.true.
        ENDIF

        do iwr=1,2

          if (iwr == 1) then
            uwrite=nim_wr
          else
            uwrite=out_unit
          endif

          write (uwrite,*)
          write (uwrite,*)
     $       'Timing statistics:'
          write (uwrite,*)
     $       '                        CPU secs    % of loop'
          write (uwrite,10)'    Seam  time =      ',time_seam_loop,
     $       time_seam_loop/time_loop*100
          write (uwrite,10)'    Seg   time =      ',time_seg_loop,
     $       time_seg_loop/time_loop*100
          write (uwrite,10)'    I/O   time =      ',time_io_loop,
     $       time_io_loop/time_loop*100
          write (uwrite,10)'    Iteration time =  ',time_iter_loop,
     $       time_iter_loop/time_loop*100
          write (uwrite,10)'    Ext solve time =  ',time_exsolv_loop,
     $       time_exsolv_loop/time_loop*100
          write (uwrite,10)'    Factoring time =  ',time_fac_loop,
     $       time_fac_loop/time_loop*100
          write (uwrite,10)'    Ext factor time = ',time_exfact_loop,
     $       time_exfact_loop/time_loop*100
          write (uwrite,10)'    Line comm time =  ',time_line_loop,
     $       time_line_loop/time_loop*100
          write (uwrite,10)'    FFT   time =      ',time_fft_loop,
     $       time_fft_loop/time_loop*100
          write (uwrite,10)'    FE matrix time =  ',time_mat_loop,
     $       time_mat_loop/time_loop*100
          write (uwrite,10)'    FE rhs time =     ',time_rhs_loop,
     $       time_rhs_loop/time_loop*100
          write (uwrite,10)'    Static con time = ',time_stcon_loop,
     $       time_stcon_loop/time_loop*100
          write (uwrite,10)'    Nim matvec time = ',time_matvec_loop,
     $       time_matvec_loop/time_loop*100
          write (uwrite,10)'    QP field time =   ',time_qpfld_loop,
     $       time_qpfld_loop/time_loop*100
          write (uwrite,*)'    # of Seamings in loop =',seamcount_loop
          write (uwrite,*)'    # of SegSeams in loop =',segcount_loop
          write (uwrite,*)
     $       '                  CPU secs   % of total'
          write (uwrite,10)'    Loop  time =',time_loop,
     $       time_loop/time_total*100
          write (uwrite,10)'    Setup time =',time_setup,
     $       time_setup/time_total*100
          write (uwrite,10)'    Total time =',time_total,100._r8

        enddo

      endif

      RETURN
      END SUBROUTINE timer_stats


      END MODULE time
