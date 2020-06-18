c-----------------------------------------------------------------------
c     program fluxgrid.
c     generates flux grids for NIMROD.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. fluxgrid.
c-----------------------------------------------------------------------
c     subprogram 1. fluxgrid.
c     initializes i/o and controls PROGRAM flow.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM fluxgrid
      USE local
      USE global
      USE input
      USE analyze
      IMPLICIT NONE
      
      TYPE(global_type) :: gt
      LOGICAL :: file_stat
c-----------------------------------------------------------------------
c     start timer, read input, and open output file.
c-----------------------------------------------------------------------
      CALL time_stat(0_i4,2_i4)
      CALL read_input
      WRITE(nim_wr,*)"Equilibrium Type: ",TRIM(eq_type)
      WRITE(nim_wr,*)"Equilibrium File: ",TRIM(filename)
      INQUIRE(FILE=TRIM(filename),EXIST=file_stat)
      IF (.NOT.file_stat) CALL nim_stop
     $  (TRIM(filename)//' does not exist.')
      OPEN(UNIT=out_unit,FILE='fluxgrid.out',STATUS='UNKNOWN')
c-----------------------------------------------------------------------
c     process equilibrium
      WRITE(nim_wr,*) 'Processing Equilibrium'
c-----------------------------------------------------------------------
      NULLIFY(gt%eigvec)
      SELECT CASE(eq_type)
      CASE("efit")
         CALL read_efit(gt)
      CASE("galkin")
         CALL read_galkin(gt)
      CASE("miller")
         CALL read_miller(gt)
      CASE("miller8")
         CALL read_miller8(gt)
      CASE("millasc")
         CALL read_millasc(gt)
      CASE("chum")
         CALL read_chum(gt)
      CASE("soloviev")
         CALL read_soloviev(gt)
      CASE("rsteq")
         CALL read_rsteq(gt)
      CASE("chease")
         CALL read_chease(gt)
      CASE default
         CALL nim_stop('Not a valid equilibrium type')
      END SELECT
c-----------------------------------------------------------------------
c     find useful information about equilibrium  (see analyze.f)
      write(nim_wr,*) 'Finished processing.  Writing output.'
c-----------------------------------------------------------------------
      CALL get_global(gt)
      IF(stability .OR. out_neoclassical .OR. out_pest)    
     &            CALL get_stability(gt)

c-----------------------------------------------------------------------
c     write typical output. (output.f)
c-----------------------------------------------------------------------
      CALL write_out(gt,out_unit,binary_unit)		!fluxgrid.out, ...
      CALL write_2d(gt,ascii_unit,binary_unit)		!2d.out, 2d.bin
      CALL draw_eig(gt,binary_unit)
      IF(out_neoclassical)    CALL write_neoclassical(gt,ascii_unit)
      IF(out_tecplot) CALL write_tecplot(gt)

c-----------------------------------------------------------------------
c     write input file for other codes (codes.f)
c-----------------------------------------------------------------------
      CALL write_nimrod(gt,binary_unit)		        !fluxgrid.dat
      IF(out_pies)    CALL write_pies(gt,binary_unit)
      IF(out_pest)    CALL write_pest(gt,binary_unit)
      IF(out_far)     CALL write_far(gt,binary_unit)

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL time_stat(1_i4,2_i4)
      CLOSE(UNIT=out_unit)
      CALL nim_stop('Normal termination')
      END PROGRAM fluxgrid
