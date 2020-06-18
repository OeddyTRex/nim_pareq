      SUBROUTINE onefl(d_eta,ib_start,x_start,y_start,
     $                 per_start1,nstep1,plane_type1,any_cross)
      USE local
      USE io
      USE dumpc
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: d_eta,x_start,y_start,per_start1
      INTEGER(i4), INTENT(IN) :: nstep1,ib_start
      CHARACTER(*), INTENT(IN) :: plane_type1
      LOGICAL, INTENT(INOUT) :: any_cross

      REAL(r8), DIMENSION(3) :: dep_start,dep_end
      REAL(r8) :: x,y,rc,zc,perc
      INTEGER(i4) :: ib

c     						Initialize position.
      CALL lagr_quad_eval(rbc(ib_start)%rz,x_start,y_start,0_i4)
      dep_start(1:2)=rbc(ib_start)%rz%f
      dep_start(3)=per_start1
      x=x_start
      y=y_start
      ib=ib_start

c						Integrate the field
c						line.
      rc=HUGE(rc)
      zc=HUGE(rc)
      perc=HUGE(rc)
      CALL int_segment(0._r8,d_eta,nstep1,x,y,ib,
     $                 dep_start,dep_end,plane_type1,
     $                 0_i4,rc,zc,perc,any_cross)

      RETURN
      END SUBROUTINE onefl
