c-----------------------------------------------------------------------
c     file nim_locate.f
c     contains a module for locating a position in terms of logical
c     coordinates using block-based searching.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  nim_locate
c     1.  nim_rb_locate_init
c     2.  nim_rb_loc_dealloc
c     3.  nim_rb_locate
c     4.  nim_rb_contrz
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     0. module declaration for nim_locate
c-----------------------------------------------------------------------
      MODULE nim_locate
      USE local
      USE edge_type_mod
      IMPLICIT NONE

      TYPE :: nim_rb_loc_type
        REAL(r8), DIMENSION(:,:), POINTER :: rzmdp
        REAL(r8), DIMENSION(:,:,:), POINTER :: rzedg
        REAL(r8), DIMENSION(:,:,:), POINTER :: tanvc
      END TYPE nim_rb_loc_type

      INTEGER(i4), PARAMETER :: nrbloc=12
      REAL(r8), DIMENSION(0:nrbloc) :: dedge
      TYPE(nim_rb_loc_type), DIMENSION(:), POINTER :: nim_rb_locdata

c!$OMP THREADPRIVATE(dedge,nim_rb_locdata)

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. nim_rb_locate_init
c     initializes and saves data along the border of an rblock that
c     are used to test whether a point is located in the block.  more
c     specifically, this routine sets:
c
c     nedge = number of points to check along the edge of each element
c     dedge = array of element coordinate positions along an edge
c     rzmdp = midpoint between ends of an edge
c     tanvc = tangent vector at the dedge element coordinates along
c             an edge, pointed in the counter-clockwise direction
c             around the block.
c     rzedg = physical coordinates at the dedge element coordinates
c
c     the above information is saved in a nim_rb_loc_type data
c     structure.
c-----------------------------------------------------------------------
      SUBROUTINE nim_rb_locate_init(nrb,rb)
      USE rblock_type_mod

      INTEGER(i4), INTENT(IN) :: nrb
      TYPE(rblock_type), DIMENSION(nrb), INTENT(INOUT) :: rb

      REAL(r8), DIMENSION(0:nrbloc) :: wedge
      REAL(r8), PARAMETER :: doffset=1.e-8_r8

      INTEGER(i4) :: ibl,ix,iy,itst,iv,nvert
      REAL(r8) :: x,y
c-----------------------------------------------------------------------
c     use lobatto spacing along each element edge for testing the
c     cross-product condition.
c-----------------------------------------------------------------------
      CALL lobleg(0._r8,1._r8,dedge,wedge,nrbloc+1_i4)
      dedge(0)=doffset
      dedge(nrbloc)=1._r8-doffset
c-----------------------------------------------------------------------
c     allocate the data structure.
c-----------------------------------------------------------------------
      ALLOCATE(nim_rb_locdata(nrb))
      DO ibl=1,nrb
        iv=0
        nvert=2*(rb(ibl)%mx+rb(ibl)%my)
        ALLOCATE(nim_rb_locdata(ibl)%rzmdp(2,nvert))
        ALLOCATE(nim_rb_locdata(ibl)%rzedg(2,0:nrbloc,nvert))
        ALLOCATE(nim_rb_locdata(ibl)%tanvc(2,0:nrbloc,nvert))
c-----------------------------------------------------------------------
c       block bottom:
c-----------------------------------------------------------------------
        DO ix=1,rb(ibl)%mx
          iv=iv+1
          nim_rb_locdata(ibl)%rzmdp(:,iv)=
     $      (/0.5_r8*(rb(ibl)%rz%fs(1,ix,0)+rb(ibl)%rz%fs(1,ix-1,0)),
     $        0.5_r8*(rb(ibl)%rz%fs(2,ix,0)+rb(ibl)%rz%fs(2,ix-1,0))/)
          y=0._r8
          DO itst=0,nrbloc
            x=REAL(ix-1_i4,r8)+dedge(itst)
            CALL lagr_quad_eval(rb(ibl)%rz,x,y,1_i4)
            nim_rb_locdata(ibl)%rzedg(:,itst,iv)=rb(ibl)%rz%f
            nim_rb_locdata(ibl)%tanvc(:,itst,iv)=rb(ibl)%rz%fx
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       right side:
c-----------------------------------------------------------------------
        DO iy=1,rb(ibl)%my
          iv=iv+1
          nim_rb_locdata(ibl)%rzmdp(:,iv)=
     $      (/0.5_r8*(rb(ibl)%rz%fs(1,rb(ibl)%mx,iy)+
     $                rb(ibl)%rz%fs(1,rb(ibl)%mx,iy-1)),
     $        0.5_r8*(rb(ibl)%rz%fs(2,rb(ibl)%mx,iy)+
     $                rb(ibl)%rz%fs(2,rb(ibl)%mx,iy-1))/)
          x=REAL(rb(ibl)%mx,r8)
          DO itst=0,nrbloc
            y=REAL(iy-1_i4,r8)+dedge(itst)
            CALL lagr_quad_eval(rb(ibl)%rz,x,y,1_i4)
            nim_rb_locdata(ibl)%rzedg(:,itst,iv)=rb(ibl)%rz%f
            nim_rb_locdata(ibl)%tanvc(:,itst,iv)=rb(ibl)%rz%fy
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       block top:
c-----------------------------------------------------------------------
        DO ix=rb(ibl)%mx,1,-1
          iv=iv+1
          nim_rb_locdata(ibl)%rzmdp(:,iv)=
     $      (/0.5_r8*(rb(ibl)%rz%fs(1,ix,rb(ibl)%my)+
     $                rb(ibl)%rz%fs(1,ix-1,rb(ibl)%my)),
     $        0.5_r8*(rb(ibl)%rz%fs(2,ix,rb(ibl)%my)+
     $                rb(ibl)%rz%fs(2,ix-1,rb(ibl)%my))/)
          y=REAL(rb(ibl)%my,r8)
          DO itst=0,nrbloc
            x=REAL(ix,r8)-dedge(itst)
            CALL lagr_quad_eval(rb(ibl)%rz,x,y,1_i4)
            nim_rb_locdata(ibl)%rzedg(:,itst,iv)= rb(ibl)%rz%f
            nim_rb_locdata(ibl)%tanvc(:,itst,iv)=-rb(ibl)%rz%fx
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       left side:
c-----------------------------------------------------------------------
        DO iy=rb(ibl)%my,1,-1
          iv=iv+1
          nim_rb_locdata(ibl)%rzmdp(:,iv)=
     $      (/0.5_r8*(rb(ibl)%rz%fs(1,0,iy)+rb(ibl)%rz%fs(1,0,iy-1)),
     $        0.5_r8*(rb(ibl)%rz%fs(2,0,iy)+rb(ibl)%rz%fs(2,0,iy-1))/)
          x=0._r8
          DO itst=0,nrbloc
            y=REAL(iy,r8)-dedge(itst)
            CALL lagr_quad_eval(rb(ibl)%rz,x,y,1_i4)
            nim_rb_locdata(ibl)%rzedg(:,itst,iv)= rb(ibl)%rz%f
            nim_rb_locdata(ibl)%tanvc(:,itst,iv)=-rb(ibl)%rz%fy
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE nim_rb_locate_init
c-----------------------------------------------------------------------
c     subprogram 2. nim_rb_loc_dealloc
c     deallocate the data structure used to hold boundary-location
c     information for nim_rb searches.
c-----------------------------------------------------------------------
      SUBROUTINE nim_rb_loc_dealloc

      INTEGER(i4) :: ibl
c-----------------------------------------------------------------------
c     loop over the blocks and deallocate the arrays in each, then
c     deallocate the structure.
c-----------------------------------------------------------------------
      DO ibl=1,SIZE(nim_rb_locdata)
        DEALLOCATE(nim_rb_locdata(ibl)%rzmdp)
        DEALLOCATE(nim_rb_locdata(ibl)%rzedg)
        DEALLOCATE(nim_rb_locdata(ibl)%tanvc)
      ENDDO
      DEALLOCATE(nim_rb_locdata)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE nim_rb_loc_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. nim_rb_locate
c     finds the logical coordinates of a point within an rblock using
c     the rz mapping structure.
c
c     parameter list:
c     rpos (real scalar input) = R position of point
c     zpos (real scalar input) = Z position of point
c     scale (real scalar input) = scale for relative error check
c     rb (rblock structure type) = rblock to search
c     x (real scalar input/output) = guess at logical-x on input
c                                    converged logical-x on output
c     y (real scalar input/output) = guess at logical-y on input
c                                    converged logical-y on output
c     ierr (integer output) = 0 if point is converged, nonzero otherwise
c-----------------------------------------------------------------------
      SUBROUTINE nim_rb_locate(rpos,zpos,scale,rb,x,y,ierr)
      USE rblock_type_mod

      REAL(r8), INTENT(IN) :: rpos,zpos,scale
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(INOUT) :: x,y
      INTEGER(i4), INTENT(OUT) :: ierr

      INTEGER(i4), PARAMETER :: nmax=30
      REAL(r8), PARAMETER :: rztol=1.e-11_r8,dyd=0.125_r8,dgtol=1.e-8_r8
      LOGICAL, PARAMETER :: refinedg=.true.

      TYPE(lagr_quad_2D_type) :: laq
      INTEGER(i4) :: irb
      REAL(r8) :: jac,dxdr,dxdz,dydr,dydz,d2,dr,dz,eps,xmx,ymx,xmn,ymn
      REAL(r8) :: rmag,rmag0,rmag1,th,th0,th1,rc,zc,dr0,dz0,dr1,dz1
c-----------------------------------------------------------------------
c     initialize limits.
c-----------------------------------------------------------------------
      eps=0.25_r8
      xmx=rb%mx+0.1_r8
      ymx=rb%my+0.1_r8
      xmn=-0.1_r8
      ymn=-0.1_r8
c-----------------------------------------------------------------------
c     if the input x is less than 0, start from the center of the block,
c     otherwise, use the provided (x,y) as an initial guess.
c
c     for a degenerate block, start with an approximation based on
c     trigonometry.  rays 0 and 1 should bracket the position vector
c     from the degenerate point to (rpos,zpos), and the angle of that
c     vector is computed relative to theta0.
c-----------------------------------------------------------------------
c-CHG
c This assignment is intended to avoid changing rb%rz so that this
c routine is threadsafe. Now all lagr_quad_eval calls should take
c laq as an argument -A.P.S
c-----------------------------------------------------------------------
      
      CALL lagr_quad_2D_alloc(laq,rb%rz%mx,rb%rz%my,rb%rz%nqty,
     $                        rb%rz%n_side+1)
      laq = rb%rz
 
      IF (rb%degenerate) THEN
        rc=laq%fs(1,0,0)
        zc=laq%fs(2,0,0)
        dr0=laq%fs(1,rb%mx,0)-rc
        dz0=laq%fs(2,rb%mx,0)-zc
        dr1=laq%fs(1,rb%mx,rb%my)-rc
        dz1=laq%fs(2,rb%mx,rb%my)-zc
        th0=ATAN2(dz0,dr0)
        th1=ATAN2(dz1,dr1)
        IF (th1-th0<=dgtol) th1=th1+twopi

        rmag0=SQRT(dr0**2+dz0**2)
        rmag1=SQRT(dr1**2+dz1**2)
        rmag =SQRT((rpos-rc)**2+(zpos-zc)**2)
        th=ATAN2(zpos-zc,rpos-rc)
        IF (th<0.5_r8*(th0+th1)-pi) th=th+twopi
        x=2._r8*rb%mx*rmag/(rmag0+rmag1)
        y=rb%my*(th-th0)/(th1-th0)
        IF (rb%periodicy) CALL nim_rb_pershift
c-TMP
c       WRITE(73,*) "rz",REAL((/rpos,zpos,rc,zc/),r4)
c       WRITE(73,*) "  ",REAL(rmag,r4)
c       WRITE(73,*) REAL((/th0,th1,th,x,y/),r4)

        IF (refinedg) THEN
          CALL lagr_quad_eval(laq,x,y-dyd,0_i4)
          dr0=laq%f(1)-rc
          dz0=laq%f(2)-zc
          CALL lagr_quad_eval(laq,x,y+dyd,0_i4)
          dr1=laq%f(1)-rc
          dz1=laq%f(2)-zc
          th0=ATAN2(dz0,dr0)
          th1=ATAN2(dz1,dr1)
          IF (th1-th0<=dgtol) th1=th1+twopi
          th=ATAN2(zpos-zc,rpos-rc)
          IF (th<0.5_r8*(th0+th1)-pi) th=th+twopi


          CALL lagr_quad_eval(laq,REAL(rb%mx,r8),y-dyd,0_i4)
          dr0=laq%f(1)-rc
          dz0=laq%f(2)-zc
          CALL lagr_quad_eval(laq,REAL(rb%mx,r8),y+dyd,0_i4)
          dr1=laq%f(1)-rc
          dz1=laq%f(2)-zc
          rmag0=SQRT(dr0**2+dz0**2)
          rmag1=SQRT(dr1**2+dz1**2)
          x=2._r8*rb%mx*rmag/(rmag0+rmag1)
          y=y-dyd+2._r8*dyd*(th-th0)/(th1-th0)
          IF (rb%periodicy) CALL nim_rb_pershift
c-TMP
c       WRITE(73,*) REAL((/th0,th1,th,x,y/),r4)
        ENDIF
      ELSE IF (x<0._r8) THEN
        x=0.5_r8*rb%mx
        y=0.5_r8*rb%my
      ENDIF
c-----------------------------------------------------------------------
c     use Newton's method to (x,y) such that r(x,y)-rpos=z(x,y)-zpos=0.
c-----------------------------------------------------------------------
      DO irb=0,nmax
c-----------------------------------------------------------------------
c       check error level.
c-----------------------------------------------------------------------
        IF (rb%periodicy) CALL nim_rb_pershift
        CALL lagr_quad_eval(laq,x,y,1_i4)
        dr=rpos-laq%f(1)
        dz=zpos-laq%f(2)
        d2=dr**2+dz**2
        IF (SQRT(d2)<=rztol*scale) THEN
          ierr=0
c-TMP
c        WRITE(73,*) "fnd ",rb%id,REAL((/rpos,zpos,x,y/),r4)
          RETURN
        ENDIF
c-----------------------------------------------------------------------
c       evaluate the Jacobian of the map and the elements of the
c       Jacobian matrix for the inverse map.
c-----------------------------------------------------------------------
        jac=laq%fx(1)*laq%fy(2)-laq%fx(2)*laq%fy(1)
        dxdr= laq%fy(2)/jac
        dxdz=-laq%fy(1)/jac
        dydr=-laq%fx(2)/jac
        dydz= laq%fx(1)/jac
        x=MIN(xmx,MAX(xmn,x+eps*(dxdr*dr+dxdz*dz)))
        y=MIN(ymx,MAX(ymn,y+eps*(dydr*dr+dydz*dz)))
        eps=MIN(2._r8*eps,1._r8)
      ENDDO
      WRITE(nim_wr,*) laq%f(1),rb%rz%f(1)
      ierr=irb
c-TMP
c        WRITE(73,*) "lst ",rb%id,REAL((/rpos,zpos,x,y,eps/),r4)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN

      CONTAINS

c-----------------------------------------------------------------------
c       internal function for periodic shifting.
c-----------------------------------------------------------------------
        SUBROUTINE nim_rb_pershift

        IF (y<0._r8) THEN
          y=y+rb%my
        ELSE IF (y>rb%my) THEN
          y=y-rb%my
        ENDIF

        RETURN
        END SUBROUTINE nim_rb_pershift

      END SUBROUTINE nim_rb_locate
c-----------------------------------------------------------------------
c     subprogram 4. nim_rb_contrz
c     checks whether an rblock contains a physical coordinate pair.
c     the method tests the sign of the third component of
c
c        (R_b-R_rz)XdR_b
c
c     over the block border, where R_rz is the position vector to the
c     (r,z) point, and R_b is the position vector to any point along
c     the block border.
c
c     if [(R_b-R_rz)XdR_b]_3 is positive where |R_b-R_rz| is smallest,
c     the point is enclosed.  unlike the rz_to_cell routine, curved
c     boundaries are considered here.
c
c     parameter list:
c     rpos (real scalar input) = R position of point
c     zpos (real scalar input) = Z position of point
c     rb (rblock structure type) = rblock to check
c     seam (edge structure type) = seam for the block
c     ibl (integer output) = current block id if enclosed, otherwise
c                            neighbor block id closest to (rpos,zpos)
c     iv (integer output) = neighbor seam segment index closest to
c                           (rpos,zpos) if not enclosed.
c-----------------------------------------------------------------------
      SUBROUTINE nim_rb_contrz(rpos,zpos,rb,seam,ibl,iv)
      USE rblock_type_mod
      USE edge_type_mod

      REAL(r8), INTENT(IN) :: rpos,zpos
      TYPE(rblock_type), INTENT(INOUT) :: rb
      TYPE(edge_type), INTENT(INOUT) :: seam
      INTEGER(i4), INTENT(OUT) :: ibl

      REAL(r8), PARAMETER :: doffset=1.e-8_r8

      INTEGER(i4) :: ix,iy,itst,iv,ivmn,iv0,ivp,ivtl,ivtr,ivs
      REAL(r8) :: x,y,relr,relz,dx,dy
      REAL(r4) :: dist2,d2mn     !  lower precision is intentional
      REAL(r4) :: hg4=HUGE(hg4)
      LOGICAL :: outside,outlast,outcur
c-----------------------------------------------------------------------
c     set ibl to the block id, to be changed if the point is not
c     within this block.  also set the number of points to test along
c     each border-element edge.
c-----------------------------------------------------------------------
      ibl=rb%id
      outside=.false.
      d2mn=HUGE(d2mn)
      ivtl=2*rb%mx+rb%my  !  top left corner
      ivtr=  rb%mx+rb%my  !  top right corner
c-----------------------------------------------------------------------
c     first identify the boundary segment that is closest to the point
c     in question.  the midpoints are stored as rzmdp in the 
c     nim_rb_locdata structure.
c-----------------------------------------------------------------------
      IF (rb%periodicy) THEN
        iv=rb%mx+1
      ELSE
        iv=1
      ENDIF
      IF (.NOT.rb%periodicy) THEN
        DO ix=1,rb%mx
          dist2=(nim_rb_locdata(ibl)%rzmdp(1,iv)-rpos)**2+
     $          (nim_rb_locdata(ibl)%rzmdp(2,iv)-zpos)**2
          IF (dist2<=d2mn) THEN
            d2mn=dist2
            iv0=iv
          ENDIF
          iv=iv+1
        ENDDO
      ENDIF
      DO iy=1,rb%my
        dist2=(nim_rb_locdata(ibl)%rzmdp(1,iv)-rpos)**2+
     $        (nim_rb_locdata(ibl)%rzmdp(2,iv)-zpos)**2
        IF (dist2<=d2mn) THEN
          d2mn=dist2
          iv0=iv
        ENDIF
        iv=iv+1
      ENDDO
      IF (.NOT.rb%periodicy) THEN
        DO ix=rb%mx,1,-1
          dist2=(nim_rb_locdata(ibl)%rzmdp(1,iv)-rpos)**2+
     $          (nim_rb_locdata(ibl)%rzmdp(2,iv)-zpos)**2
          IF (dist2<=d2mn) THEN
            d2mn=dist2
            iv0=iv
          ENDIF
          iv=iv+1
        ENDDO
      ELSE
        iv=ivtl+1
      ENDIF
      IF (.NOT.rb%degenerate) THEN
        DO iy=rb%my,1,-1
          dist2=(nim_rb_locdata(ibl)%rzmdp(1,iv)-rpos)**2+
     $          (nim_rb_locdata(ibl)%rzmdp(2,iv)-zpos)**2
          IF (dist2<=d2mn) THEN
            d2mn=dist2
            iv0=iv
          ENDIF
          iv=iv+1
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     determine the segment before the previous one.
c-----------------------------------------------------------------------
      ivs=iv0-2
      IF (ivs<1) ivs=ivs+seam%nvert
      IF (rb%periodicy) THEN
        IF (ivs<=rb%mx) THEN
          ivs=ivtr-rb%mx+ivs
        ELSE IF (ivs<=ivtl.AND.ivs>ivtr) THEN
          ivs=seam%nvert-ivtl+ivs
        ENDIF
      ELSE IF (rb%degenerate) THEN  !  no need for per & degen condition
        IF (ivs>ivtl) ivs=ivtl-seam%nvert+ivs
      ENDIF
c-TMP
c     ivp=ivs
c-----------------------------------------------------------------------
c     check all five sides centered about iv0.
c-----------------------------------------------------------------------
      outlast=.false.
      d2mn=hg4
      DO iv=1,4
        DO itst=0,nrbloc
          relr=nim_rb_locdata(ibl)%rzedg(1,itst,ivs)-rpos
          relz=nim_rb_locdata(ibl)%rzedg(2,itst,ivs)-zpos
          dist2=relr**2+relz**2
          IF (dist2<=d2mn) THEN
            outcur=relr*nim_rb_locdata(ibl)%tanvc(2,itst,ivs)-
     $             relz*nim_rb_locdata(ibl)%tanvc(1,itst,ivs)<0._r8
            outside=outcur.OR.outlast
            IF (outcur.OR..NOT.outlast) THEN
              d2mn=dist2
              ivmn=ivs
            ENDIF
          ENDIF
          outlast=.false.
        ENDDO
        IF (dist2<=d2mn)
     $    outlast=(relr*nim_rb_locdata(ibl)%tanvc(2,nrbloc,ivs)-
     $             relz*nim_rb_locdata(ibl)%tanvc(1,nrbloc,ivs)<0._r8)
        ivs=ivs+1
        IF (ivs>seam%nvert) ivs=1
        IF (rb%periodicy) THEN
          IF (ivs==1) THEN
            ivs=ivtl+1
          ELSE IF (ivs>ivtr.AND.ivs<=ivtl) THEN
            ivs=rb%mx+1
          ENDIF
        ELSE IF (rb%degenerate) THEN
          IF (ivs>ivtl) ivs=1
        ENDIF
      ENDDO

      DO itst=0,nrbloc  !  last side doesn't need the ivs update
        relr=nim_rb_locdata(ibl)%rzedg(1,itst,ivs)-rpos
        relz=nim_rb_locdata(ibl)%rzedg(2,itst,ivs)-zpos
        dist2=relr**2+relz**2
        IF (dist2<=d2mn) THEN
          outcur=relr*nim_rb_locdata(ibl)%tanvc(2,itst,ivs)-
     $           relz*nim_rb_locdata(ibl)%tanvc(1,itst,ivs)<0._r8
          outside=outcur.OR.outlast
          IF (outcur.OR..NOT.outlast) THEN
            d2mn=dist2
            ivmn=ivs
          ENDIF
        ENDIF
        outlast=.false.
      ENDDO
c-----------------------------------------------------------------------
c     set neighbor block and seam segment index if (rpos,zpos) is not
c     enclosed.
c
c     if the connected block is the same as this block, the block is
c     periodic, and the outside test may be true at a corner.  in this
c     case, the previous seam connection will be outside this block.
c-----------------------------------------------------------------------
c-TMP
c     WRITE(78,*) ibl,outside,ivmn,d2mn,ivp,iv0,
c    $  REAL((/rpos,zpos/),r4)

      IF (outside) THEN
        ibl=seam%segment(ivmn)%ptr(1)
        IF (ibl==rb%id) THEN
          ivmn=ivmn-1
          IF (ivmn==0) ivmn=seam%nvert
          ibl=seam%segment(ivmn)%ptr(1)
        ENDIF
        iv =seam%segment(ivmn)%ptr(2)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE nim_rb_contrz
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE nim_locate
