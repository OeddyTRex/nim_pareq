c      Routines for generating poly input files for the pie and rim regions
c      in formats suitable for the triangle code.
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  write_pie.
c     2.  write_rim.
c-----------------------------------------------------------------------
      SUBROUTINE write_pie
c-----------------------------------------------------------------------
c     subprogram 1. write_pie
c     Writes out the boundary of the pie structure in
c     a format for the triangle code.
c-----------------------------------------------------------------------
      USE local
      USE input
      USE fields
      IMPLICIT NONE
      INTEGER(i4) :: i,ibl,iv
       
c
c						Write the pie file
      OPEN(UNIT=pie_unit,FILE='pie.poly',STATUS='UNKNOWN')
c
c						Number of boundary points.
      WRITE(pie_unit,*)my,' 2 0 0'
c
c						Node id, r, z,block,seam vertex
      iv=0
      DO ibl=1,nybl
        DO i=1, rb(ibl)%my
          iv=iv+1
          WRITE(pie_unit,*)iv,
     &                   rb(ibl)%rz%fs(1,0,i),
     &                   rb(ibl)%rz%fs(2,0,i),
     &                   ibl,
     &                   (rb(ibl)%my+rb(ibl)%mx)*2-i
        ENDDO
      ENDDO
c
c						Number of line segments
      WRITE(pie_unit,*)my,' 1'
c
c						Segment id, node id, node id
      DO i=1, iv-1
        WRITE(pie_unit,fmt='(x,i4,",",i4,",",i4)')i,i,i+1
      ENDDO
      WRITE(pie_unit,fmt='(x,i4,",",i4,", 1")')iv,iv
c
c						No holes
      WRITE(pie_unit,*)'0'
      CLOSE(UNIT=pie_unit)
      END SUBROUTINE write_pie


      SUBROUTINE write_rim(numref)
c-----------------------------------------------------------------------
c     subprogram 2. write_rim
c	Reads the wall.dat file.
c       Creates a set of wedges numbering nbl_rim between the last
c       rblock and the wall and then writes the poly input files for
c       each annuli in a format for the triangle code.
c-----------------------------------------------------------------------
      USE local
      USE input
      USE fields
      IMPLICIT NONE
      INTEGER(i4),INTENT(IN) :: numref
      INTEGER(i4) :: ibl,i,j,k,iv,jv,nbsvert,nisvert,mvert
      INTEGER :: ios
      REAL(r8) :: dmin,dval,deltar,deltaz,dmax
      REAL(r8) :: bmin,bmax,bval,h,g,areamin,areatest
      REAL(r8),DIMENSION(:),ALLOCATABLE :: rbs,zbs   ! outer boundary
      REAL(r8),DIMENSION(:),ALLOCATABLE :: ris,zis   ! inner boundary
      INTEGER(i4),DIMENSION(:,:),ALLOCATABLE :: nis  ! inner boundary
      INTEGER(i4),DIMENSION(:),ALLOCATABLE :: nisref,nbsref,nref
      REAL(r8),DIMENSION(:,:),ALLOCATABLE :: rpoly,zpoly
      INTEGER(i4),DIMENSION(:,:,:),ALLOCATABLE :: ppoly
      INTEGER(i4),DIMENSION(:),ALLOCATABLE :: npoly
      CHARACTER(64) :: tmpfile
      CHARACTER(8) :: kstring
c
c						Read the wall file
      OPEN(UNIT=rim_unit,FILE='wall.dat',STATUS='old',iostat=ios)
      IF(ios  /=  0)RETURN
      READ(rim_unit,*)nbsvert
      ALLOCATE(rbs(nbsvert+1))
      ALLOCATE(zbs(nbsvert+1))
      ALLOCATE(nisref(nbsvert+1))
      DO i=1,nbsvert
        READ(rim_unit,*)rbs(i),zbs(i)
      ENDDO
      CLOSE(rim_unit)
c
c						Make sure the last point
c						repeats the first point.
      IF((rbs(1) /=  rbs(nbsvert)).AND.(zbs(1).NE.zbs(nbsvert)))THEN
        rbs(nbsvert+1)=rbs(nbsvert)
        zbs(nbsvert+1)=zbs(nbsvert)
        nbsvert=nbsvert+1
      ENDIF
c
c						Determine the number of
c						inner boundary vertici.
      nisvert=0
      DO i=1,nybl
        ibl=(nxbl-1)*nybl+i
        nisvert=nisvert+rb(ibl)%my
      ENDDO
      ALLOCATE(ris(nisvert+1))
      ALLOCATE(zis(nisvert+1))
      ALLOCATE(nis(nisvert+1,2))
c
c						Build the list of points
c						for the inner boundary.
      iv=0
      areamin=HUGE(0)
      DO i=1,nybl
        ibl=(nxbl-1)*nybl+i
        DO j=0,rb(ibl)%my-1
          iv=iv+1
          ris(iv)=rb(ibl)%rz%fs(1,rb(ibl)%mx,j)
          zis(iv)=rb(ibl)%rz%fs(2,rb(ibl)%mx,j)
          nis(iv,1)=ibl
          nis(iv,2)=rb(ibl)%mx+j
          areatest=
     &            (rb(ibl)%rz%fs(1,rb(ibl)%mx  ,j+1)-
     &             rb(ibl)%rz%fs(1,rb(ibl)%mx  ,j  ))*
     &            (rb(ibl)%rz%fs(2,rb(ibl)%mx-1,j  )-
     &             rb(ibl)%rz%fs(2,rb(ibl)%mx  ,j  ))-
     &            (rb(ibl)%rz%fs(2,rb(ibl)%mx  ,j+1)-
     &             rb(ibl)%rz%fs(2,rb(ibl)%mx  ,j  ))*
     &            (rb(ibl)%rz%fs(1,rb(ibl)%mx-1,j  )
     &            -rb(ibl)%rz%fs(1,rb(ibl)%mx  ,j  ))
          areatest=ABS(areatest)
          IF(areatest < areamin)areamin=areatest
        ENDDO
      ENDDO
      nisvert=iv
      ris(nisvert+1)=ris(1)
      zis(nisvert+1)=zis(1)
      nis(nisvert+1,1)=nis(1,1)
      nis(nisvert+1,2)=nis(1,2)
      WRITE(out_unit,*)" "
      WRITE(out_unit,*)"Recommended minimum triangle area: ",areamin
      WRITE(out_unit,*)" "
      WRITE(nim_wr,*)" "
      WRITE(nim_wr,*)"Recommended minimum triangle area: ",areamin
      WRITE(nim_wr,*)" "
c-TMP
c     WRITE(6,*)nisvert+1
c     DO i=1,nisvert+1
c       WRITE(6,*)ris(i),zis(i)
c     ENDDO
c
c						Find the minimum length
c						around the inner section.
      dmin=rim_length
      IF(dmin <= 0 )THEN
        dmin=(ris(1)-ris(nisvert))**2+(zis(1)-zis(nisvert))**2
        dmax=0
        DO i=1, nisvert-1
          dval=(ris(i)-ris(i+1))**2+(zis(i)-zis(i+1))**2
          if(dval <  dmin.and.dval.ne.0)dmin=dval
          if(dval >  dmax)dmax=dval
        ENDDO
        dmin=dmin**0.5_r8
        dmax=dmax**0.5_r8
c
c						dmin seems to overrefine
c						the wall so instead use:
        dmin=2.0*dmin
      ENDIF
c
c						Refine the wall to satisfy that
c						each segment is less than dmin.
      mvert=0
      DO i=1,nbsvert-1
        dval=((rbs(i)-rbs(i+1))**2+(zbs(i)-zbs(i+1))**2)**0.5_r8
        nisref(i)=INT(dval/dmin+1.0_r8)
        mvert=mvert+nisref(i)
      ENDDO
      ALLOCATE(npoly(nbl_rim))
      ALLOCATE(rpoly(mvert+nisvert+1,nbl_rim))
      ALLOCATE(zpoly(mvert+nisvert+1,nbl_rim))
      ALLOCATE(ppoly(mvert+nisvert+1,nbl_rim,2))
      ALLOCATE(nbsref(nbl_rim+1))
      ALLOCATE(nref(nbl_rim+1))
      ppoly(:,:,:)=0
c
c						Define the r,z vertici for the
c						wall boundary and use rpoly
c						and zpoly for temporary storage.
      iv=1
      DO i=1,nbsvert-1
        deltar=(rbs(i+1)-rbs(i))/nisref(i)
        deltaz=(zbs(i+1)-zbs(i))/nisref(i)
        DO j=1,nisref(i)
          rpoly(iv,nbl_rim)=rbs(i)+deltar*j
          zpoly(iv,nbl_rim)=zbs(i)+deltaz*j
          iv=iv+1
        ENDDO
      ENDDO
      nbsvert=mvert
      DEALLOCATE(rbs,zbs,nisref)
      ALLOCATE(rbs(nbsvert))
      ALLOCATE(zbs(nbsvert))
      ALLOCATE(nisref(nbl_rim+1))
c
c						Divide the inner boundary
c						into nbl_rim sections and
c						identify the vertex.
      mvert=nisvert/nbl_rim
      DO i=1,nbl_rim
        nisref(i)=(i-1)*mvert+1
      ENDDO
      nisref(nbl_rim+1)=nisvert+1
c                                               Find the closest wall point
c						to the reference starting 
c						point on the wall.
      jv=1
      nbsref(1)=jv
      bmin=(rpoly(jv,nbl_rim)-ris(nisref(1)))**2+
     &       (zpoly(jv,nbl_rim)-zis(nisref(1)))**2

      DO jv=2,nbsvert
        bval=(rpoly(jv,nbl_rim)-ris(nisref(1)))**2+
     &         (zpoly(jv,nbl_rim)-zis(nisref(1)))**2
        IF(bval.LT.bmin)THEN
           bmin=bval
           nbsref(1)=jv
        ENDIF
      ENDDO
c
c						Determine wall rotation
c						direction relative to rblocks.
      g = 0
      DO iv=1,nisvert
        g=g+ris(iv+1)*zis(iv)
     &     -ris(iv)*zis(iv+1)
     &      +0.5*(ris(iv+1)-ris(iv))*(zis(iv+1)-zis(iv))
      ENDDO
      h = 0
      DO iv=1,nbsvert-1
        h=h+rpoly(iv+1,nbl_rim)*zpoly(iv,nbl_rim)-
     &      zpoly(iv+1,nbl_rim)*rpoly(iv,nbl_rim)+
     &      0.5*(zpoly(iv+1,nbl_rim)-zpoly(iv,nbl_rim))*
     &          (rpoly(iv+1,nbl_rim)-rpoly(iv,nbl_rim))
      ENDDO
      IF(g*h < 0)THEN
c						Shift and counter rotate
c						wall vertici.
        jv=0
        DO iv=nbsref(1),1,-1
          jv=jv+1
          rbs(jv)=rpoly(iv,nbl_rim)
          zbs(jv)=zpoly(iv,nbl_rim)
        ENDDO
        DO iv=nbsvert,nbsref(1)+1,-1
          jv=jv+1
          rbs(jv)=rpoly(iv,nbl_rim)
          zbs(jv)=zpoly(iv,nbl_rim)
        ENDDO
      ELSE
c						Shift wall vertici.
        jv=0
        DO iv=nbsref(1),nbsvert
          jv=jv+1
          rbs(jv)=rpoly(iv,nbl_rim)
          zbs(jv)=zpoly(iv,nbl_rim)
        ENDDO
        DO iv=1,nbsref(1)-1
          jv=jv+1
          rbs(jv)=rpoly(iv,nbl_rim)
          zbs(jv)=zpoly(iv,nbl_rim)
        ENDDO
      ENDIF
      rpoly(:,:)=0
      zpoly(:,:)=0
      npoly(:)=0_i4
c                                               Find the closest wall point
c                                               and determine number of
c                                               refinement sections in segment.
      DO i=1,nbl_rim
        nbsref(i)=1
        bmin=(rbs(1)-ris(nisref(i)))**2+
     &       (zbs(1)-zis(nisref(i)))**2
        DO jv=2,nbsvert
          bval=(rbs(jv)-ris(nisref(i)))**2+
     &         (zbs(jv)-zis(nisref(i)))**2
          IF(bval.LT.bmin)THEN
             bmin=bval
             nbsref(i)=jv
          ENDIF
        ENDDO
        dval=bmin**0.5_r8
        nref(i)=INT(dval/dmin+1.0)
        IF(nref(i).LT.2)nref(i)=2_i4
      ENDDO
      nbsref(nbl_rim+1)=nbsvert+1
      nref(nbl_rim+1)=nref(1)
      IF(nbsref(1).NE.1)STOP 10
c
c						Define the first section 
c						between the rblock and the wall.
      IF(nbl_rim.ne.1)THEN
        DO i=1,nbl_rim
          iv=0
          deltar=(rbs(nbsref(i))-ris(nisref(i)))/nref(i)
          deltaz=(zbs(nbsref(i))-zis(nisref(i)))/nref(i)
          DO j=1,nref(i)
            iv=iv+1
            rpoly(iv,i)=ris(nisref(i))+deltar*j
            zpoly(iv,i)=zis(nisref(i))+deltaz*j
          ENDDO
          npoly(i)=nref(i)
        ENDDO
      ELSE
        rpoly(1,1)=rbs(1)
        zpoly(1,1)=zbs(1)
        ppoly(1,1,1)=0
        ppoly(1,1,2)=1
        npoly(1)=1
      ENDIF
c
c                                               Add the wall points.
      DO i=1,nbl_rim
        iv=npoly(i)
        DO jv=nbsref(i)+1,nbsref(i+1)-1
            iv=iv+1
            rpoly(iv,i)=rbs(jv)
            zpoly(iv,i)=zbs(jv)
            ppoly(iv,i,1)=0
            ppoly(iv,i,2)=jv
            npoly(i)=npoly(i)+1
        ENDDO
      ENDDO
c
c                                               Copy the wall to rblock
c						segment from the rblock to
c						wall segment.
      IF(nbl_rim.ne.1)THEN
        DO i=1,nbl_rim-1
          iv=npoly(i)
          DO jv=nref(i+1),1,-1
            iv=iv+1
            rpoly(iv,i)=rpoly(jv,i+1)
            zpoly(iv,i)=zpoly(jv,i+1)
            ppoly(iv,i,1)=i+1+numref
            ppoly(iv,i,2)=jv
            ppoly(jv,i+1,1)=i+numref
            ppoly(jv,i+1,2)=iv
            npoly(i)=npoly(i)+1
            

          ENDDO
        ENDDO
        i=nbl_rim
        iv=npoly(i)
        DO jv=nref(1),1,-1
          iv=iv+1
          rpoly(iv,i)=rpoly(jv,1)
          zpoly(iv,i)=zpoly(jv,1)
          ppoly(iv,i,1)=1+numref
          ppoly(iv,i,2)=jv
          ppoly(jv,1,1)=i+numref
          ppoly(jv,1,2)=iv
          npoly(i)=npoly(i)+1
        ENDDO
      ENDIF
c
c						Define the r,z vertici from
c						inner boundary.
      DO i=1,nbl_rim
        iv=npoly(i)
        DO jv=nisref(i+1),nisref(i),-1
          iv=iv+1
          rpoly(iv,i)=ris(jv)
          zpoly(iv,i)=zis(jv)
          ppoly(iv,i,1)=nis(jv,1)
          ppoly(iv,i,2)=nis(jv,2)
          npoly(i)=npoly(i)+1
        ENDDO
      ENDDO
      IF(nbl_rim.EQ.1)npoly(1)=npoly(1)-1
      
c---------------------------------------------------------------------------
c
c						Write the rim.poly file.
      DO i=1,nbl_rim

c
c                                               Create the filename
        WRITE(kstring,fmt='(i4)')i
        DO
          k= INDEX(TRIM(kstring)," ")
          IF(k == 0)EXIT
          kstring(k:)=kstring(k+1:)
        ENDDO
        tmpfile="rim"//TRIM(kstring)//".poly"
c
c                                               Open the rim?.poly file 
        OPEN(UNIT=rim_unit,FILE=tmpfile, STATUS='replace',iostat=ios)
c
c                                               Number of boundary points.
        WRITE(rim_unit,*)npoly(i),' 2 0 0'
c
c                                               Node id, r, z , ibl, iv
c						Seam information (ibl,iv) is
c						not complete for corner.
        DO jv=1, npoly(i)
          WRITE(rim_unit,fmt='(x,i4,x,e18.12,x,e18.12,x,i4,x,i4)')
     &                     jv,
     &                     rpoly(jv,i),
     &                     zpoly(jv,i),
     &                     ppoly(jv,i,1),
     &                     ppoly(jv,i,2)
        ENDDO
c
c                                               Number of line segments
c                                               Segment id, node id, node id
        IF(nbl_rim.NE.1)THEN
          WRITE(rim_unit,*)npoly(i),' 1'
          DO jv=1, npoly(i)-1
            WRITE(rim_unit,fmt='(x,i4,x,i4,x,i4)')jv,jv,jv+1
          ENDDO
          WRITE(rim_unit,fmt='(x,i4,x,i4,x,"1")')npoly(i),npoly(i)
        ELSE
          WRITE(rim_unit,*)npoly(i),' 1'
          DO jv=1, nbsvert-1
            WRITE(rim_unit,fmt='(x,i4,x,i4,x,i4)')jv,jv,jv+1
          ENDDO
          WRITE(rim_unit,fmt='(x,i4,x,i4," 1")')nbsvert,nbsvert
          DO jv=nbsvert+1, nbsvert+nisvert-1
            WRITE(rim_unit,fmt='(x,i4,x,i4,x,i4)')jv,jv,jv+1
          ENDDO
          WRITE(rim_unit,fmt='(x,i4,x,i4,x,i4)')
     &        nisvert+nbsvert,nisvert+nbsvert,nbsvert+1
        ENDIF
c
c                                               Number of holes
        IF(nbl_rim.NE.1)THEN
          WRITE(rim_unit,fmt='(x,"0")')
          CLOSE(UNIT=rim_unit)
        ELSE
          WRITE(rim_unit,fmt='(x,"1")')
          SELECT CASE(TRIM(pieflag))
          CASE('rblock')
            WRITE(rim_unit,*)" 1",
     &        rb(1)%rz%fs(1,1,1),rb(1)%rz%fs(2,1,1)
          CASE('tblock0')
            WRITE(rim_unit,fmt='(x,"1,",f16.8,",",f16.8)')
     &        tb(nrbl+1)%tgeom%xs(0),tb(nrbl+1)%tgeom%ys(0)
          CASE('tblock1')
            WRITE(rim_unit,fmt='(x,"1,",f16.8,",",f16.8)')
     &        tb(nrbl+1)%tgeom%xs(0),tb(nrbl+1)%tgeom%ys(0)
          END SELECT
        ENDIF
      ENDDO
c
c						Deallocate work arrays
      DEALLOCATE(rbs,zbs,rpoly,zpoly,npoly)
      DEALLOCATE(nbsref,nisref,ppoly)

      RETURN
      END SUBROUTINE write_rim
