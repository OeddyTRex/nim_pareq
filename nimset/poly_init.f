c-----------------------------------------------------------------------
c     file poly_init.f
c     performs initialization for triangulated polygon regions. 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. poly_init
c     1. poly_pie_init
c     2. poly_pie_seam_init
c     3. poly_rim_init
c     4. poly_rim_seam_init
c-----------------------------------------------------------------------
c     subprogram 0. poly_init.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE poly_init
      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE edge_type_mod
      USE input
      USE physdat
      IMPLICIT NONE

      CONTAINS

      SUBROUTINE poly_pie_init(tb)
c-----------------------------------------------------------------------
c     subprogram 1. poly_pie_init.
c
c						Intended to read in the
c						node and element files
c						for a pie structure for a 
c						format from the triangle code.
c-----------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(tblock_type), DIMENSION(:), POINTER :: tb

      INTEGER(i4) :: ivert,mvert,icell,mcell,kb,iseg,mseg
      INTEGER(i4) :: idum,i,ivert1,ivert2
      INTEGER(i4) :: ix,iy,iter_loop
      INTEGER(i4) :: iverta,ivertb,ivertc
      INTEGER(i4) :: max_loop=100
      INTEGER(i4),ALLOCATABLE,DIMENSION(:) :: ibnd1
      REAL(r8) :: r,z,xval,yval,dr,dz,dx,dy,jacobian,err
      REAL(r8) :: tol=1.e-08
      TYPE(bicube_type) :: work
      LOGICAL :: file_stat

 
c
c						Get the number of vertex
c						and the number of cells.
      INQUIRE(FILE='pie.1.node',EXIST=file_stat)
      IF (.NOT.file_stat) CALL nim_stop('File pie.1.node does not '//
     $  'exist.  First, run nimset with pieflag=tblock0, '//
     $  'run triangle, rerun fluxgrid, then nimset with '//
     $  'pieflag=tblock1.')
      OPEN(UNIT=pie_unit,FILE='pie.1.node',STATUS='OLD')
      READ(pie_unit,*)mvert,idum,idum,idum
      mvert=mvert-1_i4
      CLOSE(UNIT=pie_unit)
      INQUIRE(FILE='pie.1.ele',EXIST=file_stat)
      IF (.NOT.file_stat) CALL nim_stop('File pie.1.ele does not '//
     $  'exist.  First, run nimset with pieflag=tblock0, '//
     $  'run triangle, rerun fluxgrid, then nimset with '//
     $  'pieflag=tblock1.')
      OPEN(UNIT=pie_unit,FILE='pie.1.ele',STATUS='OLD')
      READ(pie_unit,*)mcell,idum,idum
      CLOSE(UNIT=pie_unit)
c
c						Allocate new pie block
      kb=nxbl*nybl+1
      CALL tri_linear_geom_alloc(tb(kb)%tgeom,mvert,mcell)
      tb(kb)%mvert=mvert
      tb(kb)%mcell=mcell
      CALL tri_linear_alloc(tb(kb)%be_eq,mvert,3_i4)
      tb(kb)%be_eq%name='bq'
      tb(kb)%be_eq%title=' be_eq'
      CALL tri_linear_alloc(tb(kb)%ja_eq,mvert,3_i4)
      tb(kb)%ja_eq%name='jq'
      tb(kb)%ja_eq%title=' ja_eq'
      CALL tri_linear_alloc(tb(kb)%ve_eq,mvert,3_i4)
      tb(kb)%ve_eq%name='vq'
      tb(kb)%ve_eq%title=' ve_eq'
      CALL tri_linear_alloc(tb(kb)%pres_eq,mvert,1_i4)
      tb(kb)%pres_eq%name='pq'
      tb(kb)%pres_eq%title=' pr_eq'
      CALL tri_linear_alloc(tb(kb)%prese_eq,mvert,1_i4)
      tb(kb)%prese_eq%name='pq'
      tb(kb)%prese_eq%title='pre_eq'
      CALL tri_linear_alloc(tb(kb)%nd_eq,mvert,1_i4)
      tb(kb)%nd_eq%name='nq'
      tb(kb)%nd_eq%title=' nd_eq'
      CALL tri_linear_alloc(tb(kb)%diff_shape,mvert,1_i4)
      tb(kb)%diff_shape%name='ds'
      tb(kb)%diff_shape%title='dif_sh'

c
c
c						Get the coordinates of
c						each vertex
c						id,r,z,block,seam vert.
c						it's a boundary cell.
c     tb(kb)%tgeom%xs(0)=xo
c     tb(kb)%tgeom%ys(0)=yo
      OPEN(UNIT=pie_unit,FILE='pie.1.node',STATUS='OLD')
      READ(pie_unit,*)idum,idum,idum,idum
      DO i=0,mvert
        READ(pie_unit,*)ivert,
     &                  tb(kb)%tgeom%xs(ivert-1),
     &                  tb(kb)%tgeom%ys(ivert-1),
     &                  idum
      ENDDO
      CLOSE(UNIT=pie_unit)
c
c						Get the verticies for
c						each cell.
      OPEN(UNIT=pie_unit,FILE='pie.1.ele',STATUS='OLD')
      READ(pie_unit,*)idum,idum,idum
      DO i=1, mcell
        READ(pie_unit,*)icell,iverta,ivertb,ivertc
        tb(kb)%tgeom%vertex(icell,1)=iverta-1_i4
        tb(kb)%tgeom%vertex(icell,2)=ivertb-1_i4
        tb(kb)%tgeom%vertex(icell,3)=ivertc-1_i4
      ENDDO
      CLOSE(UNIT=pie_unit)
c
c						Build the list of neighbors.
c						from the list of line segments.
c							determine vertex
c							bounds on first pass
      INQUIRE(FILE='pie.1.edge',EXIST=file_stat)
      IF (.NOT.file_stat) CALL nim_stop('File pie.1.edge does not '//
     $  'exist.  First, run nimset with pieflag=tblock0, '//
     $  'run triangle, rerun fluxgrid, then nimset with '//
     $  'pieflag=tblock1.')
      OPEN(UNIT=pie_unit,FILE='pie.1.edge',STATUS='OLD')
      READ(pie_unit,*)mseg,idum
      ALLOCATE(ibnd1(0:mvert))
      ibnd1(:)=0_i4
      DO i=1,mseg
        READ(pie_unit,*)iseg,ivert1,ivert2,idum
        ivert1=ivert1-1_i4
        ivert2=ivert2-1_i4
        ibnd1(ivert1)=ibnd1(ivert1)+1
        ibnd1(ivert2)=ibnd1(ivert2)+1
      ENDDO
      CLOSE(UNIT=pie_unit)
c
c							Second pass sets
c							vertexes
      OPEN(UNIT=pie_unit,FILE='pie.1.edge',STATUS='OLD')
      READ(pie_unit,*)mseg,idum
      DO i=0,mvert
        ALLOCATE(tb(kb)%tgeom%neighbor(i)%vertex(0:ibnd1(i)))
        tb(kb)%tgeom%neighbor(i)%vertex(0)=i
        ibnd1(i)=0_i4
      ENDDO
      DO i=1,mseg
        READ(pie_unit,*)iseg,ivert1,ivert2,idum
        ivert1=ivert1-1_i4
        ivert2=ivert2-1_i4
        ibnd1(ivert1)=ibnd1(ivert1)+1_i4
        ibnd1(ivert2)=ibnd1(ivert2)+1_i4
        tb(kb)%tgeom%neighbor(ivert1)%vertex(ibnd1(ivert1))=ivert2
        tb(kb)%tgeom%neighbor(ivert2)%vertex(ibnd1(ivert2))=ivert1
      ENDDO
      CLOSE(UNIT=pie_unit)
c
c						Read the equilibrium magnetic
c						field as previously processed
c						by fluxgrid.
      tb(kb)%ve_eq=0
      tb(kb)%prese_eq=0
      tb(kb)%ja_eq=0
      tb(kb)%nd_eq=ndens
      tb(kb)%diff_shape=1
      INQUIRE(FILE='pie.dat',EXIST=file_stat)
      IF (.NOT.file_stat) CALL nim_stop('File pie.dat does not '//
     $  'exist.  First, run nimset with pieflag=tblock0, '//
     $  'run triangle, rerun fluxgrid, then nimset with '//
     $  'pieflag=tblock1.')
      OPEN(UNIT=pie_unit,FILE='pie.dat',STATUS='OLD')
      DO i=0,mvert
          READ(pie_unit,*)tb(kb)%be_eq%fs(1,i,0),
     &                    tb(kb)%be_eq%fs(2,i,0),
     &                    tb(kb)%be_eq%fs(3,i,0),
     &                    tb(kb)%pres_eq%fs(1,i,0)
          tb(kb)%be_eq%fs(3,i,0)=tb(kb)%be_eq%fs(3,i,0)
     &                          *tb(kb)%tgeom%xs(i)
      ENDDO
      CLOSE(pie_unit)
      DEALLOCATE(ibnd1)
      RETURN
      END SUBROUTINE poly_pie_init


      SUBROUTINE poly_pie_seam_init(rb,tb,seam)
c-----------------------------------------------------------------------
c     subprogram 2. poly_pie_seam_init.
c     initializes the seam for the pie central block.
c-----------------------------------------------------------------------

      TYPE(rblock_type), DIMENSION(:), POINTER :: rb
      TYPE(tblock_type), DIMENSION(:), POINTER :: tb
      TYPE(edge_type), DIMENSION(:), POINTER :: seam

      INTEGER(i4) :: ixbl,iybl,ib,lx,ly,iv,ix,iy,ip,nv,mx1,mx2
      INTEGER(i4) :: kv,kb,nqty,jvold,jbold,mvert,is,lv
      INTEGER(i4) :: i,idum,ibold,ivold
      REAL(r8) :: rdum
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: np
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: jb,jv
      REAL(r8), DIMENSION(:), ALLOCATABLE :: ixv,iyv
      LOGICAL :: file_stat
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 2010 FORMAT('rblock(',i2.2,'), lx = ',i3,', ly = ',i3,', mseam = ',i3)
 2020 FORMAT(/4x,'iv',4x,'ip',4x,'ix',4x,'iy',3x,'outb',2x,'outv'/)
 2030 FORMAT(6i6)
 2040 FORMAT('tblock(',i2.2,'), mvert = ',i3,', mseam = ',i3)
c-----------------------------------------------------------------------
c     start loops over rblocks.
c-----------------------------------------------------------------------
      ib=0
      DO ixbl=1,nxbl
         DO iybl=1,nybl
            ib=ib+1
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
            lx=rb(ib)%mx
            ly=rb(ib)%my
            nv=2*lx+2*ly
            seam(ib)%nvert=nv
            ALLOCATE(seam(ib)%vertex(nv))
            ALLOCATE(jb(nv,3))
            ALLOCATE(jv(nv,3))
            ALLOCATE(np(nv))
            ALLOCATE(ixv(2*(lx+ly)))
            ALLOCATE(iyv(2*(lx+ly)))
            np=(/(1,ix=1,lx-1),3,(1,iy=1,ly-1),3,
     $           (1,ix=1,lx-1),3,(1,iy=1,ly-1),3/)
            ixv=(/(ix,ix=1,lx),(lx,iy=1,ly),
     $           (ix,ix=lx-1,0,-1),(0_i4,iy=ly-1,0,-1)/)
            iyv=(/(0_i4,ix=1,lx),(iy,iy=1,ly),
     $           (ly,ix=lx-1,0,-1),(iy,iy=ly-1,0,-1)/)
c-----------------------------------------------------------------------
c     fill block pointers.
c-----------------------------------------------------------------------
            jb=0
            jb(:,1)=(/(ib-1_i4,ix=1,lx),(ib+nybl,iy=1,ly),
     $           (ib+1_i4,ix=1,lx),(ib-nybl,iy=1,ly)/)
            jb(lx,2)=ib-1+nybl
            jb(lx,3)=ib+nybl
            jb(lx+ly,2)=ib+nybl+1
            jb(lx+ly,3)=ib+1
            jb(2*lx+ly,2)=ib+1-nybl
            jb(2*lx+ly,3)=ib-nybl
            jb(2*lx+2*ly,2)=ib-nybl-1
            jb(2*lx+2*ly,3)=ib-1
c-----------------------------------------------------------------------
c     trim left pointers.
c-----------------------------------------------------------------------
            IF(ixbl.EQ.1)THEN
               jb(2*lx+ly+1:2*(lx+ly),1)=0
               jb(2*lx+ly,2:3)=0
               jb(2*(lx+ly),2)=0
            ENDIF
c-----------------------------------------------------------------------
c     trim right pointers.
c-----------------------------------------------------------------------
            IF(ixbl.EQ.nxbl)THEN
               jb(lx+1:lx+ly,1)=0
               jb(lx,2:3)=0
               jb(lx+ly,2)=0
            ENDIF
c-----------------------------------------------------------------------
c     trim bottom pointers.
c-----------------------------------------------------------------------
            IF(iybl.EQ.1)THEN
               jb(1:lx,1)=0
               jb(2*(lx+ly),2:3)=0
               jb(lx,2)=0
            ENDIF
c-----------------------------------------------------------------------
c     trim top pointers.
c-----------------------------------------------------------------------
            IF(iybl.EQ.nybl)THEN
               jb(1+lx+ly:2*lx+ly,1)=0
               jb(lx+ly,2:3)=0
               jb(2*lx+ly,2)=0
            ENDIF
c-----------------------------------------------------------------------
c     fill vertex pointers, non-corners.
c-----------------------------------------------------------------------
            jv=0
            DO iv=1,lx
               IF(jb(iv,1).NE.0)
     $              jv(iv,1)=2*lx+rb(jb(iv,1))%my-iv
               IF(jb(iv+lx+ly,1).NE.0)
     $              jv(iv+lx+ly,1)=lx-iv
            ENDDO
            DO iv=1,ly
               IF(jb(iv+lx,1).NE.0)
     $              jv(iv+lx,1)=2*rb(jb(iv+lx,1))%mx+2*ly-iv
               IF(jb(iv+2*lx+ly,1).NE.0)
     $              jv(iv+2*lx+ly,1)=rb(jb(iv+2*lx+ly,1))%mx+ly-iv
            ENDDO
c-----------------------------------------------------------------------
c     fill vertex pointers, corners.
c-----------------------------------------------------------------------
            IF(jb(lx,2).NE.0)
     $           jv(lx,2)=2*rb(jb(lx,2))%mx+rb(jb(lx,2))%my
            IF(jb(lx,3).NE.0)
     $           jv(lx,3)=2*rb(jb(lx,3))%mx+2*rb(jb(lx,3))%my
            IF(jb(lx+ly,2).NE.0)
     $           jv(lx+ly,2)=2*rb(jb(lx+ly,2))%mx+2*rb(jb(lx+ly,2))%my
            IF(jb(lx+ly,3).NE.0)
     $           jv(lx+ly,3)=rb(jb(lx+ly,3))%mx
            IF(jb(2*lx+ly,1).NE.0)
     $           jv(2*lx+ly,1)=2*(rb(jb(2*lx+ly,1))%mx
     $           +rb(jb(2*lx+ly,1))%my)
            IF(jb(2*lx+ly,2).NE.0)
     $           jv(2*lx+ly,2)=rb(jb(2*lx+ly,2))%mx
            IF(jb(2*lx+ly,3).NE.0)
     $           jv(2*lx+ly,3)=rb(jb(2*lx+ly,3))%mx+rb(jb(2*lx+ly,3))%my
            IF(jb(2*(lx+ly),2).NE.0)
     $           jv(2*(lx+ly),2)=rb(jb(2*(lx+ly),2))%mx
     $           +rb(jb(2*(lx+ly),2))%my
            IF(jb(2*(lx+ly),3).NE.0)
     $           jv(2*(lx+ly),3)=2*rb(jb(2*(lx+ly),3))%mx
     $           +rb(jb(2*(lx+ly),3))%my
c-----------------------------------------------------------------------
c     fill seam pointers.
c-----------------------------------------------------------------------
            DO iv=1,nv
               ALLOCATE(seam(ib)%vertex(iv)%ptr(2,np(iv)))
               seam(ib)%vertex(iv)%intxy=(/ixv(iv),iyv(iv)/)
               DO ip=1,np(iv)
                  seam(ib)%vertex(iv)%ptr(1,ip)=jb(iv,ip)
                  seam(ib)%vertex(iv)%ptr(2,ip)=jv(iv,ip)
               ENDDO
            ENDDO
c-----------------------------------------------------------------------
c     finish loop over blocks.
c-----------------------------------------------------------------------
            DEALLOCATE(np)
            DEALLOCATE(jb)
            DEALLOCATE(jv)
            DEALLOCATE(ixv)
            DEALLOCATE(iyv)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     periodic boundary conditions in theta.
c-----------------------------------------------------------------------
      ib=1
      kb=nybl
      kv=2*rb(kb)%mx+rb(kb)%my
      DO ixbl=1,nxbl
         iv=SIZE(seam(ib)%vertex)
         seam(ib)%vertex(iv)%ptr(1,3)=kb
         seam(ib)%vertex(iv)%ptr(2,3)=kv
         seam(kb)%vertex(kv)%ptr(1,1)=ib
         seam(kb)%vertex(kv)%ptr(2,1)=iv
         DO iv=1,rb(ib)%mx-1
            kv=kv-1
            seam(ib)%vertex(iv)%ptr(1,1)=kb
            seam(ib)%vertex(iv)%ptr(2,1)=kv
            seam(kb)%vertex(kv)%ptr(1,1)=ib
            seam(kb)%vertex(kv)%ptr(2,1)=iv
         ENDDO
         iv=rb(ib)%mx
         kv=kv-1
         seam(ib)%vertex(iv)%ptr(1,1)=kb
         seam(ib)%vertex(iv)%ptr(2,1)=kv
         seam(kb)%vertex(kv)%ptr(1,3)=ib
         seam(kb)%vertex(kv)%ptr(2,3)=iv
         IF(ixbl.lt.nxbl)THEN
            seam(ib+nybl)%vertex(SIZE(seam(ib+nybl)%vertex))%ptr(1,2)=kb
            seam(ib+nybl)%vertex(SIZE(seam(ib+nybl)%vertex))%ptr(2,2)=kv
            seam(kb)%vertex(kv)%ptr(1,2)=ib+nybl
            seam(kb)%vertex(kv)%ptr(2,2)=SIZE(seam(ib+nybl)%vertex)
            kb=kb+nybl
            kv=2*rb(kb)%mx+rb(kb)%my
            seam(ib)%vertex(iv)%ptr(1,2)=kb
            seam(ib)%vertex(iv)%ptr(2,2)=kv
            seam(kb)%vertex(kv)%ptr(1,2)=ib
            seam(kb)%vertex(kv)%ptr(2,2)=iv
            ib=ib+nybl
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     diagnose rblock seams before pie.
c-----------------------------------------------------------------------
      IF(detflag)THEN
         IF(SIZE(tb).NE.0)WRITE(out_unit,'(a)')"Before pie:"
         DO ib=1,SIZE(rb)
            WRITE(out_unit,2010)ib,rb(ib)%mx,
     &               rb(ib)%my,SIZE(seam(ib)%vertex)
            WRITE(out_unit,2020)
            DO iv=1,SIZE(seam(ib)%vertex)
               DO ip=1,SIZE(seam(ib)%vertex(iv)%ptr,2)
                  WRITE(out_unit,2030)iv,ip,
     &                 seam(ib)%vertex(iv)%intxy,
     &                 seam(ib)%vertex(iv)%ptr(1,ip),
     &                 seam(ib)%vertex(iv)%ptr(2,ip)
               ENDDO
            ENDDO
            WRITE(out_unit,2020)
         ENDDO
      ENDIF
c
c     						set up pointers between 
c						rblocks and pie, non-corners.
c						Everything but the corners
c						is in the pie.poly file.
      kb=nxbl*nybl+1
      INQUIRE(FILE='pie.poly',EXIST=file_stat)
      IF (.NOT.file_stat) CALL nim_stop('File pie.poly does not '//
     $  'exist.  First, run nimset with pieflag=tblock0, '//
     $  'run triangle, rerun fluxgrid, then nimset with '//
     $  'pieflag=tblock1.')
      OPEN(UNIT=pie_unit,FILE='pie.poly',STATUS='OLD')
      READ(pie_unit,*)mvert,idum,idum,idum
      ALLOCATE(seam(kb)%vertex(mvert))
      seam(kb)%nvert=mvert

      DO i=1,mvert
        READ(pie_unit,fmt=*)kv,rdum,rdum,ib,iv
        seam(kb)%vertex(kv)%intxy=(/kv-1_i4,0_i4/)
        ALLOCATE(seam(kb)%vertex(kv)%ptr(2,1))
c                                               connections on seam being read.
        seam(kb)%vertex(kv)%ptr(1,1)=ib
        seam(kb)%vertex(kv)%ptr(2,1)=iv
c						connections on the rblock seam
        seam(ib)%vertex(iv)%ptr(1,1)=kb
        seam(ib)%vertex(iv)%ptr(2,1)=kv
      ENDDO
      CLOSE(pie_unit)
c
c						Identify the corner seams.
c						(I search but I don't really
c						need to.)
      DO kv=1,mvert-1
        IF(seam(kb)%vertex(kv)%ptr(1,1).NE.
     &     seam(kb)%vertex(kv+1)%ptr(1,1))THEN
           ibold=seam(kb)%vertex(kv)%ptr(1,1)
           ivold=seam(kb)%vertex(kv)%ptr(2,1)
           jbold=seam(kb)%vertex(kv+1)%ptr(1,1)
           jvold=seam(kb)%vertex(kv+1)%ptr(2,1)+1
           DEALLOCATE(seam(kb)%vertex(kv)%ptr)
           DEALLOCATE(seam(ibold)%vertex(ivold)%ptr)
           DEALLOCATE(seam(jbold)%vertex(jvold)%ptr)
           ALLOCATE(seam(kb)%vertex(kv)%ptr(2,2))
           ALLOCATE(seam(ibold)%vertex(ivold)%ptr(2,2))
           ALLOCATE(seam(jbold)%vertex(jvold)%ptr(2,2))

           seam(kb)%vertex(kv)%ptr(1,1)=ibold
           seam(kb)%vertex(kv)%ptr(2,1)=ivold
           seam(kb)%vertex(kv)%ptr(1,2)=jbold
           seam(kb)%vertex(kv)%ptr(2,2)=jvold

           seam(ibold)%vertex(ivold)%ptr(1,1)=kb
           seam(ibold)%vertex(ivold)%ptr(2,1)=kv
           seam(ibold)%vertex(ivold)%ptr(1,2)=jbold
           seam(ibold)%vertex(ivold)%ptr(2,2)=jvold

           seam(jbold)%vertex(jvold)%ptr(1,1)=kb
           seam(jbold)%vertex(jvold)%ptr(2,1)=kv
           seam(jbold)%vertex(jvold)%ptr(1,2)=ibold
           seam(jbold)%vertex(jvold)%ptr(2,2)=ivold
        ENDIF
      ENDDO
      kv=mvert
      ibold=seam(kb)%vertex(kv)%ptr(1,1)
      ivold=seam(kb)%vertex(kv)%ptr(2,1)
      jbold=seam(kb)%vertex(1)%ptr(1,1)
      jvold=seam(kb)%vertex(1)%ptr(2,1)+1

      DEALLOCATE(seam(kb)%vertex(kv)%ptr)
      DEALLOCATE(seam(ibold)%vertex(ivold)%ptr)
      DEALLOCATE(seam(jbold)%vertex(jvold)%ptr)
      ALLOCATE(seam(kb)%vertex(kv)%ptr(2,2))
      ALLOCATE(seam(ibold)%vertex(ivold)%ptr(2,2))
      ALLOCATE(seam(jbold)%vertex(jvold)%ptr(2,2))

      seam(kb)%vertex(kv)%ptr(1,1)=ibold
      seam(kb)%vertex(kv)%ptr(2,1)=ivold
      seam(kb)%vertex(kv)%ptr(1,2)=jbold
      seam(kb)%vertex(kv)%ptr(2,2)=jvold

      seam(ibold)%vertex(ivold)%ptr(1,1)=kb
      seam(ibold)%vertex(ivold)%ptr(2,1)=kv
      seam(ibold)%vertex(ivold)%ptr(1,2)=jbold
      seam(ibold)%vertex(ivold)%ptr(2,2)=jvold

      seam(jbold)%vertex(jvold)%ptr(1,1)=kb
      seam(jbold)%vertex(jvold)%ptr(2,1)=kv
      seam(jbold)%vertex(jvold)%ptr(1,2)=ibold
      seam(jbold)%vertex(jvold)%ptr(2,2)=ivold


c
c     						diagnose rblock seams after pie.
      IF(detflag.AND.SIZE(tb).NE.0)THEN
         WRITE(out_unit,'(a)')"After pie, before rim."
         is=0
         DO ib=1,SIZE(rb)
            is=is+1
            WRITE(out_unit,2010)ib,rb(ib)%mx,
     &            rb(ib)%my,SIZE(seam(is)%vertex)
            WRITE(out_unit,2020)
            DO iv=1,SIZE(seam(is)%vertex)
               DO ip=1,SIZE(seam(is)%vertex(iv)%ptr,2)
                  WRITE(out_unit,2030)iv,ip,
     &                 seam(is)%vertex(iv)%intxy,
     &                 seam(is)%vertex(iv)%ptr(1,ip),
     &                 seam(is)%vertex(iv)%ptr(2,ip)
               ENDDO
            ENDDO
            WRITE(out_unit,2020)
         ENDDO
c
c     						diagnose tblock seams after pie.
         DO ib=SIZE(rb)+1_i4,SIZE(rb)+1_i4
            is=is+1
            WRITE(out_unit,2040)ib,tb(ib)%mvert,SIZE(seam(is)%vertex)
            WRITE(out_unit,2020)
            DO iv=1,SIZE(seam(is)%vertex)
               DO ip=1,SIZE(seam(is)%vertex(iv)%ptr,2)
                  WRITE(out_unit,2030)iv,ip,
     &                 seam(is)%vertex(iv)%intxy,
     &                 seam(is)%vertex(iv)%ptr(1,ip),
     &                 seam(is)%vertex(iv)%ptr(2,ip)
               ENDDO
            ENDDO
            WRITE(out_unit,2020)
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE poly_pie_seam_init


      SUBROUTINE poly_rim_init(nblref,tb,nv_ext)
c-----------------------------------------------------------------------
c     subprogram 3. poly_rim_init.
c
c						Intended to read in the
c						node and element files
c						for a rim structure in a 
c						format from the triangle code.
c-----------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(tblock_type), DIMENSION(:), POINTER :: tb
      INTEGER(i4), INTENT(OUT) :: nv_ext
      INTEGER(i4), INTENT(IN) :: nblref

      INTEGER(i4) :: ivert,mvert,icell,mcell,kb,iseg,mseg,k,ibl
      INTEGER(i4) :: idum,i,ivert1,ivert2
      INTEGER(i4) :: iverta,ivertb,ivertc
      INTEGER(i4) :: ix,iy
      INTEGER(i4),ALLOCATABLE,DIMENSION(:) :: ibnd,ibnd1
      CHARACTER(64) :: polyfile,nodefile,elefile,edgefile,datfile
      CHARACTER(8) :: kstring
      LOGICAL :: file_stat

c
c						Loop through each rim tblock
      DO ibl=1,nbl_rim
c
c						Set filenames
        WRITE(kstring,fmt='(i4)')ibl
        DO
          k= INDEX(TRIM(kstring)," ")
          IF(k == 0)EXIT
          kstring(k:)=kstring(k+1:)
        ENDDO
        polyfile="rim"//TRIM(kstring)//".poly"
        nodefile="rim"//TRIM(kstring)//".1.node"
        elefile="rim"//TRIM(kstring)//".1.ele"
        edgefile="rim"//TRIM(kstring)//".1.edge"
        datfile="rim"//TRIM(kstring)//".dat"
c
c						Get the number of vertex
c						and the number of cells.
        INQUIRE(FILE=TRIM(nodefile),EXIST=file_stat)
        IF (.NOT.file_stat) CALL nim_stop('File '//TRIM(nodefile)//
     $    ' does not exist.  First, run nimset with rimflag=none, '//
     $    'run triangle, rerun fluxgrid, then nimset with '//
     $    'rimflag=tblock.')
        OPEN(UNIT=rim_unit,FILE=TRIM(nodefile),STATUS='OLD')
        READ(rim_unit,*)mvert,idum,idum,idum
        mvert=mvert-1_i4
        CLOSE(UNIT=rim_unit)
        INQUIRE(FILE=TRIM(elefile),EXIST=file_stat)
        IF (.NOT.file_stat) CALL nim_stop('File '//TRIM(elefile)//
     $    ' does not exist.  First, run nimset with rimflag=none, '//
     $    'run triangle, rerun fluxgrid, then nimset with '//
     $    'rimflag=tblock.')
        OPEN(UNIT=rim_unit,FILE=TRIM(elefile),STATUS='OLD')
        READ(rim_unit,*)mcell,idum,idum
        CLOSE(UNIT=rim_unit)
c
c						Allocate contents of new
c						rim tblock
        kb=nblref+ibl
        CALL tri_linear_geom_alloc(tb(kb)%tgeom,mvert,mcell)
        tb(kb)%mvert=mvert
        tb(kb)%mcell=mcell
        CALL tri_linear_alloc(tb(kb)%be_eq,mvert,3_i4)
        tb(kb)%be_eq%name='bq'
        tb(kb)%be_eq%title=' be_eq'
        CALL tri_linear_alloc(tb(kb)%ja_eq,mvert,3_i4)
        tb(kb)%ja_eq%name='jq'
        tb(kb)%ja_eq%title=' ja_eq'
        CALL tri_linear_alloc(tb(kb)%ve_eq,mvert,3_i4)
        tb(kb)%ve_eq%name='vq'
        tb(kb)%ve_eq%title=' ve_eq'
        CALL tri_linear_alloc(tb(kb)%pres_eq,mvert,1_i4)
        tb(kb)%pres_eq%name='pq'
        tb(kb)%pres_eq%title=' pr_eq'
        CALL tri_linear_alloc(tb(kb)%prese_eq,mvert,1_i4)
        tb(kb)%prese_eq%name='pq'
        tb(kb)%prese_eq%title='pre_eq'
        CALL tri_linear_alloc(tb(kb)%nd_eq,mvert,1_i4)
        tb(kb)%nd_eq%name='nq'
        tb(kb)%nd_eq%title=' nd_eq'
        CALL tri_linear_alloc(tb(kb)%diff_shape,mvert,1_i4)
        tb(kb)%diff_shape%name='ds'
        tb(kb)%diff_shape%title='dif_sh'

        ALLOCATE(ibnd(0:mvert))
c
c
c						Get the coordinates of
c						each vertex
c						id,r,z,boundary
c						if ibnd(ivert) = 1
c						it's a boundary cell.
        OPEN(UNIT=rim_unit,FILE=TRIM(nodefile),STATUS='OLD')
        READ(rim_unit,*)idum,idum,idum,idum
        nv_ext=0
        DO i=0,mvert
          READ(rim_unit,*)ivert,
     &                  tb(kb)%tgeom%xs(ivert-1),
     &                  tb(kb)%tgeom%ys(ivert-1),
     &                  ibnd(ivert-1)
          if(ibnd(ivert-1).eq.1_i4)nv_ext=nv_ext+1_i4
        ENDDO
        CLOSE(UNIT=rim_unit)
c
c						Get the verticies for
c						each cell.
        OPEN(UNIT=rim_unit,FILE=TRIM(elefile),STATUS='OLD')
        READ(rim_unit,*)idum,idum,idum
        DO i=1, mcell
          READ(rim_unit,*)icell,iverta,ivertb,ivertc
          tb(kb)%tgeom%vertex(icell,1)=iverta-1_i4
          tb(kb)%tgeom%vertex(icell,2)=ivertb-1_i4
          tb(kb)%tgeom%vertex(icell,3)=ivertc-1_i4
        ENDDO
        CLOSE(UNIT=rim_unit)
c
c						Build the list of neighbors.
c						from the list of line segments.
c							determine vertex
c							bounds on first pass
        INQUIRE(FILE=TRIM(edgefile),EXIST=file_stat)
        IF (.NOT.file_stat) CALL nim_stop('File '//TRIM(edgefile)//
     $    ' does not exist.  First, run nimset with rimflag=none, '//
     $    'run triangle, rerun fluxgrid, then nimset with '//
     $    'rimflag=tblock.')
        OPEN(UNIT=rim_unit,FILE=TRIM(edgefile),STATUS='OLD')
        READ(rim_unit,*)mseg,idum
        ALLOCATE(ibnd1(0:mvert))
        ibnd1(:)=0_i4
        DO i=1,mseg
          READ(rim_unit,*)iseg,ivert1,ivert2,idum
          ivert1=ivert1-1_i4
          ivert2=ivert2-1_i4
          ibnd1(ivert1)=ibnd1(ivert1)+1
          ibnd1(ivert2)=ibnd1(ivert2)+1
        ENDDO
        CLOSE(UNIT=rim_unit)
c
c							Second pass sets
c							vertexes
        OPEN(UNIT=rim_unit,FILE=TRIM(edgefile),STATUS='OLD')
        READ(rim_unit,*)mseg,idum
        DO i=0,mvert
          ALLOCATE(tb(kb)%tgeom%neighbor(i)%vertex(0:ibnd1(i)))
          tb(kb)%tgeom%neighbor(i)%vertex(0)=i
          ibnd1(i)=0_i4
        ENDDO
        DO i=1,mseg
          READ(rim_unit,*)iseg,ivert1,ivert2,idum
          ivert1=ivert1-1_i4
          ivert2=ivert2-1_i4
          ibnd1(ivert1)=ibnd1(ivert1)+1_i4
          ibnd1(ivert2)=ibnd1(ivert2)+1_i4
          tb(kb)%tgeom%neighbor(ivert1)%vertex(ibnd1(ivert1))=ivert2
          tb(kb)%tgeom%neighbor(ivert2)%vertex(ibnd1(ivert2))=ivert1
        ENDDO
        CLOSE(UNIT=rim_unit)
c						
c
c						Read equilibrium magnetic field
c						and pressure at each vertex.
c						These can be computed by
c						reruning fluxgrid, which
c						will notice a rim.node file
c						and compute for the listed
c						r,z values.
        INQUIRE(FILE=TRIM(datfile),EXIST=file_stat)
        IF (.NOT.file_stat) CALL nim_stop('File '//TRIM(datfile)//
     $    ' does not exist.  First, run nimset with rimflag=none, '//
     $    'run triangle, rerun fluxgrid, then nimset with '//
     $    'rimflag=tblock.')
        OPEN(UNIT=rim_unit,FILE=TRIM(datfile),STATUS='OLD')
        tb(kb)%ve_eq=0
        tb(kb)%prese_eq=0
        tb(kb)%ja_eq=0
        tb(kb)%nd_eq=ndens
        tb(kb)%diff_shape=dvac
        DO i=0,mvert
          READ(rim_unit,*)tb(kb)%be_eq%fs(1,i,0),
     &                    tb(kb)%be_eq%fs(2,i,0),
     &                    tb(kb)%be_eq%fs(3,i,0),
     &                    tb(kb)%pres_eq%fs(1,i,0)
          tb(kb)%be_eq%fs(3,i,0)=tb(kb)%be_eq%fs(3,i,0)
     &                          *tb(kb)%tgeom%xs(i)
        ENDDO
        tb(kb)%pres_eq=0
        CLOSE(rim_unit)
        DEALLOCATE(ibnd1,ibnd)
      ENDDO
      RETURN
      END SUBROUTINE poly_rim_init


      SUBROUTINE poly_rim_seam_init(rb,tb,seam,nblref,nrbl)
c-----------------------------------------------------------------------
c     subprogram 4. poly_rim_seam_init.
c	initializes the seams around the rim tblocks.
c			(In retrospect, the seams start and orient
c			around tblocks differently than rblock seams. TAG)
c-----------------------------------------------------------------------

      TYPE(rblock_type), DIMENSION(:), POINTER :: rb
      TYPE(tblock_type), DIMENSION(:), POINTER :: tb
      TYPE(edge_type), DIMENSION(:), POINTER :: seam
      INTEGER(i4),INTENT(in) :: nrbl,nblref

      INTEGER(i4) :: mvert,kb,kv,ip,ib,iv,is,ibl,k,i,idum
      INTEGER(i4) :: ivold,ibold,jbold,jvold,kbold,kvold
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: np
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: jb,jv
      REAL(r8), DIMENSION(:), ALLOCATABLE :: ixv,iyv
      REAL(r8) :: rdum
      CHARACTER(64) :: polyfile,nodefile,elefile,edgefile,datfile
      CHARACTER(8) :: kstring
      LOGICAL :: file_stat
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 2010 FORMAT('rblock(',i2.2,'), lx = ',i3,', ly = ',i3,', mseam = ',i3)
 2020 FORMAT(/4x,'iv',4x,'ip',4x,'ix',4x,'iy',3x,'outb',2x,'outv'/)
 2030 FORMAT(6i6)
 2040 FORMAT('tblock(',i2.2,'), mvert = ',i3,', mseam = ',i3)

c
c                                               Loop through each rim tblock
      DO ibl=1,nbl_rim
c
c                                               Set filenames
        WRITE(kstring,fmt='(i4)')ibl
        DO
          k= INDEX(TRIM(kstring)," ")
          IF(k == 0)EXIT
          kstring(k:)=kstring(k+1:)
        ENDDO
        polyfile="rim"//TRIM(kstring)//".poly"
        nodefile="rim"//TRIM(kstring)//".1.node"
        elefile="rim"//TRIM(kstring)//".1.ele"
        edgefile="rim"//TRIM(kstring)//".1.edge"
        datfile="rim"//TRIM(kstring)//".dat"
c
c						Open the poly file which
c						contains some seam information
c						minus the corners.
        INQUIRE(FILE=TRIM(polyfile),EXIST=file_stat)
        IF (.NOT.file_stat) CALL nim_stop('File '//TRIM(polyfile)//
     $    ' does not exist.  First, run nimset with rimflag=none, '//
     $    'run triangle, rerun fluxgrid, then nimset with '//
     $    'rimflag=tblock.')
        OPEN(UNIT=rim_unit,FILE=TRIM(polyfile),STATUS='OLD')
        READ(rim_unit,*)mvert,idum,idum,idum

        kb=nblref+ibl
        ALLOCATE(seam(kb)%vertex(mvert))
        seam(kb)%nvert=mvert

        DO i=1,mvert
          READ(rim_unit,fmt=*)kv,rdum,rdum,ib,iv
          seam(kb)%vertex(kv)%intxy=(/kv-1_i4,0_i4/)
          ALLOCATE(seam(kb)%vertex(kv)%ptr(2,1))
c						connections on seam being read. 
          seam(kb)%vertex(kv)%ptr(1,1)=ib
          seam(kb)%vertex(kv)%ptr(2,1)=iv
c						connections to an rblock seam
          IF((ib.LE.nblref).AND.(ib.NE.0_i4))THEN
            seam(ib)%vertex(iv)%ptr(1,1)=kb
            seam(ib)%vertex(iv)%ptr(2,1)=kv
          ENDIF
        ENDDO
        CLOSE(rim_unit)
      ENDDO
c
c						Work on the corner seams.
      IF(nbl_rim.EQ.1)THEN
c						Find the first corner
        kb=nblref+1
        DO kv=1,seam(kb)%nvert
          IF(seam(kb)%vertex(kv)%ptr(1,1).NE.0)EXIT
        ENDDO

        ibold=seam(kb)%vertex(kv)%ptr(1,1)
        ivold=seam(kb)%vertex(kv)%ptr(2,1)
        jbold=seam(kb)%vertex(kv+1)%ptr(1,1)
        jvold=seam(kb)%vertex(kv+1)%ptr(2,1)+1

        DEALLOCATE(seam(kb)%vertex(kv)%ptr)
        DEALLOCATE(seam(ibold)%vertex(ivold)%ptr)
        DEALLOCATE(seam(jbold)%vertex(jvold)%ptr)
        ALLOCATE(seam(kb)%vertex(kv)%ptr(2,2))
        ALLOCATE(seam(ibold)%vertex(ivold)%ptr(2,2))
        ALLOCATE(seam(jbold)%vertex(jvold)%ptr(2,2))

        seam(kb)%vertex(kv)%ptr(1,1)=ibold
        seam(kb)%vertex(kv)%ptr(2,1)=ivold
        seam(kb)%vertex(kv)%ptr(1,2)=jbold
        seam(kb)%vertex(kv)%ptr(2,2)=jvold

        seam(ibold)%vertex(ivold)%ptr(1,1)=kb
        seam(ibold)%vertex(ivold)%ptr(2,1)=kv
        seam(ibold)%vertex(ivold)%ptr(1,2)=jbold
        seam(ibold)%vertex(ivold)%ptr(2,2)=jvold

        seam(jbold)%vertex(jvold)%ptr(1,1)=kb
        seam(jbold)%vertex(jvold)%ptr(2,1)=kv
        seam(jbold)%vertex(jvold)%ptr(1,2)=ibold
        seam(jbold)%vertex(jvold)%ptr(2,2)=ivold
c
c						All other corner points are
c						2 rblocks and 1 tblock types.
        DO iv=kv+2,seam(kb)%nvert
          IF(seam(kb)%vertex(iv-1)%ptr(1,1).NE.
     &       seam(kb)%vertex(iv)%ptr(1,1))THEN
c						corner is at jv-1
             ibold=seam(kb)%vertex(iv)%ptr(1,1)
             ivold=seam(kb)%vertex(iv)%ptr(2,1)+1
             jbold=seam(kb)%vertex(iv-1)%ptr(1,1)
             jvold=seam(kb)%vertex(iv-1)%ptr(2,1)

             DEALLOCATE(seam(kb)%vertex(iv-1)%ptr)
             DEALLOCATE(seam(ibold)%vertex(ivold)%ptr)
             DEALLOCATE(seam(jbold)%vertex(jvold)%ptr)
             ALLOCATE(seam(kb)%vertex(iv-1)%ptr(2,2))
             ALLOCATE(seam(ibold)%vertex(ivold)%ptr(2,2))
             ALLOCATE(seam(jbold)%vertex(jvold)%ptr(2,2))
c
             seam(kb)%vertex(iv-1)%ptr(1,1)=ibold
             seam(kb)%vertex(iv-1)%ptr(2,1)=ivold
             seam(kb)%vertex(iv-1)%ptr(1,2)=jbold
             seam(kb)%vertex(iv-1)%ptr(2,2)=jvold

             seam(ibold)%vertex(ivold)%ptr(1,1)=kb
             seam(ibold)%vertex(ivold)%ptr(2,1)=iv-1
             seam(ibold)%vertex(ivold)%ptr(1,2)=jbold
             seam(ibold)%vertex(ivold)%ptr(2,2)=jvold

             seam(jbold)%vertex(jvold)%ptr(1,1)=kb
             seam(jbold)%vertex(jvold)%ptr(2,1)=iv-1
             seam(jbold)%vertex(jvold)%ptr(1,2)=ibold
             seam(jbold)%vertex(jvold)%ptr(2,2)=ivold
          ENDIF
        ENDDO
      ELSE
c						Corner seam between rim tblocks
c						needs to point to seam0 on the 
c						wall.
        DO ibl=1,nbl_rim
          kb=nblref+ibl
c						Search for the corner.
          DO kv=1,seam(kb)%nvert
            IF(seam(kb)%vertex(kv)%ptr(1,1).EQ.0)EXIT
          ENDDO
c						kv points to the wall seam,
          ibold=seam(kb)%vertex(kv)%ptr(1,1)
          ivold=seam(kb)%vertex(kv)%ptr(2,1)-1
c						but kv-1 is the corner and
c						points to a tblock
          kv=kv-1
          jbold=seam(kb)%vertex(kv)%ptr(1,1)
          jvold=seam(kb)%vertex(kv)%ptr(2,1)

          DEALLOCATE(seam(kb)%vertex(kv)%ptr)
          DEALLOCATE(seam(jbold)%vertex(jvold)%ptr)
          ALLOCATE(seam(kb)%vertex(kv)%ptr(2,2))
          ALLOCATE(seam(jbold)%vertex(jvold)%ptr(2,2))

          seam(kb)%vertex(kv)%ptr(1,1)=ibold
          seam(kb)%vertex(kv)%ptr(2,1)=ivold
          seam(kb)%vertex(kv)%ptr(1,2)=jbold
          seam(kb)%vertex(kv)%ptr(2,2)=jvold

          seam(jbold)%vertex(jvold)%ptr(1,1)=ibold
          seam(jbold)%vertex(jvold)%ptr(2,1)=ivold
          seam(jbold)%vertex(jvold)%ptr(1,2)=kb
          seam(jbold)%vertex(jvold)%ptr(2,2)=kv
        ENDDO

c                                               The third corner of the tblock
c						seam points to at least one
c						rblock.
        DO ibl=1,nbl_rim
          kb=nblref+ibl
c                                               Search for the corner.
          DO kv=seam(kb)%nvert-1,1,-1
            IF(seam(kb)%vertex(kv)%ptr(1,1).GT.nrbl)EXIT
          ENDDO
          kv=kv+1
          ibold=seam(kb)%vertex(kv)%ptr(1,1)
          ivold=seam(kb)%vertex(kv)%ptr(2,1)
c						Identify the type of corner.
          IF(ivold.EQ.rb(ibold)%mx)THEN
c							2 rblock & 2 tblocks
            jbold=seam(kb)%vertex(kv-1)%ptr(1,1)
            jvold=seam(jbold)%nvert
c
            kbold=seam(kb)%vertex(kv+1)%ptr(1,1)
            kvold=seam(kb)%vertex(kv+1)%ptr(2,1)+1
           
            DEALLOCATE(seam(kb)%vertex(kv)%ptr)
            DEALLOCATE(seam(kbold)%vertex(kvold)%ptr)
            DEALLOCATE(seam(ibold)%vertex(ivold)%ptr)
            DEALLOCATE(seam(jbold)%vertex(jvold)%ptr)
            ALLOCATE(seam(kb)%vertex(kv)%ptr(2,3))
            ALLOCATE(seam(ibold)%vertex(ivold)%ptr(2,3))
            ALLOCATE(seam(jbold)%vertex(jvold)%ptr(2,3))
            ALLOCATE(seam(kbold)%vertex(kvold)%ptr(2,3))

            seam(kb)%vertex(kv)%ptr(1,1)=ibold
            seam(kb)%vertex(kv)%ptr(2,1)=ivold
            seam(kb)%vertex(kv)%ptr(1,2)=jbold
            seam(kb)%vertex(kv)%ptr(2,2)=jvold
            seam(kb)%vertex(kv)%ptr(1,3)=kbold
            seam(kb)%vertex(kv)%ptr(2,3)=kvold

            seam(jbold)%vertex(jvold)%ptr(1,1)=ibold
            seam(jbold)%vertex(jvold)%ptr(2,1)=ivold
            seam(jbold)%vertex(jvold)%ptr(1,2)=kbold
            seam(jbold)%vertex(jvold)%ptr(2,2)=kvold
            seam(jbold)%vertex(jvold)%ptr(1,3)=kb
            seam(jbold)%vertex(jvold)%ptr(2,3)=kv

            seam(ibold)%vertex(ivold)%ptr(1,1)=kbold
            seam(ibold)%vertex(ivold)%ptr(2,1)=kvold
            seam(ibold)%vertex(ivold)%ptr(1,2)=jbold
            seam(ibold)%vertex(ivold)%ptr(2,2)=jvold
            seam(ibold)%vertex(ivold)%ptr(1,3)=kb
            seam(ibold)%vertex(ivold)%ptr(2,3)=kv

            seam(kbold)%vertex(kvold)%ptr(1,1)=ibold
            seam(kbold)%vertex(kvold)%ptr(2,1)=ivold
            seam(kbold)%vertex(kvold)%ptr(1,2)=jbold
            seam(kbold)%vertex(kvold)%ptr(2,2)=jvold
            seam(kbold)%vertex(kvold)%ptr(1,3)=kb
            seam(kbold)%vertex(kvold)%ptr(2,3)=kv
          ELSE
c							1 rblock & 2 tblock
            jbold=seam(kb)%vertex(kv-1)%ptr(1,1)
            jvold=seam(jbold)%nvert

            DEALLOCATE(seam(kb)%vertex(kv)%ptr)
            DEALLOCATE(seam(ibold)%vertex(ivold)%ptr)
            DEALLOCATE(seam(jbold)%vertex(jvold)%ptr)
            ALLOCATE(seam(kb)%vertex(kv)%ptr(2,2))
            ALLOCATE(seam(ibold)%vertex(ivold)%ptr(2,2))
            ALLOCATE(seam(jbold)%vertex(jvold)%ptr(2,2))

            seam(kb)%vertex(kv)%ptr(1,1)=ibold
            seam(kb)%vertex(kv)%ptr(2,1)=ivold
            seam(kb)%vertex(kv)%ptr(1,2)=jbold
            seam(kb)%vertex(kv)%ptr(2,2)=jvold

            seam(ibold)%vertex(ivold)%ptr(1,1)=kb
            seam(ibold)%vertex(ivold)%ptr(2,1)=kv
            seam(ibold)%vertex(ivold)%ptr(1,2)=jbold
            seam(ibold)%vertex(ivold)%ptr(2,2)=jvold

            seam(jbold)%vertex(jvold)%ptr(1,1)=kb
            seam(jbold)%vertex(jvold)%ptr(2,1)=kv
            seam(jbold)%vertex(jvold)%ptr(1,2)=ibold
            seam(jbold)%vertex(jvold)%ptr(2,2)=ivold
          ENDIF
c						Search for 2 rblock & 1 tblock
c						corners.
          DO iv=kv+2,seam(kb)%nvert-1
            IF(seam(kb)%vertex(iv-1)%ptr(1,1).NE.
     &         seam(kb)%vertex(iv)%ptr(1,1))THEN
c						corner is at jv-1
               ibold=seam(kb)%vertex(iv)%ptr(1,1)
               ivold=seam(kb)%vertex(iv)%ptr(2,1)+1
               jbold=seam(kb)%vertex(iv-1)%ptr(1,1)
               jvold=seam(kb)%vertex(iv-1)%ptr(2,1)

               DEALLOCATE(seam(kb)%vertex(iv-1)%ptr)
               DEALLOCATE(seam(ibold)%vertex(ivold)%ptr)
               DEALLOCATE(seam(jbold)%vertex(jvold)%ptr)
               ALLOCATE(seam(kb)%vertex(iv-1)%ptr(2,2))
               ALLOCATE(seam(ibold)%vertex(ivold)%ptr(2,2))
               ALLOCATE(seam(jbold)%vertex(jvold)%ptr(2,2))

               seam(kb)%vertex(iv-1)%ptr(1,1)=ibold
               seam(kb)%vertex(iv-1)%ptr(2,1)=ivold
               seam(kb)%vertex(iv-1)%ptr(1,2)=jbold
               seam(kb)%vertex(iv-1)%ptr(2,2)=jvold

               seam(ibold)%vertex(ivold)%ptr(1,1)=kb
               seam(ibold)%vertex(ivold)%ptr(2,1)=iv-1
               seam(ibold)%vertex(ivold)%ptr(1,2)=jbold
               seam(ibold)%vertex(ivold)%ptr(2,2)=jvold

               seam(jbold)%vertex(jvold)%ptr(1,1)=kb
               seam(jbold)%vertex(jvold)%ptr(2,1)=iv-1
               seam(jbold)%vertex(jvold)%ptr(1,2)=ibold
               seam(jbold)%vertex(jvold)%ptr(2,2)=ivold
            ENDIF
          ENDDO
        ENDDO
      ENDIF
c
c     						diagnose rblock seams after rim.
      IF(detflag)THEN
         WRITE(out_unit,'(a)')"After pie, After rim."
         DO ib=1,SIZE(rb)
            WRITE(out_unit,2010)ib,rb(ib)%mx,
     &          rb(ib)%my,SIZE(seam(ib)%vertex)
            WRITE(out_unit,2020)
            DO iv=1,SIZE(seam(ib)%vertex)
               DO ip=1,SIZE(seam(ib)%vertex(iv)%ptr,2)
                  WRITE(out_unit,2030)iv,ip,
     &                 seam(ib)%vertex(iv)%intxy,
     &                 seam(ib)%vertex(iv)%ptr(1,ip),
     &                 seam(ib)%vertex(iv)%ptr(2,ip)
               ENDDO
            ENDDO
            WRITE(out_unit,2020)
         ENDDO
c
c     						diagnose tblock seams after rim.
           DO ib=nrbl+1,nblref+nbl_rim
            WRITE(out_unit,2040)ib,tb(ib)%mvert,seam(ib)%nvert
            WRITE(out_unit,2020)
            DO iv=1,seam(ib)%nvert
               DO ip=1,SIZE(seam(ib)%vertex(iv)%ptr,2)
                  WRITE(out_unit,2030)iv,ip,
     &                 seam(ib)%vertex(iv)%intxy,
     &                 seam(ib)%vertex(iv)%ptr(1,ip),
     &                 seam(ib)%vertex(iv)%ptr(2,ip)
               ENDDO
            ENDDO
            WRITE(out_unit,2020)
         ENDDO
      ENDIF
c
      RETURN
      END SUBROUTINE poly_rim_seam_init

c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE poly_init
