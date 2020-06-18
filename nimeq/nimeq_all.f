c-----------------------------------------------------------------------
c     file nimeq_all.f
c     this module holds block and seam structures that are used for
c     assembling the full domain when running in parallel.  they allow
c     field-line tracing without communication during each trace
c     although multiple copies of data are kept in memory.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  nimeq_all
c     1.  nimeq_all_init
c     2.  nimeq_all_dealloc
c     3.  nimeq_all_update
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     0. module declaration for nimeq_btr
c-----------------------------------------------------------------------
      MODULE nimeq_all
      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE edge_type_mod
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     blocks for collecting data and for retaining pointers.
c-----------------------------------------------------------------------
      TYPE(rblock_type), DIMENSION(:), POINTER :: rball,rbtmp
      TYPE(tblock_type), DIMENSION(:), POINTER :: tball,tbtmp
c-----------------------------------------------------------------------
c     seam structures for collecting and for pointers.
c-----------------------------------------------------------------------
      TYPE(edge_type), DIMENSION(:), POINTER :: seamall,seamtmp
c-----------------------------------------------------------------------
c     arrays are needed for the allgather calls.
c     nalloff(0:nprocs) is the node count offset and
c     nallcount(0:nprocs-1) is the number of nodes for each process.
c     rlocal and rglobal are buffers for the send and receive data.
c-----------------------------------------------------------------------
      INTEGER(i4), DIMENSION(:), ALLOCATABLE, PRIVATE :: nalloff,
     $  nallcount
      REAL(r8), DIMENSION(:), ALLOCATABLE, PRIVATE :: rlocal,rglobal

      CONTAINS
c-----------------------------------------------------------------------
c     1. subprogram nimeq_all_init
c     this subroutine allocates and global rb, tb, and seam structures
c     and loads them with the information for tracing that does not
c     change from one step to the next.
c-----------------------------------------------------------------------
      SUBROUTINE nimeq_all_init(poly_deg)
      USE mpi_nim
      USE pardata
      USE fields
      USE seam_storage_mod

      INTEGER(i4), INTENT(IN) :: poly_deg

      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: segpt
      INTEGER(i4) :: ibl,ncnt,ierr,ipr,ibgl,iv,icnt,ncpy
c-----------------------------------------------------------------------
c     start with allocations.
c-----------------------------------------------------------------------
      ALLOCATE(rball(nrbl_total))
      ALLOCATE(tball(nbl_total-nrbl_total))
      ALLOCATE(seamall(nbl_total))
c-----------------------------------------------------------------------
c     count the number of finite element nodes in each block on the
c     current process 
c-----------------------------------------------------------------------
      ALLOCATE(nalloff(0:nprocs),nallcount(0:nprocs-1))
      ncnt=0
      DO ibl=1,nrbl
        ncnt=ncnt+SIZE(rb(ibl)%rz%fs)
        IF (ALLOCATED(rb(ibl)%rz%fsh))
     $    ncnt=ncnt+SIZE(rb(ibl)%rz%fsh)
        IF (ALLOCATED(rb(ibl)%rz%fsv))
     $    ncnt=ncnt+SIZE(rb(ibl)%rz%fsv)
        IF (ALLOCATED(rb(ibl)%rz%fsi))
     $    ncnt=ncnt+SIZE(rb(ibl)%rz%fsi)
      ENDDO
      DO ibl=nrbl+1,nbl
        ncnt=ncnt+2_i4*(tb(ibl)%mvert+1_i4)
      ENDDO
      CALL mpi_allgather(ncnt,1,mpi_nim_int,nallcount(0),1,
     $                   mpi_nim_int,mpi_comm_world,ierr)
      nalloff(0)=0
      DO ipr=1,nprocs
        nalloff(ipr)=nalloff(ipr-1)+nallcount(ipr-1)
      ENDDO
c-----------------------------------------------------------------------
c     allocate structures using previously communicated information.
c-----------------------------------------------------------------------
      DO ibgl=1,nrbl_total
        rball(ibgl)%id=ibgl
        rball(ibgl)%mx=block_sizes(2,ibgl)
        rball(ibgl)%my=block_sizes(3,ibgl)
        CALL lagr_quad_2D_alloc(rball(ibgl)%rz,rball(ibgl)%mx,
     $                          rball(ibgl)%my,2_i4,poly_deg)
        CALL lagr_quad_2D_alloc(rball(ibgl)%be_n0,rball(ibgl)%mx,
     $                          rball(ibgl)%my,4_i4,poly_deg)
        ipr=block2proc(ibgl)
        ibl=global2local(ibgl)
        IF (ipr==node) THEN
          rball(ibgl)%degenerate=rb(ibl)%degenerate
          rball(ibgl)%periodicy=rb(ibl)%periodicy
        ENDIF
        CALL mpi_bcast(rball(ibgl)%degenerate,1,mpi_nim_logical,ipr,
     $       mpi_comm_world,ierr)
        CALL mpi_bcast(rball(ibgl)%periodicy,1,mpi_nim_logical,ipr,
     $       mpi_comm_world,ierr)
      ENDDO

      DO ibgl=nrbl_total+1,nbl_total
        tball(ibgl)%id=ibgl
        tball(ibgl)%mvert=block_sizes(2,ibgl)
        tball(ibgl)%mcell=block_sizes(3,ibgl)
        ALLOCATE(tball(ibgl)%tgeom%xs(0:tball(ibgl)%mvert))
        ALLOCATE(tball(ibgl)%tgeom%ys(0:tball(ibgl)%mvert))
        CALL tri_linear_alloc(tball(ibgl)%be_n0,tball(ibgl)%mvert,4_i4)
      ENDDO
c-----------------------------------------------------------------------
c     communicate the node coordinates.  the local process loads its
c     blocks in the order set by the global data structures.
c-----------------------------------------------------------------------
      ALLOCATE(rlocal(nallcount(node)),rglobal(nalloff(nprocs)))
      ncnt=0
      DO ibgl=1,nrbl_total
        IF (block2proc(ibgl)/=node) CYCLE
        ibl=global2local(ibgl)
        icnt=ncnt+1
        ncpy=SIZE(rb(ibl)%rz%fs)
        ncnt=ncnt+ncpy
        CALL real_copy(rlocal(icnt),rb(ibl)%rz%fs(1,0,0),ncpy)
        IF (ALLOCATED(rb(ibl)%rz%fsh)) THEN
          icnt=ncnt+1
          ncpy=SIZE(rb(ibl)%rz%fsh)
          ncnt=ncnt+ncpy
          CALL real_copy(rlocal(icnt),rb(ibl)%rz%fsh(1,1,1,0),ncpy)
        ENDIF
        IF (ALLOCATED(rb(ibl)%rz%fsv)) THEN
          icnt=ncnt+1
          ncpy=SIZE(rb(ibl)%rz%fsv)
          ncnt=ncnt+ncpy
          CALL real_copy(rlocal(icnt),rb(ibl)%rz%fsv(1,1,0,1),ncpy)
        ENDIF
        IF (ALLOCATED(rb(ibl)%rz%fsi)) THEN
          icnt=ncnt+1
          ncpy=SIZE(rb(ibl)%rz%fsi)
          ncnt=ncnt+ncpy
          CALL real_copy(rlocal(icnt),rb(ibl)%rz%fsi(1,1,1,1),ncpy)
        ENDIF
      ENDDO
      DO ibgl=nrbl_total+1,nbl_total
        IF (block2proc(ibgl)/=node) CYCLE
        ibl=global2local(ibgl)
        icnt=ncnt+1
        ncnt=ncnt+tb(ibl)%mvert+1
        rlocal(icnt:ncnt)=tb(ibl)%tgeom%xs
        icnt=ncnt+1
        ncnt=ncnt+tb(ibl)%mvert+1
        rlocal(icnt:ncnt)=tb(ibl)%tgeom%ys
      ENDDO
c-----------------------------------------------------------------------
c     collect information with a gather communication operation.
c-----------------------------------------------------------------------
      CALL mpi_allgatherv(rlocal(1),nallcount(node),mpi_nim_real,
     $     rglobal(1),nallcount(0),nalloff(0),mpi_nim_real,
     $     mpi_comm_world,ierr)
c-----------------------------------------------------------------------
c     extract the information in the blocks, paying attention to the
c     order of the packing in rglobal.
c-----------------------------------------------------------------------
      DO ipr=0,nprocs-1
        ncnt=nalloff(ipr)
        DO ibgl=1,nrbl_total
          IF (block2proc(ibgl)==ipr) THEN
            icnt=ncnt+1
            ncpy=SIZE(rball(ibgl)%rz%fs)
            ncnt=ncnt+ncpy
            CALL real_copy(rball(ibgl)%rz%fs(1,0,0),rglobal(icnt),ncpy)
            IF (ALLOCATED(rball(ibgl)%rz%fsh)) THEN
              icnt=ncnt+1
              ncpy=SIZE(rball(ibgl)%rz%fsh)
              ncnt=ncnt+ncpy
              CALL real_copy(rball(ibgl)%rz%fsh(1,1,1,0),
     $                       rglobal(icnt),ncpy)
            ENDIF
            IF (ALLOCATED(rball(ibgl)%rz%fsv)) THEN
              icnt=ncnt+1
              ncpy=SIZE(rball(ibgl)%rz%fsv)
              ncnt=ncnt+ncpy
              CALL real_copy(rball(ibgl)%rz%fsv(1,1,0,1),
     $                       rglobal(icnt),ncpy)
            ENDIF
            IF (ALLOCATED(rball(ibgl)%rz%fsi)) THEN
              icnt=ncnt+1
              ncpy=SIZE(rball(ibgl)%rz%fsi)
              ncnt=ncnt+ncpy
              CALL real_copy(rball(ibgl)%rz%fsi(1,1,1,1),
     $                       rglobal(icnt),ncpy)
            ENDIF
          ENDIF
        ENDDO

        DO ibgl=nrbl_total+1,nbl_total
          IF (block2proc(ibgl)==ipr) THEN
            icnt=ncnt+1
            ncpy=SIZE(tball(ibgl)%tgeom%xs)
            ncnt=ncnt+ncpy
            CALL real_copy(tball(ibgl)%tgeom%xs(0),rglobal(icnt),ncpy)
            icnt=ncnt+1
            ncpy=SIZE(tball(ibgl)%tgeom%ys)
            ncnt=ncnt+ncpy
            CALL real_copy(tball(ibgl)%tgeom%ys(0),rglobal(icnt),ncpy)
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     reset communication arrays for subsequent updates.
c-----------------------------------------------------------------------
      DEALLOCATE(rlocal,rglobal)
      nalloff=2_i4*nalloff
      nallcount=2_i4*nallcount
      ALLOCATE(rlocal(nallcount(node)),rglobal(nalloff(nprocs)))
c-----------------------------------------------------------------------
c     traces only use nvert, segment%ptr and segment%intxys from the
c     seam data structures.
c-----------------------------------------------------------------------
      DO ibgl=1,nbl_total
        seamall(ibgl)%nvert=block_sizes(1,ibgl)
        ALLOCATE(seamall(ibgl)%segment(seamall(ibgl)%nvert))
        ALLOCATE(segpt(4,seamall(ibgl)%nvert))
        ipr=block2proc(ibgl)
        IF (ipr==node) THEN
          ibl=global2local(ibgl)
          DO iv=1,seam(ibl)%nvert
            segpt(1:2,iv)=seam(ibl)%segment(iv)%ptr
            segpt(3:4,iv)=seam(ibl)%segment(iv)%intxys
          ENDDO
        ENDIF
        CALL mpi_bcast(segpt(1,1),4_i4*seamall(ibgl)%nvert,mpi_nim_int,
     $                 ipr,mpi_comm_world,ierr)
        DO iv=1,seamall(ibgl)%nvert
          seamall(ibgl)%segment(iv)%ptr=segpt(1:2,iv)
          seamall(ibgl)%segment(iv)%intxys=segpt(3:4,iv)
        ENDDO
        DEALLOCATE(segpt)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE nimeq_all_init
c-----------------------------------------------------------------------
c     2. subprogram nimeq_all_dealloc
c     deallocate global rb, tb, and seam structures.
c-----------------------------------------------------------------------
      SUBROUTINE nimeq_all_dealloc
      USE mpi_nim
      USE pardata
      USE fields

      INTEGER(i4) :: ibgl

      DO ibgl=1,nbl_total
        IF (ibgl<=nrbl_total) THEN
          CALL lagr_quad_dealloc(rball(ibgl)%rz)
          CALL lagr_quad_dealloc(rball(ibgl)%be_n0)
        ELSE
          DEALLOCATE(tball(ibgl)%tgeom%xs)
          DEALLOCATE(tball(ibgl)%tgeom%ys)
          CALL tri_linear_dealloc(tball(ibgl)%be_n0)
        ENDIF
        DEALLOCATE(seamall(ibgl)%segment)
      ENDDO

      DEALLOCATE(rball,tball,seamall)
      DEALLOCATE(nalloff,nallcount)
      DEALLOCATE(rglobal,rlocal)

      RETURN
      END SUBROUTINE nimeq_all_dealloc
c-----------------------------------------------------------------------
c     3. subprogram nimeq_all_update
c     this subroutine updates the flux function in the global blocks
c     for use by the field-line tracing.  the field is in the be_n0
c     data structure.
c-----------------------------------------------------------------------
      SUBROUTINE nimeq_all_update
      USE mpi_nim
      USE pardata
      USE fields

      INTEGER(i4) :: ibl,ncnt,ierr,ipr,ibgl,icnt,ncpy
c-----------------------------------------------------------------------
c     update the be_n0 data structure.  the local process loads its
c     blocks in the order set by the global data structures.
c-----------------------------------------------------------------------
      ncnt=0
      DO ibgl=1,nrbl_total
        IF (block2proc(ibgl)/=node) CYCLE
        ibl=global2local(ibgl)
        icnt=ncnt+1
        ncpy=SIZE(rb(ibl)%be_n0%fs)
        ncnt=ncnt+ncpy
        CALL real_copy(rlocal(icnt),rb(ibl)%be_n0%fs(1,0,0),ncpy)
        IF (ALLOCATED(rb(ibl)%be_n0%fsh)) THEN
          icnt=ncnt+1
          ncpy=SIZE(rb(ibl)%be_n0%fsh)
          ncnt=ncnt+ncpy
          CALL real_copy(rlocal(icnt),rb(ibl)%be_n0%fsh(1,1,1,0),ncpy)
        ENDIF
        IF (ALLOCATED(rb(ibl)%be_n0%fsv)) THEN
          icnt=ncnt+1
          ncpy=SIZE(rb(ibl)%be_n0%fsv)
          ncnt=ncnt+ncpy
          CALL real_copy(rlocal(icnt),rb(ibl)%be_n0%fsv(1,1,0,1),ncpy)
        ENDIF
        IF (ALLOCATED(rb(ibl)%be_n0%fsi)) THEN
          icnt=ncnt+1
          ncpy=SIZE(rb(ibl)%be_n0%fsi)
          ncnt=ncnt+ncpy
          CALL real_copy(rlocal(icnt),rb(ibl)%be_n0%fsi(1,1,1,1),ncpy)
        ENDIF
      ENDDO
      DO ibgl=nrbl_total+1,nbl_total
        IF (block2proc(ibgl)/=node) CYCLE
        ibl=global2local(ibgl)
        icnt=ncnt+1
        ncpy=4_i4*(tb(ibl)%mvert+1)
        ncnt=ncnt+ncpy
        CALL real_copy(rlocal(icnt),tb(ibl)%be_n0%fs(1,0,0),ncpy)
      ENDDO
c-----------------------------------------------------------------------
c     collect information with a gather communication operation.
c-----------------------------------------------------------------------
      CALL mpi_allgatherv(rlocal(1),nallcount(node),mpi_nim_real,
     $     rglobal(1),nallcount(0),nalloff(0),mpi_nim_real,
     $     mpi_comm_world,ierr)
c-----------------------------------------------------------------------
c     extract the information in the blocks, paying attention to the
c     order of the packing in rglobal.
c-----------------------------------------------------------------------
      DO ipr=0,nprocs-1
        ncnt=nalloff(ipr)
        DO ibgl=1,nrbl_total
          IF (block2proc(ibgl)==ipr) THEN
            icnt=ncnt+1
            ncpy=SIZE(rball(ibgl)%be_n0%fs)
            ncnt=ncnt+ncpy
            CALL real_copy(rball(ibgl)%be_n0%fs(1,0,0),
     $                     rglobal(icnt),ncpy)
            IF (ALLOCATED(rball(ibgl)%be_n0%fsh)) THEN
              icnt=ncnt+1
              ncpy=SIZE(rball(ibgl)%be_n0%fsh)
              ncnt=ncnt+ncpy
              CALL real_copy(rball(ibgl)%be_n0%fsh(1,1,1,0),
     $                       rglobal(icnt),ncpy)
            ENDIF
            IF (ALLOCATED(rball(ibgl)%be_n0%fsv)) THEN
              icnt=ncnt+1
              ncpy=SIZE(rball(ibgl)%be_n0%fsv)
              ncnt=ncnt+ncpy
              CALL real_copy(rball(ibgl)%be_n0%fsv(1,1,0,1),
     $                       rglobal(icnt),ncpy)
            ENDIF
            IF (ALLOCATED(rball(ibgl)%be_n0%fsi)) THEN
              icnt=ncnt+1
              ncpy=SIZE(rball(ibgl)%be_n0%fsi)
              ncnt=ncnt+ncpy
              CALL real_copy(rball(ibgl)%be_n0%fsi(1,1,1,1),
     $                       rglobal(icnt),ncpy)
            ENDIF
          ENDIF
        ENDDO

        DO ibgl=nrbl_total+1,nbl_total
          IF (block2proc(ibgl)==ipr) THEN
            icnt=ncnt+1
            ncpy=SIZE(tball(ibgl)%be_n0%fs)
            ncnt=ncnt+ncpy
            CALL real_copy(tball(ibgl)%be_n0%fs(1,0,0),
     $                     rglobal(icnt),ncpy)
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE nimeq_all_update
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE nimeq_all

c-----------------------------------------------------------------------
c     subprogram 4. real_copy
c     the following external subroutine avoids reshape operations.
c-----------------------------------------------------------------------
      SUBROUTINE real_copy(rout,rin,ncount)
      USE local
      IMPLICIT NONE

      REAL(r8), DIMENSION(*), INTENT(OUT) :: rout
      REAL(r8), DIMENSION(*), INTENT(IN) :: rin
      INTEGER(i4), INTENT(IN) :: ncount

      rout(1:ncount)=rin(1:ncount)

      RETURN
      END SUBROUTINE real_copy
